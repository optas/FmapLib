classdef Patch < dynamicprops
    
    properties (GetAccess = public, SetAccess = private)
        % Basic properties that every instance of the Patch class has.        
        source_image    %   (Image)           -    Image over which the patch was sampled.
        
        corners         %   (1 x 4 uint)    -    x-y coordinates wrt. source image, corresponding to the corners 
                        %                          of the rectangular patch. They are [xmin, ymin, xmax, ymax]
    end
       
    methods (Access = public)
        % Class constructror
        function obj = Patch(src_img, corners)
            if nargin == 0
                obj.source_image = Image();
                obj.corners      = zeros(1,4);
            else
                if ~ Patch.valid_corners(corners, src_image)
                    error('The corners of the Patch do not follow the xmin, ymin, xmax, ymax protocol.')
                end
                obj.source_image = src_img;                
                obj.corners      = corners;
                
            end
        end
        
        function [xmin, ymin, xmax, ymax] = get_corners(obj)
            % Getter of object's property 'corners' corresponding to the 4 extema of
            % the x-y coordinates of the patch.
            xmin = obj.corners(1);
            ymin = obj.corners(2);
            xmax = obj.corners(3);
            ymax = obj.corners(4);
        end
    
        function [F] = plot(obj)
         % Plots the boundary of the patch on its source image.
            shape_inserter = vision.ShapeInserter('LineWidth', 4);         
            [xmin, ymin, xmax, ymax] = obj.get_corners();
            rectangle = int32([xmin ymin (xmax - xmax) (ymax - ymin)]);
            im_out    = step(shape_inserter, obj.source_image, rectangle);            
            image(im_out);
        end

        function a = area(obj)
                [xmin, ymin, xmax, ymax] = obj.get_corners();
                a = (ymax-ymin) * (xmax - xmin);
        end
        
        function area = area_of_intersection(obj, another_patch)
            
            [xmin1, ymin1, xmax1, ymax1] = obj.get_corners();
            [xmin2, ymin2, xmax2, ymax2] = another_patch.get_corners();
            
            xmin = max(xmin1, xmin2);
            ymin = max(ymin1, ymin2);
            xmax = min(xmax1, xmax2);
            ymax = min(ymax1, ymax2);
            
            area = max(0, (ymax-ymin) * (xmax-xmin));
        end
            
            
    end
    
    methods (Static, Access = public)
        function [b] = is_valid(corners, src_image)            
            b =  corners(1) <= src_image.width  && corners(3) <= src_image.width && ...   % Inside photo's x-dim.
                 corners(2) <= src_image.height && corners(4) <= src_image.height && ...  % Inside photo's y-dim.
                 all(corners) > 0 ;                                         
        end
        
        function [b] = is_within_limits(patch,  width_limit, height_limit, src_image)
            if width_limit <= 0  || height_limit <= 0 || height_limit > 1 || width_limit > 1
                error('Limits must be percents in (0,1] and they refer to the size of the original image where the patch comes from.')
            end
                
            b  = (patch(3) - patch(1))  <= width_limit  * src_image.width   && ...
                 (patch(4) - patch(2))  <= height_limit * src_image.height ;
            
        end
            
        
        
        function [F] = extract_patch_features(corners, features, type)
                switch type
                    case 'zero_pad'
                        F = zeros(size(features));
                        F(corners(2):corners(4), corners(1):corners(3), :) = features(corners(2):corners(4), corners(1):corners(3), :);
                    case 'tiling'
                        tile = features(corners(2):corners(4), corners(1):corners(3), :);
                        tile_length = ceil(size(features,2)/size(tile,2));
                        tile_height = ceil(size(features,1)/size(tile,1));
                        
                        big_tiled = repmat(tile, tile_height, tile_length);
                        F = big_tiled(1:size(features, 1), 1:size(features, 2),:);                 
                    otherwise
                        error('Type not implemented.')
                end
        end
                
        function new_corners = find_new_corners(old_height, old_width, new_height, new_width, old_corners)
            xmin = old_corners(1);
            ymin = old_corners(2);
            xmax = old_corners(3);
            ymax = old_corners(4);
            
            xmin_n = (double(xmin) * new_width)  / old_width;                % Linear interpolation case.
            xmax_n = (double(xmax) * new_width)  / old_width;
            ymin_n = (double(ymin) * new_height) / old_height;
            ymax_n = (double(ymax) * new_height) / old_height;
            
            new_corners = [xmin_n ymin_n xmax_n ymax_n];
            new_corners = uint16(max(1, round(new_corners)));            
        end
        
        function [top_patches, top_scores] = compute_top_patches_in_images(nns, patch_feat, fmaps, top_p)
            num_images = size(patch_feat, 1);
            top_patches = zeros(num_images, top_p);
            top_scores  = zeros(num_images, top_p);            
            for i = 1:num_images
                loc_i = find(~cellfun(@isempty, patch_feat(i,:)));
                score_i = zeros(length(loc_i), 1);
                               
                for pi = loc_i % Proposals for image_i       
                    for j = nns(i, :)
                        mis_pi_j = +Inf;           
                        loc_j = find(~cellfun(@isempty, patch_feat(j,:)));
                        for pj = loc_j      %Proposals for image_j                        
                            mis_pi_j = min(mis_pi_j, patch_transferability(i, j, patch_feat{i, pi}, patch_feat{j, pj}, fmaps));
                        end
%                         mis_pi_j = mis_pi_j / length(loc_j);
                    end
                    score_i(pi) = score_i(pi)  + mis_pi_j;
                end    

                [sort_scores, indices] = sort(score_i);
                top_patches(i, :)      = indices(1:top_p);
                top_scores(i,:)        = sort_scores(1:top_p);                
            end
        end
                
        function [new_nns, top_patches] = iterative_compute_patches(nns, patch_feat, fmaps, top_p, number_iters)                        
            num_images = size(patch_feat, 1);            
            new_nns     = zeros(num_images, top_p, number_iters + 1);
            new_nns(:, :, 1) = nns;   
            top_patches = zeros(num_images, top_p, number_iters);
                                                
            for k = 1:number_iters                
                [top_indices, top_scores] = Patch.compute_top_patches_in_images(new_nns(:,:,k), patch_feat, fmaps, top_p);                
                top_patches(:,:,k) = top_indices;
                top_feat = cell(num_images, top_p);
            
                for i = 1:num_images                    
                    top_feat(i, :) = patch_feat(i, top_indices(i, :));                     
                end
            
                new_distances = Patch.all_pairwise_distances(top_feat, fmaps);
                [sorted_dists, sorted_indices] = sort(new_distances);
                new_nns(:,:,k+1) = sorted_indices(2 : top_p + 1, :)' ;            
            end
            
        end
                        
        function distances = all_pairwise_distances(top_patches, all_fmaps)
            distances = zeros(length(top_patches));
    
            for i= 1:length(top_patches) -1
                for j=i+1:length(top_patches)
                    distances(i,j) = Patch.distance_between_two_sets_of_patches(i,j,top_patches, all_fmaps);
                    distances(j,i) = distances(i,j);
                end
            end
        end
        
        function distance = distance_between_two_sets_of_patches(img1, img2, top_patches, all_fmaps)
            distances = zeros(size(top_patches, 1), size(top_patches, 2));
            for i = 1:length(top_patches(img1, :))
                for j = 1:length(top_patches(img2, :))
                    dist = patch_transferability(img1, img2, top_patches{img1, i}, top_patches{img2, j}, all_fmaps);
                    distances(i,j) = dist;
                end                
            end                
            
            distance = sum(min(distances))
        end
        
%         function [misalignment] = patch_transferability_triplet(source_image, target_image, source_patch, target_patch, fmaps) 
%             misalignment = sum(sum(abs((fmaps{source_image, target_image} * source_patch) - target_patch), 1));
%             misalignment = misalignment + sum(sum(abs((fmaps{target_image, source_image} * target_patch) - source_patch), 1));
% 
%         end
        

         
    end
    
end