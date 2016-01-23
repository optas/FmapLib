classdef Patch < dynamicprops
    
    properties (GetAccess = public, SetAccess = private)
        % Basic properties that every instance of the Patch class has.                
        corners         %   (1 x 4 uint)      -    x-y coordinates wrt. source image, corresponding to the corners 
                        %                          of the rectangular patch. They are [xmin, ymin, xmax, ymax].
    end
       
    methods (Access = public)                
        function obj = Patch(varargin)
            % Class constructror:
            %   Input: 
            %             (1 x 4) vector describing the corners as: [xmin, ymin, xmax, ymax].
            if nargin == 0                
                obj.corners  = zeros(1,4);
            elseif nargin == 1
                in_corners = varargin{1};
                if ~ Patch.is_valid(in_corners)
                    error('The corners of the Patch do not follow the xmin, ymin, xmax, ymax protocol.')
                else                    
                    obj.corners = in_corners;
                end            
            else
                error('Wrong number of arguments.')
            end
        end
           
        function [varargout] = size(obj)                        
            if length(obj) > 1                % Array of Patches.
                varargout{1} = length(obj);               
                return
            end
            if nargout == 2
                varargout = cell(nargout);
                varargout{1} = obj.width;
                varargout{2} = obj.height;
            else
                varargout{1} = [obj.width, obj.height];
            end
        end
        
        function [varargout] = get_corners(obj)
            % Getter of object's property 'corners'.
            if nargout == 4
                varargout = cell(nargout);
                varargout{1} = obj.corners(1);
                varargout{2} = obj.corners(2);
                varargout{3} = obj.corners(3);
                varargout{4} = obj.corners(4);
            else
                varargout{1} = obj.corners;
            end
        end
            
        function [F] = plot(obj, varargin)
            % Plots the boundary of the patch.
            options = struct('image', [], 'color', 'r', 'line_width', 3);
            options = load_key_value_input_pairs(options, varargin{:});         
            
            [xmin, ymin, xmax, ymax] = obj.get_corners();
            if ~ isempty(options.image)
                image.plot();            
                hold on;
            end
            F = plot([xmin xmax xmax xmin xmin],[ymin ymin ymax ymax ymin], 'Color', options.color, ...
                                                                            'LineWidth', options.line_width);
        end

        function w = width(obj)                       
            [~, ymin, ~, ymax] = obj.get_corners();
            w = ymax - ymin + 1;
            w = double(w);
            assert(w >= 1);
        end
        
        function h = height(obj)                       
            [xmin, ~, xmax, ~] = obj.get_corners();
            h = xmax - xmin + 1;
            h = double(h);
            assert(h >= 1);
        end
                
        function a = area(obj)                
            a = obj.height() * obj.width();            
        end
        
        function a = area_of_intersection(obj, another_patch)            
            [xmin1, ymin1, xmax1, ymax1] = obj.get_corners();
            [xmin2, ymin2, xmax2, ymax2] = another_patch.get_corners();            
            xmin = max(xmin1, xmin2);
            ymin = max(ymin1, ymin2);
            xmax = min(xmax1, xmax2);
            ymax = min(ymax1, ymax2);
            if ymin <= ymax && xmin <= xmax
                a = double((ymax-ymin+1)) * double((xmax-xmin+1));
            else
                a = 0;
            end
        end
        
        function c = corloc(obj, another_patch)
            intersection = obj.area_of_intersection(another_patch);
            union        = obj.area() + another_patch.area() - intersection;
            c            = double(intersection) / double(union);
            assert(c >= 0 && c <= 1);
        end
        
        function [np] = linear_interpolation(obj, from_width, to_width, from_height, to_height)                                   
            [xmin, ymin, xmax, ymax] = obj.get_corners();                        
            n_xmin = (double(xmin) * to_width)  / from_width;
            n_xmax = (double(xmax) * to_width)  / from_width;
            n_ymin = (double(ymin) * to_height) / from_height;
            n_ymax = (double(ymax) * to_height) / from_height;                                
            new_corners = [n_xmin, n_ymin, n_xmax, n_ymax];
            new_corners = uint16(max(1, round(new_corners)));
            assert(Patch.is_valid(new_corners, Image(zeros(to_height, to_width))));
            np = Patch(new_corners);                       
        end
        
        function [F] = extract_features(obj, features, method)
                in_corners = obj.corners();
                switch method
                    case 'mask'
                        F = features(in_corners(2):in_corners(4), in_corners(1):in_corners(3), :);                    
                    case 'zero_pad'
                        F = zeros(size(features));
                        F(in_corners(2):in_corners(4), in_corners(1):in_corners(3), :) = features(in_corners(2):in_corners(4), in_corners(1):in_corners(3), :);
                    case 'tiling'
                        tile = features(in_corners(2):in_corners(4), in_corners(1):in_corners(3), :);
                        tile_length = ceil(size(features,2)/size(tile,2));
                        tile_height = ceil(size(features,1)/size(tile,1));
                        
                        big_tiled = repmat(tile, tile_height, tile_length);
                        F = big_tiled(1:size(features, 1), 1:size(features, 2),:);                                     
                    otherwise
                        error('Type not implemented.')
                end
        end
                       
    end
    
    methods (Static, Access = public)
        function [b] = is_valid(corners, src_image)
            % Checks if the corners obey the necessary conditions for being a patch.
            if exist('src_image', 'var') % Verify it fits the image.
                b =  corners(1) <= src_image.width  && corners(3) <= src_image.width && ...   
                     corners(2) <= src_image.height && corners(4) <= src_image.height ;
            else
                b = true;                               
            end            
            b = b && all(corners > 0) ...
                  && length(corners) == 4 && corners(1) <= corners(3)  && corners(2) <= corners(4);                  
        end
        
        function [P] = tightest_box_of_segment(in_mask)            
            [y, x] = find(in_mask);
            P      = Patch([min(x), min(y), max(x), max(y)]);            
            
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
                        for pj = loc_j      % Proposals for image_j                        
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
        
        
        function [b] = uodl_patch_constraint(corners, src_image)                        
            bValid1 = corners(1) > src_image.height *0.01 & corners(3) < src_image.height*0.99 ...
                    & corners(2) > src_image.width * 0.01 & corners(4) < src_image.width*0.99;
            bValid2 = corners(1) < src_image.height*0.01 & corners(3) > src_image.height*0.99 ...
                    & corners(2) < src_image.width*0.01 & corners(4) > src_image.width*0.99;
            b = bValid1 | bValid2;                           
        end 

%    function [misalignment] = patch_transferability(source_image, target_image, source_patch, target_patch, fmaps)
%     
%         misalignment = sum(sum(abs((fmaps{P, target_image} * source_patch) - target_patch), 1));
%         misalignment = misalignment + sum(sum(abs((fmaps{target_image, P} * target_patch) - source_patch), 1));
% 
%     end     


         
    end
    
end