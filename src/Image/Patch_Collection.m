classdef Patch_Collection < dynamicprops
    % TODO-s aggregate: implement subref
    
    properties (GetAccess = public, SetAccess = private)
        collection = Patch;        % (Patch Array) storing all the patches.
        image;
    end

    methods (Access = public)                
        function obj = Patch_Collection(corners, in_image)            
            % Class Constructor.
            % Input:
            %           corners  - (N x 4) matrix carrying the (x,y) coordinates of the 2 extreme corners defining a
            %                      rectangular patch. The expected format for each row is: [xmin, ymin, xmax, ymax].                
            %           in_image - (Image) The image, over which the patches were defined.
            if nargin == 0
                obj.collection = Patch();            
                obj.image      = Image();
            else
                obj.image = in_image;                
                if isa(corners, 'double')                
                    m = size(corners, 1);    
                    if m < 1 
                        error('A Patch Collection must consist of at least one Patch objects.')
                    end                
                    obj.collection(1:m) = Patch();
                    for i = 1:m
                        if ~ Patch.is_valid(corners(i,:), in_image)
                            error('Patch %d cannot go with given image.\n', i);
                            
                        else
                            obj.collection(i) = Patch(corners(i, :));
                        end                    
                    end
                else
                    error('Not correct argument.')
                end                
            end
        end
        
        function m = size(obj)
            if length(obj) > 1              % An array of Patch_Collection objects.
                m = length(obj);                
            else
                m = length(obj.collection);
            end
        end
        
        function p = get_patch(obj, p_index)
            p = obj.collection(p_index) ;
        end

        function F = plot_collection(obj, scores)
            obj.image.plot();
            hold on;            
            if exist('scores', 'var')                
                colors = vals2colormap(scores);
                for i = 1:length(scores)
                    if scores(i) ~= 0
                        F = obj.get_patch(i).plot('color', colors(i,:));
                    end
                end
            else                    
                for i = 1:size(obj)
                    F = obj.get_patch(i).plot();
                end                
            end
        end
                
        function I = masked_image(obj)
            I = obj.image;            
            mask = zeros(I.height, I.width);
            for i = 1:size(obj)
                pi = obj.get_patch(i);
                [xmin, ymin, xmax, ymax] = pi.get_corners;
                mask(ymin:ymax, xmin:xmax) = 1;
            end 
            I = I.apply_mask(mask);            
        end
        
        function [keep] = filter_patches(obj, varargin)
            % Keep patches that are not area-wise too big or too small wrt. to the obj.image.
            % Also removes line patches.
            
            options = struct('no_lines', 1, 'area_less', .99, 'area_more', .01);
            options = load_key_value_input_pairs(options, varargin{:});

            if ~ IS.percent([options.area_less, options.area_more])
                error('Limits must be percents in [0,1] that refer to the size of the original image where the patch comes from.')
            end                
            
            keep    = true(size(obj), 1);
            for i = 1:size(obj)
                pi = obj.collection(i);
                
                if options.no_lines && (pi.width() == 1 || pi.height() == 1)        % Remove line - patches.
                    keep(i) = false;
                    continue;
                end
                
                if options.area_less ~= 1                                       % Remove small area-wise patches.
                    if pi.area() >= options.area_less * obj.image.area();
                        keep(i) = false;
                        continue;
                    end
                end
                
                if options.area_more ~= 0                                       % Remove big area-wise patches.
                    if pi.area() <= options.area_more * obj.image.area();
                        keep(i) = false;
                        continue;
                    end
                end
                
            end           
        end
                
        function purge_patches(obj, purgatory_list)
            obj.collection = obj.collection(purgatory_list);
        end
        
        function set_collection(obj, new_patch_array, new_image)
            if isa(new_patch_array, 'Patch') && isa(new_image, 'Image')
                obj.collection = new_patch_array;
                obj.image      = new_image;
            else
                error('Incompatible types.')
            end            
        end          
        
        function [new_collection] = embed_in_new_image(obj, new_image)                       
            % TODO think of working with array of corners than Patches => use Patch_Collection constructor and not
            % set_collection.
            ow = obj.image.width;            
            oh = obj.image.height;
            nw = new_image.width;
            nh = new_image.height;            
            new_patches(size(obj)) = Patch();
            for i = 1:size(obj)                     
                new_patches(i) =  obj.get_patch(i).linear_interpolation(ow, nw, oh, nh);
            end          
            new_collection = Patch_Collection();
            new_collection.set_collection(new_patches, new_image);
        end
                
        function [A] = area_inside_mask(obj, bitmask)
            % Computes for every patch of the Patch_Collection, the area that resides inside the bit mask.
            %
            % Input:    
            %
            %           bitmask  - (m x n) binary matrix. Usually positions flagged with 1, correspond to ROI wrt. an image. E.g., 
            %                      they describe a segmentation of an object.
            %     
            % Output:                     
            %           A        - (N x 1) vector A(i) is the the area of the i-th patch (as returned 
            %                      from get_patch(i)) in the bitmask.                                    
            A = zeros(size(obj), 1);
                
            if ~ IS.binary(bitmask)
                error('Bitmask is not binary.');
            elseif ~ any(bitmask == 1)
                warning('Bitmask has zero elements set to 1.')            
                return                            
            end
                        
            for i = 1:size(obj)
                pi = obj.collection(i); 
                [xmin, ymin, xmax, ymax] = pi.get_corners();                
                A(i)  = sum(sum(bitmask(ymin:ymax, xmin:xmax)));
            end
        end
        
        function [overlaps] = overlap_with_patch(obj, patch)
            % Computes for every patch of the Patch_Collection, the fraction of its area that resides inside a given
            % patch.
            %
            % Input:                
            %           patch    - (Patch)
            %
            % Output:                     
            %           overlaps - (N x 1) vector overlaps(i) is the fraction of the area of the i-th patch (as returned 
            %                      from get_patch(i)) in the bitmask.                                                
            overlaps = zeros(size(obj), 1);
            for i = 1:size(obj)
                pi = obj.collection(i); 
                overlaps(i) = pi.area_of_intersection(patch) / pi.area();
                assert(overlaps(i) >= 0 && overlaps(i) <= 1);
            end
        end
        
        function [C] = corloc(obj, patch)
            C = zeros(size(obj), 1);             
            for i = 1:size(obj)                
                C(i) = obj.collection(i).corloc(patch);                
            end            
        end
        
        function [inter] = intersection(obj, patch)
            inter = zeros(size(obj), 1);             
            for i = 1:size(obj)                
                inter(i) = obj.collection(i).area_of_intersection(patch);
            end            
        end
        
        function [H] = weight_map(obj, normalized)
            [h, w] = size(obj.image);
            H = zeros(w, h);
            for i = 1:size(obj)                
                [xmin, ymin, xmax, ymax] = obj.collection(i).get_corners();
                H(ymin:ymax, xmin:xmax) = H(ymin:ymax, xmin:xmax) + 1;               
            end            
            if exist('normalized', 'var') && normalized
                H = H ./ sum(H(:)) ;
            end
        end
        
        function new_collection = keep_only(obj, keep_list)
            new_collection = Patch_Collection();
            new_collection.set_collection(obj.collection(keep_list), obj.image);
        end
        
        function [P] = patches_as_matrix(obj)
            s =  size(obj);
            P = zeros(s, 4);
            for i = 1:s
                P(i,:) = obj.collection(i).get_corners();
            end
        end
        
        function [p] = closest_to_gt(obj, top_k)            
            gt = obj.image.gt_segmentation;
            p = [];            
            if ~ isempty(gt)                
                best_per_gt = zeros(size(gt), top_k);
                for g = 1:size(gt)
                    C   = obj.corloc(gt.get_patch(g));
                    [~, ids] = sort(C, 'descend');
                    best_per_gt(g, :) = ids(1:top_k);                    
                end                
                p = obj.keep_only(unique(best_per_gt(:)));
            end        
        end
        
        
        function obj = merge_with(obj, other_collection)
            new_patches = size(other_collection);
            obj.collection(end+1: end+new_patches) = other_collection.collection(:);
        end
     
        function C = corloc_with_gt(obj)
            C = zeros(size(obj), 1);
            gt = obj.image.gt_segmentation;
            
            if size(gt) == 1                                     
                for g = 1:size(gt)
                    C = max(C, obj.corloc(gt.get_patch(g)));
                end            
            elseif size(gt) > 1            
                
                for p = 1:size(obj)                    % Measure how well you do only in those that you hit.
                    p_patch        = obj.collection(p);
                    inter_with_all = gt.intersection(p_patch);
                    has_overlap    = find(inter_with_all);
                    total_gt_area  = 0;
                    
                    for g = 1:size(has_overlap)
                        total_gt_area = total_gt_area + gt.get_patch(has_overlap(g)).area();            
                    end
                    
                    with_all_corloc   = sum(inter_with_all) ./ ((double(total_gt_area) + p_patch.area() - sum(inter_with_all)));
                    if with_all_corloc > 1
                       with_all_corloc = 1;  % Can happen when two gt' overlap (i.e. their tightest_boxes)
                    end
                    C(p) = max(C(p), with_all_corloc);
                end                               
                
                gt_all    = zeros(gt.image.height, gt.image.width);
                ones_mask = ones(gt.image.height, gt.image.width);                
                for j=1:size(gt)                    
                    gt_all = gt_all + gt.get_patch(j).extract_features(ones_mask , 'zero_pad');
                end
                gt_all(gt_all > 0) = 1;
                gt = Patch.tightest_box_of_segment(gt_all);                                
                C  = max(C, obj.corloc(gt));
            end                           

        end
        
    end
    
    
    

end
