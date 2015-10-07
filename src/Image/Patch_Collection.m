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
                obj.image = Image();
            else
                m = size(corners, 1);
                if m < 1 
                    error('A Patch Collection must consist of at least one Patch objects.')
                end                
                obj.image = in_image;
                obj.collection(1:m) = Patch();                                
                for i = 1:m
                    if ~ Patch.is_valid(corners(i,:), in_image)
                        error('Patch %d is cannot go with given image.\n', i);
                    else
                        obj.collection(i) = Patch(corners(i, :));
                    end                    
                end
            end
        end
        
        function m = size(obj)
            m = length(obj.collection);
        end
        
        function p = get_patch(obj, p_index)
            p = obj.collection(p_index) ;
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
        
        function [corloc] = corloc(obj, patch)
            corloc = zeros(size(obj), 1);             
            for i = 1:size(obj)                
                corloc(i) = obj.collection(i).corloc(patch);                
            end            
        end
        
     
    end
    
    
    

end
