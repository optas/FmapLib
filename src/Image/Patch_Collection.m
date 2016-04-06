classdef Patch_Collection < dynamicprops
    
    properties (GetAccess = public, SetAccess = private)
        collection = Patch;     %    (Patch Array) storing all the patches.
        image;                  %    (Image) the image over which the patches were defined.
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
            % TODO: implement subref.
            p = obj.collection(p_index);
        end

        function F = plot_collection(obj, scores, add_text)
            obj.image.plot();
            hold on;            
            if exist('scores', 'var')
                colors = vals2colormap(scores);
                if exist('add_text', 'var') && add_text
                    for i = 1:length(scores)                    
                        F = obj.get_patch(i).plot('color', colors(i,:), 'text', sprintf('%.2f', scores(i)));                    
                    end                
                else                    
                    for i = 1:length(scores)                    
                        F = obj.get_patch(i).plot('color', colors(i,:));                    
                    end
                end
            else                    
                for i = 1:size(obj)
                    F = obj.get_patch(i).plot();
                end                
            end
        end
           
        function F = plot_patches(obj, patch_ids, scores)
            obj.image.plot();
            hold on;                               
            if exist('scores', 'var')
                colors = vals2colormap(scores);
                for i = 1:length(scores)                                        
                    F = obj.get_patch(patch_ids(i)).plot('color', colors(i,:));                    
                end
            else
                for i = 1:length(patch_ids)                
                    F = obj.get_patch(patch_ids(i)).plot();
                end                            
            end
        end
        
        function obj = add_patch_property(obj, prop_id, prop_content)
            if size(prop_content, 1) ~= size(obj)
                error('Give content for every patch (add empty if needed).')
            end
            if ~ isprop(obj, prop_id)
                obj.addprop(prop_id);
            end
            obj.(prop_id) = prop_content;
        end
        
        function R = rects(obj, ids)
            if isprop(obj, 'rects_prop')
                R = obj.rects_prop;
            else
                n = size(obj);
                c = obj.collection;
                R = zeros(n,4);
                for i = 1:n
                    R(i,1) = c(i).corners(1);
                    R(i,2) = c(i).corners(2);                
                    R(i,3) = c(i).height;
                    R(i,4) = c(i).width;
                end
                obj.addprop('rects_prop');
                obj.('rects_prop') = R;
            end
            if nargin == 2
                R = R(ids,:);
            end
            
        end
        
        function [keep] = filter_patches(obj, varargin)
            % Keep patches that are not area-wise too big or too small wrt. to their image.
            % Also removes line patches.
            
            options = struct('no_lines', 1, 'area_less', (.99)^2, 'area_more', (.01)^2);
            options = load_key_value_input_pairs(options, varargin{:});

            if ~ IS.percent([options.area_less, options.area_more])
                error('Limits must be percents in [0,1] that refer to the size of the original image where the patch comes from.')
            end                
            
            keep    = true(size(obj), 1);
            for i = 1:size(obj)
                pi = obj.collection(i);
                
                if options.no_lines && (pi.width() == 1 || pi.height() == 1)    % Remove line - patches.
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
            % Computes for every patch of the Patch_Collection, the area (num pixels) that resides inside the bit mask.
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
        
        function I = mask_image_outside_patches(obj)
            I = obj.image;            
            mask = zeros(I.height, I.width);
            for i = 1:size(obj)
                pi = obj.get_patch(i);
                [xmin, ymin, xmax, ymax] = pi.get_corners;
                mask(ymin:ymax, xmin:xmax) = 1;
            end 
            I = I.apply_mask(mask);            
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
        
        function [H] = weight_map(obj, masses, normalized)            
            H = zeros(size(obj.image));
            for i = 1:size(obj)                
                [xmin, ymin, xmax, ymax] = obj.collection(i).get_corners();
                if exist('masses', 'var')
                    H(ymin:ymax, xmin:xmax) = H(ymin:ymax, xmin:xmax) + masses(i);               
                else
                    H(ymin:ymax, xmin:xmax) = H(ymin:ymax, xmin:xmax) + 1;               
                end
            end            
            if exist('normalized', 'var') && normalized
                H = H ./ sum(H(:)) ;
            end
        end
        
        function [H] = color_map(obj, masses)            
            H = zeros(size(obj.image));
            for i = 1:size(obj)                
                [xmin, ymin, xmax, ymax] = obj.collection(i).get_corners();
                H(ymin:ymax, xmin:xmax) = max(H(ymin:ymax, xmin:xmax),  masses(i));                              
            end                       
        end
                
        function new_collection = keep_only(obj, keep_list)
            new_collection = Patch_Collection();
            new_collection.set_collection(obj.collection(keep_list), obj.image);
        end
        
        function [P] = corners(obj)
            s =  size(obj);
            P = zeros(s, 4);
            for i = 1:s
                P(i,:) = obj.collection(i).get_corners();
            end
        end
        
        function [varargout] = closest_to_gt(obj, top_k)            
            gt = obj.image.gt_segmentation;
            p = [];            
            if ~ isempty(gt)                
                best_per_gt = zeros(size(gt), top_k);
                for g = 1:size(gt)
                    C   = obj.corloc(gt.get_patch(g));
                    [~, ids] = sort(C, 'descend');
                    best_per_gt(g, :) = ids(1:top_k);                    
                end               
                p = unique(best_per_gt(:), 'stable');                                
            end
            varargout = cell(nargout);                
            varargout{1} = p;
            if nargout == 2                
                varargout{2} = obj.keep_only(p);
            end               
        end
        
        function A = areas(obj, patch_id)
            if nargin == 1 % Return the areas of every patch.                
                if size(obj) == 1
                    A = obj.collection(1).area();
                else
                    A = arrayfun( @(p) p.area(), obj.collection);
                end                
            else
                A = arrayfun( @(p) p.area(), obj.collection(patch_id));                
            end
        end
        
        function C = centers(obj, patch_id)
            if nargin == 1 % Return the areas of every patch.                
                if size(obj) == 1
                    C = obj.collection(1).center();
                else
                    C = arrayfun( @(p) p.center(), obj.collection, 'UniformOutput', false);
                end                
            else
                C = arrayfun( @(p) p.center(), obj.collection(patch_id));                
            end
            C = cell2mat(C);
        end
        
        function C = diagonal_lengths(obj, patch_id)
            if nargin == 1 % Return the areas of every patch.                
                if size(obj) == 1
                    C = obj.collection(1).diagonal_length();
                else
                    C = arrayfun( @(p) p.diagonal_length(), obj.collection);
                end                
            else
                C = arrayfun( @(p) p.diagonal_length(), obj.collection(patch_id));                
            end
        end
        
        
        function obj = merge_with(obj, other_collection)
            new_patches = size(other_collection);
            obj.collection(end+1: end+new_patches) = other_collection.collection(:);
        end
     
        function C = corloc_with_gt(obj, options)
            if exist('options', 'var') && isfield(options, 'lax_corloc')
                lax_corloc = options.lax_corloc;
            else
                lax_corloc = false;
            end
            
            C  = zeros(size(obj), 1);
            gt = obj.image.gt_segmentation;

            if size(gt) == 1                                     
                C = obj.corloc(gt.get_patch(1));                            
            elseif size(gt) > 1 && ~lax_corloc
                for g = 1:size(gt)
                    C = max(C, obj.corloc(gt.get_patch(g)));                            
                end                
            else           
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
        
        function C = equivalence_classes_in_area(obj, num_classes)          
            if num_classes == 1
                C = ones(size(obj), 1);
                return
            end
            areas = obj.areas;
            C     = zeros(length(areas), 1);
            mass_in_class = 100 / num_classes;
            prcs = mass_in_class: mass_in_class: (100-mass_in_class);
            prcs = prctile(areas, prcs);
            prcs = [0, prcs];
            class_ind = 1;
            for i = 1:length(prcs)-1                
                in_class = (areas >= prcs(i) & areas < prcs(i+1));
                if sum(in_class) > 0    % This conditions complicates coding logic but is necessary for extreme input cases.
                    C = C + class_ind * in_class;
                    class_ind  = class_ind  + 1;
                end                    
            end            
            
            i = i + 1;          % Largest/last class defined only by being bigger than last percentile.
            C = C + class_ind * (areas >= prcs(i));
            
            assert(all(C>0) & all(C<=num_classes))               
            for i = 1:class_ind-1
                assert(max(areas(C==i)) <= min(areas(C==i+1)))                
            end            
        end
              
        function I = inclusions(obj, contain, more_than, how_many)
            n = size(obj);            
            I = cell(n,1);
            areas = obj.areas;                        
            R     = obj.rects;
            for p = 1:n                
                bigger_boxes = (areas > more_than .* areas(p));                  
                inter_areas  = rectint(R(p,:), R);                
                bigger_boxes = bigger_boxes & (inter_areas > (contain * areas(p)))';  % Discard big boxes that do not contain enough of p.                                                  
                I{p}         = find(bigger_boxes);
                if isempty(I{p})
                    continue
                end
                [~, ids]     = sort(areas(I{p}), 'ascend');
                if exist('how_many', 'var') && ~strcmp(how_many, 'all')                    
                    keep = min(how_many, length(ids)); 
                else
                    keep = length(ids);                                               % Discover all inclusions.
                end                
                I{p}         = I{p}(ids(1:keep));
            end                                               
        end
        
        function S = most_similar_inclusions(obj, descriptors, topk, inclusions)
            S = zeros(size(obj), topk);                                
            all_dists = pdist2_vec(descriptors, descriptors);            % Due to vectorization of pdist2_vec computing
            for p = 1:size(obj)                                          % all distances is faster than otherwise.
                if isempty(inclusions{p})
                    continue;
                end
                [~, ids]      = sort(all_dists(p, inclusions{p}));
                keep          = min(topk, length(ids));                
                S(p,1:keep)   = inclusions{p}(ids(1:keep));
            end
        end
        
        function I = pw_area_intersections(obj)
            % Computes all pairwise intersections between the patches stored in the collection.
            I = single(rectint(obj.rects, obj.rects));            
        end
        
        
        function S = symmetries(obj, descriptors, topk, intersection_thres, side_thres)
            if ~ exist('intersection_thres', 'var')
                intersection_thres = 0.1;
            end
            if ~ exist('side_thres', 'var')
                side_thres         = 0.05; 
            end                         
            R = obj.rects;        
            A = obj.areas;    
            
            [~, small_patch] = sort(A, 'ascend');               % TODO: parameterize on function call.
            small_patch      = small_patch(1: ceil(0.1*size(obj))); 
            
            S = zeros(size(obj), topk);
            all_pairs = pdist2_vec(descriptors, descriptors);            
            for p = 1:size(obj)
                if ~any(small_patch == p)
                    continue
                end
                proper_inter   = (rectint(R(p,:), R) ./ A(p)) < intersection_thres;             % Intersection constraint.                 
                proper_side    = proper_inter' & (abs(R(p,3) - R(:,3)) ./ R(p,3) < side_thres); % Area constraint.               
                proper_side    = proper_side   & (abs(R(p,4) - R(:,4)) ./ R(p,4) < side_thres);
                proper_side(p) = 0;                   
                proper_side    = find(proper_side);
              
                if isempty(proper_side)
                    continue;
                end                   
                
                [~, ids]     = sort(all_pairs(p, proper_side)); 
                keep         = min(topk, length(ids));
                S(p,1:keep)  = proper_side(ids(1:keep));
            end                                
            clear all_pairs;            
        end
                       

        function P = inside_gt(obj, topk)
            % topk patches that are covered the most by the gt(s) of th collection.            
            gt = obj.image.gt_segmentation;
            areas = obj.areas();            
            overlaps = zeros(size(obj), 1);                        
            for g = 1 : size(gt)
                bitmask  = obj.image.patch_indicator(gt.get_patch(g));                
                overlaps = max(overlaps, obj.area_inside_mask(bitmask));
            end            
            overlaps = overlaps ./ areas;   % Normalize to get fraction of area inside gt.                      
            [~, top_ids] = sort(overlaps, 'Descend');
            P = top_ids(1:topk);
        end
        
        function P = parts_of(obj, frame_ids)
            pw_int = rectint(obj.rects, obj.rects(frame_ids));
            areas  = obj.areas;
            P      = cell(length(frame_ids), 1);
            for i = 1:length(frame_ids)                
                P{i} = setdiff(find(pw_int(:, i) ./ areas == 1), frame_ids(i));
            end            
            P = unique(cell2mat(P(~cellfun('isempty', P))));
        end
        
        
        
        %%  %% %%  %% %%  %% %%  %%  %%  %% %%  %%  %%  %% %%  %% %%  %% %%  %%  %%  %% %%  %%
        %%              Functions deriving Features on Patches.                             %%
        %%  %% %%  %% %%  %% %%  %%  %%  %% %%  %%  %%  %% %%  %% %%  %% %%  %%  %%  %% %%  %%
        
        function [F] = hog_features(obj, cell_size, n_orients, scaling)
            % Fill in unset optional values.
            if nargin == 1                
                cell_size  = 8;
                n_orients  = 9;
                scaling    = 64;
            end                                          
            I        = obj.image;
            patches  = obj.collection;
            hog_size = (floor(scaling / cell_size))^2 * 4 * n_orients;                       
%             temp =  extractHOGFeatures(imResample(single(I.content_in_patch(patches(1))), [scaling scaling]), 'CellSize', [cell_size, cell_size], 'NumBins', n_orients);            
%             hog_size = length(temp);
            F        = zeros(size(obj), hog_size);
            for p = 1:size(obj)
                content = I.content_in_patch(patches(p));
                box     = imResample(single(content), [scaling scaling]) / 255; 
%                 sig     = vl_hog(box, cell_size, 'numOrientations', n_orients, 'variant', 'dalaltriggs');
%                 sig     = extractHOGFeatures(box, 'CellSize', [cell_size, cell_size], 'NumBins', n_orients);
                sig     = hog(box, cell_size, n_orients);
                F(p,:)  = vec(sig);
            end
        end
        
        function [F] = neural_net_features(obj, net, layer_id)
            I           = obj.image; 
            I           = I.im2single();
            nn_im_size  = net.meta.normalization.imageSize(1:2);
            net.layers  = net.layers(1:end-layer_id);
            num_patches = size(obj);
            patches     = obj.collection;
            im          = single(zeros(nn_im_size(1), nn_im_size(2), 3, num_patches));
            av          = repmat(net.meta.normalization.averageImage, nn_im_size(1), nn_im_size(2));            
            for i = 1:size(obj)                
                im(:,:,:,i) = imresize(I.content_in_patch(patches(i)), nn_im_size) - av;
            end                        
            net.info.opts.batchSize = 40;
            temp = vl_simplenn(net, im);            
            F = squeeze(gather(temp(end).x))';
            clear temp;
        end
              
    end
    
    
    

end
