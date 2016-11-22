classdef Patch_Collection < dynamicprops
    
    properties (GetAccess = public, SetAccess = private)
        collection = Patch;     %    (Patch Array) 1-dimensional array storing Patch objects.
        image;                  %    (Image) the image over which the patches were defined.
    end

    methods (Access = public)                
        function obj = Patch_Collection(corners, in_image)
            % Class Constructor.
            % Input:
            %           corners  - (N x 4 int Matrix) carrying the (x,y) coordinates of the 2 extreme corners defining N
            %                      rectangular patches. The expected format for each row is: [xmin, ymin, xmax, ymax].
            %           in_image - (Image) The image over which the patches were defined.
            if nargin == 0
                obj.collection = Patch();
                obj.image = Image();
            else
                obj.image = in_image;                
                if isa(corners, 'single') || isa(corners, 'double') 
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
                    error('Not correct arguments.')
                end
            end
        end
        
        function m = size(self)
            if length(self) > 1 % An array of Patch_Collection objects
                m = length(self);                
            elseif size(self.collection) == 0 % Empty patch
                m = 0;
            else
                m = length(self.collection);
            end
        end
        
        function p = get_patch(obj, p_index)
            % TODO: implement subref.
            p = obj.collection(p_index);
        end
        

        function [S] = intersection_over_union(self, patch)            
            S = bboxOverlapRatio(self.rects, patch.as_rectangle);           
            % If bboxOverlapRatio is not defined.
            % r = patch.as_rectangle();
            % pw_int = rectint(self.rects, r);
            % S = pw_int ./ (self.areas() + patch.area() - pw_int);        
        end
        
        function [inter] = intersection(self, patch)
            r = patch.as_rectangle();
            inter = rectint(self.rects, r);                      
        end
        
        function [overlaps] = overlap_with_patch(self, patch)
            % Computes for every patch of the Patch_Collection, the fraction of its area that resides inside a given
            % patch.
            %
            % Input:                
            %           patch    - (Patch)
            %
            % Output:                     
            %           overlaps - (N x 1) vector overlaps(i) is the fraction of the area of the i-th patch in the  patch.                                                                         
            overlaps = self.intersection(patch) ./ self.areas;           
            assert(all(overlaps >= 0) && all(overlaps <= 1));
        end
                
        function [P] = closest_to_gt(self, top_k)
            % Returns the indices of the patches that have biggest IOU statistic
            % with any of the groundtruth bounding boxes of the image.
            S = self.iou_with_gt();
            [~, ids] = sort(S, 'descend');            
            effective_k = min(length(ids), top_k);
            P = ids(1:effective_k );            
        end
        
        function S = iou_with_gt(self)
            % Intersection Over Union (IOU) statistic for every patch of the collection           
            % wrt. to the ground-truth bounding box of the underlying image. If more
            % than one groun-truth boxes exist, it reports the maximum IOU a patch has with 
            % any of them.
            gt = self.image.gt_segmentation;            
            S = [];            
            if size(gt) == 1
                S = self.intersection_over_union(gt.get_patch(1));
            elseif size(gt) > 1
                S = zeros(size(self), 1);
                for g = 1:size(gt)
                    S = max(S, self.intersection_over_union(gt.get_patch(g)));                            
                end                            
            end
        end
                        
        function F = plot_collection(obj, scores, add_text)
            obj.image.plot();
            hold on;            
            if exist('scores', 'var')
                colors = vals2colormap(scores, 'jet');                
%                 colors = ceil(colors .* 255);                
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
                
        function [H, M] = weight_map(obj, masses)            
            H = zeros(size(obj.image));
            M = zeros(size(obj.image));
            for i = 1:size(obj)
                [xmin, ymin, xmax, ymax] = obj.collection(i).get_corners();
                if exist('masses', 'var')
                    M(ymin:ymax, xmin:xmax) = M(ymin:ymax, xmin:xmax) + masses(i);
                end
                H(ymin:ymax, xmin:xmax) = H(ymin:ymax, xmin:xmax) + 1;
            end            
        end
        
        function [H] = color_map(obj, masses)            
            H = zeros(size(obj.image));            
            for i = 1:size(obj)                
                [xmin, ymin, xmax, ymax] = obj.collection(i).get_corners();
                H(ymin:ymax, xmin:xmax) = H(ymin:ymax, xmin:xmax) +  masses(i);                              
            end                       
        end
                
        function set_collection(obj, new_patch_array, new_image)
            if isa(new_patch_array, 'Patch') && isa(new_image, 'Image')
                obj.collection = new_patch_array;
                obj.image = new_image;
            else
                error('Incompatible types.')
            end            
        end          
        
        function new_collection = keep_only(obj, keep_list)
            new_collection = Patch_Collection();
            new_collection.set_collection(obj.collection(keep_list), obj.image);
        end
        
        function new_collection = keep_top_scoring(self, scores, top_k)
            if length(scores) ~= self.size
                error('Provide a score for every proposal.');
            end
            [~, top_ids] = sort(scores, 'descend');
            new_collection  = self.keep_only(top_ids(1:top_k));
        end
                        
        function obj = merge_with(obj, other_collection)
            new_patches = size(other_collection);
            obj.collection(end+1: end+new_patches) = other_collection.collection(:);
        end
     
                
        function C = area_equipartition(obj, num_classes)          
            if num_classes == 1
                C = ones(size(obj), 1);
                return
            end
            areas = obj.areas;
            C = zeros(length(areas), 1);
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
        
        function C = area_geometric_partition(self, separation)
            [s_areas, s_ids] = sort(self.areas, 'ascend');
            membership = zeros(size(self), 1);
            membership(1) = 1;
            min_of_class = 1;
            class_id = 1;
            while true
                bound = s_areas(min_of_class) * separation;
                last_of_class = find(s_areas > bound, 1);
                                
                if isempty(last_of_class)                    
                    membership(min_of_class:end) = class_id;                    
                    break
                end
                
                membership(min_of_class:last_of_class-1) = class_id;
                min_of_class = last_of_class;
                class_id = class_id + 1;                
            end
            assert(all(membership > 0))
            C = zeros(length(membership),1);
            C(s_ids) = membership;            
            areas = self.areas;
            for i = 1:max(C)-1
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
            % topk patches that are covered the most by the gt(s) of the collection.
            gt = obj.image.gt_segmentation;
            areas = obj.areas();
            overlaps = zeros(size(obj), 1);
            for g = 1 : size(gt)
                bitmask  = gt.get_patch(g).indicator_on_image(obj.image);                
                overlaps = max(overlaps, obj.area_inside_mask(bitmask));
            end            
            overlaps = overlaps ./ areas;   % Normalize to get fraction of area inside gt.                      
            [~, top_ids] = sort(overlaps, 'Descend');
            P = top_ids(1:topk);
        end
        
        function P = parts_of(obj, patch_id)
            pw_int = rectint(obj.rects, obj.rects(patch_id));
            areas  = obj.areas;            
            P      = setdiff(find(pw_int ./ areas == 1), patch_id);            
            if size(P,1) ~= 1
                P = P';
            end
        end
        
        
        function [top_ids, top_scores] = k_best_non_overlaping_patches(obj, scores, k, overlap_ratio)
            if ~ exist('overlap_ratio', 'var')
                overlap_ratio = 0.8;
            end            
            R = obj.rects;
            A = obj.areas;                              
            top_ids = []; top_scores = [];
            [max_score, top_id] = max(scores);
            while numel(top_ids)<k && ~(max_score==-inf)
                top_ids        = [top_ids, top_id];
                top_scores     = [top_scores, max_score];
                scores(top_id) = -inf;
                id_valid       = find(isfinite(scores));
                area_int       = rectint(R(top_id, :), R(id_valid,:))';                
                IOU            = area_int ./ (A(top_id) + A(id_valid) - area_int);
                id_nm          = id_valid(find(IOU > overlap_ratio));
                scores(id_nm)  = -inf;
                [max_score, top_id] = max(scores);    
            end
        end

        
        %%  %% %%  %% %%  %% %%  %%  %%  %% %%  %%  %%  %% %%  %% %%  %% %%  %%  %%  %% %%  %%
        %%              Functions on basic properties of the Patch collection.              %%
        %%  %% %%  %% %%  %% %%  %%  %%  %% %%  %%  %%  %% %%  %% %%  %% %%  %%  %%  %% %%  %%        
        
        function [P] = corners(self)
            n = size(self);
            patches = self.collection;
            P = single(zeros(n, 4));          
            for i = 1:n
                P(i,:) = patches(i).corners;
            end
        end
        
        function [R] = rects(self, ids)
            if isprop(self, 'rects_prop')
                R = self.rects_prop;
            else
                n = size(self);
                c = self.collection;
                R = single(zeros(n, 4));
                for i = 1:n
                    R(i,:) = c(i).as_rectangle();
                end
                self.addprop('rects_prop');
                self.('rects_prop') = R;
            end
            if nargin == 2
                R = R(ids,:);
            end            
        end
        
        function [I] = pw_area_intersections(self)
            % Computes all pairwise intersections between the patches stored in the collection.
            I = single(rectint(self.rects, self.rects));            
        end
        
        
        function [S] = pw_iou(self)
            % Computes all pairwise intersections between the patches stored in the collection.
            I = self.pw_area_intersections();
            A = self.areas;
            U = LA.outersum(A, A);
            U = U - I;
            S = I ./ U;            
        end
        
        
        function [A] = areas(self, ids)            
            if nargin == 1 % Return the areas of every patch.                
                R = self.rects;                
            else
                R = self.rects(ids);
            end
            A = R(:,3) .* R(:,4);           
        end
        
        function [C] = centers(obj, patch_id)
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
        
        function [C] = diagonal_lengths(obj, patch_id)
            if nargin == 1
                if size(obj) == 1
                    C = obj.collection(1).diagonal_length();
                else
                    C = arrayfun( @(p) p.diagonal_length(), obj.collection);
                end                
            else
                C = arrayfun( @(p) p.diagonal_length(), obj.collection(patch_id));
            end
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
              
    end % End of public methods.
    
    methods (Static, Access = public)
        function [kept, thrown] = cho_proposal_filter(proposals, image)
            % Cho et al. CVPR 2015 apply this filtering to the proposals generated via the RP method on the dataset of VOC_2007x2.            
            width = image.width;
            height = image.height;
            c1 = proposals(:,1) > width  * 0.01 & proposals(:,3) < width  * 0.99 ...    % Proposal has to be relatively centered.
               & proposals(:,2) > height * 0.01 & proposals(:,4) < height * 0.99;
            c2 = proposals(:,1) < width  * 0.01 & proposals(:,3) > width  * 0.99 ...    % Or it has to be very big.
               & proposals(:,2) < height * 0.01 & proposals(:,4) > height * 0.99;           
            valid = c1 | c2;
            kept = proposals(valid, :);
            thrown = proposals(~valid, :);
        end
    end
    
    
end % End of class.
