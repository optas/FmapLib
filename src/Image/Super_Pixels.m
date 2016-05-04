classdef Super_Pixels
    properties (SetAccess = private)
        I;          % (Image) Image object where Super_Pixels are defined.
        mask;       % Matrix indicating for each pixel of I in which super-pixel it belongs.
        states;     % (int) number of Super_Pixels.
    end
            
    methods (Access = public)
        % Class Constructor.
        function obj = Super_Pixels(in_image, mask)
            if  nargin == 0
                obj.I      = Image();
                obj.mask   = [];
                obj.states = 0;
            else
                obj.I      = in_image;
                obj.mask   = mask;
                numbering  = unique(mask)';
                if ~all(numbering  == 1:length(numbering))
                    error('A super pixel mask is expected to be comprised by contiguous integers, starting from 1.')
                end                    
                obj.states = length(numbering);
            end
        end
        
        function h = plot(obj)            
            h = imagesc(obj.mask);            
        end
        
        function [V] = super_pixels_to_pixels(obj, values)
            % Given a real value associated with every super-pixel, computes a mask over the image,
            % where each location/pixel is equal with the real value of the corresponding SP.
            if length(values) ~= obj.states
                error('Input values do not match the dimensions of the super pixels.')
            end            
            V = zeros(size(obj.mask));            
            for i = 1:length(values)
              V(obj.mask == i) = values(i);
            end            
        end
        
        function C = image_content_on_sp(obj, sp_id)
            sp_mask = obj.mask == sp_id;              
            C = obj.I.apply_mask(sp_mask);            
            p = Patch.tightest_box_of_segment(sp_mask);
            C = p.extract_features(C, 'mask');                        
        end
        
        function [F] = mean_features_on_super_pixels(obj, pixel_features)
            % Given feature vectors defined on the pixels of an image, define a feature vector for each super-pixel
            % by considering the mean of the features defined on the pixels that constitute itself.                
            [h, w, fdim]  = size(pixel_features);
            if w ~= obj.I.width || h ~= obj.I.height
                error('Pixel Features are not defined on an image of the same size.');
            end
            F = zeros(obj.states, fdim);                    
            for s = 1:obj.states
                [x, y] = find(obj.mask == s);                         
                aggregator = zeros(length(x), fdim);
                for m = 1:length(x)
                    aggregator(m,:) =  pixel_features(x(m), y(m), :);
                end
                F(s,:) = mean(aggregator);
            end                
        end
                
        function H = color_histograms(obj, num_bins)                                       
               bounds = 1: (256/num_bins) : 256;   % TODO. formalize.
               bounds(end+1) = 256;               
               H = zeros(obj.states, num_bins);
               for s = 1:obj.states
                   all_vals = obj.I.CData(obj.mask == s);                         
                   H(s,:) = histcounts(all_vals, bounds, 'Normalization', 'probability');               
               end
        end
        
        
        
    end % Object's Public Methods.
end