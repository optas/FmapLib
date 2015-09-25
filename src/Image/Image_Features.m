classdef Image_Features < dynamicprops      
    methods (Static, Access = public)
        
        function [signatures] = hog_signature(in_image)                      
            h = in_image.height;
            w = in_image.width;
            pad_size  = 10;      % TODO - remove fixed size.

            padded_im = padarray(in_image.CData, [pad_size, pad_size] );
            pad_mask  = padarray(ones(size(in_image.CData)), [pad_size, pad_size] );
            [x, y]    = find(pad_mask);                                                    % Returns only points that lie inside original image.

            [hog_feats, valid_points] = extractHOGFeatures(padded_im, [y x]);            % Need to swap for matrix coordinates to represent image's one.                                    

            feat_carrier = nan([h, w , size(hog_feats, 2)]);                                
            for p = 1:length(valid_points)                   
                feat_carrier(valid_points(p, 2), valid_points(p, 1), :) = hog_feats(p, :);
            end

            signatures = feat_carrier(pad_size+1:pad_size+h, pad_size+1:pad_size+w, :);  % Only keep those inside original image.

            if any(isnan(signatures(:)))
                warning('Some pixels were not assigned a hog descriptor.')
            end                    
        end
        
        function [signatures] = sift_signature(in_image)
            SIFTparam.grid_spacing = 1; % Distance between grid centers.            
            SIFTparam.patch_size = 12;
            pad_size = SIFTparam.patch_size/2 - 1;
            padded_im = padarray(in_image.CData, [pad_size, pad_size]);
            signatures = LMdenseSift(padded_im, '', SIFTparam);            
        end
        
        function [signatures] = color_histogram(in_image, num_bins, varargin)
            % Histogram of color values of input image. 
            % Parameters
            %               num_bins   - number of bins of the resulting histogram.
            %               normalize  - if true, histogram is a probability distribution.
            %               grayscale  - if true, the image first is converted to an intensity one.
            options = struct('normalize', 1, 'grayscale', 1);
            options = load_key_value_input_pairs(options, varargin{:});
            
            I = in_image.CData;            
            if options.grayscale
                I          = rgb2gray(I);                
                signatures = imhist(I, num_bins);
            else                
                signatures        = zeros(num_bins, 3);
                signatures(:, 1)  = imhist(I(:, :, 1), num_bins);    % Red   channel.
                signatures(:, 2)  = imhist(I(:, :, 2), num_bins);    % Green channel.
                signatures(:, 3)  = imhist(I(:, :, 3), num_bins);    % Blue  channel.
            end            
            
            if options.normalize                
                signatures = divide_columns(signatures, sum(signatures, 1));
            end
            
        end
        
        
    end
end


