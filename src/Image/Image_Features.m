classdef Image_Features < dynamicprops      
    
    methods (Static, Access = public)        
                
        function [signatures] = hog_signature(in_image)
            h = in_image.height;
            w = in_image.width;
            block_size = [1,1];
            [hog_feats] = extractHOGFeatures(in_image.CData);
            [hogfeatures, visualization] = extractHOGFeatures(images(1).resize(100,100).CData, [10,10], 'BlockSize', block_size);
        
        end
                       
        function [signatures] = hog_signature_2(in_image)
            h = in_image.height;
            w = in_image.width;
            pad_size  = 10;      % TODO - remove fixed size.

            padded_im = padarray(in_image.CData, [pad_size, pad_size] );
            pad_mask  = padarray(ones(size(in_image.CData)), [pad_size, pad_size] );
            [x, y]    = find(pad_mask);                                                    % Returns only points that lie inside original image.

            [hog_feats, valid_points] = extractHOGFeatures(padded_im, [y x]);              % Need to swap for matrix coordinates to represent image's one.                                    

            feat_carrier = nan([h, w , size(hog_feats, 2)]);                                                            
            for p = 1:length(valid_points)                   
                feat_carrier(valid_points(p, 2), valid_points(p, 1), :) = hog_feats(p, :);
            end

            signatures = feat_carrier(pad_size+1:pad_size+h, pad_size+1:pad_size+w, :);  % Only keep those inside original image.

            if any(isnan(signatures(:)))
                warning('Some pixels were not assigned a hog descriptor.')
            end                    
        end
        
        function [signatures] = local_binary_pattern_signatures(in_image)
            signatures = efficientLBP(in_image.CData);
            signatures  = 1 - double(signatures) / 255;                         % "1 minus signature" to put high values on white pixels (as sift).
        end

        function [signatures] = sift_signature_2(in_image)
            cellsize    = 3; 
            gridspacing = 1;
            I = im2double(in_image.CData);
            signatures = im2double(mexDenseSIFT(I, cellsize, gridspacing));
        end
        
        function [signatures] = sift_signature(in_image)
            SIFTparam.grid_spacing = 1;     % Distance between grid centers.
            patch_sizes = [8, 12, 16];
            [w, h] = size(in_image);
            signatures = [];
            for k = 1:length(patch_sizes)
                SIFTparam.patch_size = patch_sizes(k); % Size of patch which is used to compute SIFT (it has to be a factor of 4).
                pad_size = SIFTparam.patch_size/2 - 1;
                if in_image.is_rgb
                    gray = rgb2gray(in_image.CData); 
                else
                    gray = in_image.CData;
                end
                I1 = [zeros(pad_size, 2*pad_size + w); [zeros(h, pad_size) gray zeros(h, pad_size)]; zeros(pad_size, 2 * pad_size + w)];
                v  = LMdenseSift(I1, '', SIFTparam);
                signatures = cat(3, signatures, double(v));
            end
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
        
        function [E] = gist_embedding(in_image)
            param.imageSize = [256 256];
            param.orientationsPerScale = [8 8 8 8];
            param.numberBlocks = 4;
            param.fc_prefilt = 4;
            E = LMgist(in_image.CData, '', param);
        end
        
        function [features] = compute_default_features(in_image, feature_types)
            if ~ isa(feature_types, 'cell')
                feature_types = {feature_types};
            end
            features = [];
            for i = 1:length(feature_types)
                switch feature_types{i}
                    case 'sift'
                        f = Image_Features.sift_signature(in_image);                        
                    case 'hog'
                        f = Image_Features.hog_signature(in_image);                        
                    case 'color'
                        f = in_image.color();
                    case 'lbp'
                        f = Image_Features.local_binary_pattern_signatures(in_image);
                    otherwise
                        error('Unknown type of feature was requested.')                    
                end                                
                features = cat(3, features, double(f));
            end
        end
        
    end
end


