classdef Image_Features < dynamicprops      
%     TODO-Z: is it working nicely on pixel - level ?
    methods (Static, Access = public)
        
        function [signatures] = hog_signature(in_image)
% TODO- only send into extractHOGF the points inside the original image.

%                 if in_image.is_rgb()    TODO-P: 
%                      content = rgb2gray(in_image.CData);                             
%                 else
%                     content = in_image.CData;
%                 end
                              
                h = in_image.height;
                w = in_image.width;
                pad_size  = 10;
                
%                 I1 = [zeros(pad_size, 2*pad_size + w); [zeros(h, pad_size) content zeros(h, pad_size)]; zeros(pad_size, 2 * pad_size + w)];                
                I1 = padarray(in_image.CData, [pad_size, pad_size] );
                                
                [y, x] = find(ones(size(I1, 1), size(I1, 2)));              % nope - only find original pixels.
                                                              
                [v, valid_points] = extractHOGFeatures(I1, [x y]);                                
                
                f = nan([h, w , size(v, 2)]);
                                
                for p = 1:length(valid_points)                   
                    f(valid_points(p, 2), valid_points(p, 1), :) = v(p,:);  % swap coordinates                
                end
                                                
                signatures = f(pad_size+1:pad_size+h, pad_size+1:pad_size+w, :);
%                 assert( any(isnan(signatures(:))) == 0)
                if any(isnan(signatures(:)))
                    warning('oops. bad hog!')
                end
                    
        end
        
    end
end


