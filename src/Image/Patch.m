classdef Patch < dynamicprops
    
    properties (GetAccess = public, SetAccess = private)
        % Basic properties that every instance of the Patch class has.                
        
        corners     % (1 x 4 positive int) x,y coordinates of two extreme corner points
                    % defining a rectangular image patch. They are [xmin, ymin, xmax, ymax].                    
    end
       
    methods (Access = public)                
        function obj = Patch(varargin)
            % Class constructror:
            %   Input: 
            %             (1 x 4) vector describing the corners as: [xmin, ymin, xmax, ymax].
            if nargin == 0
                obj.corners = [];
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
           
        function [varargout] = size(self)                        
            if length(self) > 1     % Array of Patches.
                varargout{1} = length(self);               
                return
            end
            
            if isempty(self.corners)    % Empty patch.
                varargout{1} = 0;
                return 
            end
            
            if nargout == 2
                varargout = cell(nargout);
                varargout{1} = self.width;
                varargout{2} = self.height;
            else
                varargout{1} = [self.width, self.height];
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
        
        function w = width(self)                       
            [xmin, ~, xmax, ~] = self.get_corners();
            w = xmax - xmin + 1;            
        end
        
        function h = height(self)                       
            [~, ymin, ~, ymax] = self.get_corners();
            h = ymax - ymin + 1;            
        end
        
        function a = area(self)            
            [xmin, ymin, xmax, ymax] = self.get_corners();
            a = (xmax-xmin+1) * (ymax-ymin+1);
%             a = self.height() * self.width();
        end
        
        function R = as_rectangle(self)
            R = single(zeros(1,4));
            R(1) = self.corners(1);
            R(2) = self.corners(2);
            R(3) = self.width;
            R(4) = self.height;
        end
            
        function c = center(self)
            [xmin, ymin, xmax, ymax] = self.get_corners();
            c = [(xmax+xmin)/2, (ymax+ymin)/2];
        end
        
        function d = diagonal_length(self)
            [xmin, ymin, xmax, ymax] = self.get_corners();
            d = norm(xmax - xmin, ymax - ymin);
        end
        
        function [F] = plot(obj, varargin)
            % Plots the boundary of the patch.
            options = struct('image', [], 'color', 'r', 'line_width', 3, 'text', []);
            options = load_key_value_input_pairs(options, varargin{:});         
            
            [xmin, ymin, xmax, ymax] = obj.get_corners();
            if ~ isempty(options.image)
                image.plot();            
                hold on;
            end
            F = plot([xmin xmax xmax xmin xmin], [ymin ymin ymax ymax ymin], ...
                     'Color', options.color, 'LineWidth', options.line_width);
            
            if ~isempty(options.text)
                text(xmin+3, ymin+10, options.text, 'FontSize', 15);                 
            end
                                                                            
        end

        function a = area_of_intersection(obj, another_patch)            
            a = single(rectint(obj.as_rectangle, another_patch.as_rectangle));
%             [xmin1, ymin1, xmax1, ymax1] = obj.get_corners();
%             [xmin2, ymin2, xmax2, ymax2] = another_patch.get_corners();            
%             xmin = max(xmin1, xmin2);
%             ymin = max(ymin1, ymin2);
%             xmax = min(xmax1, xmax2);
%             ymax = min(ymax1, ymax2);
%             if ymin <= ymax && xmin <= xmax
%                 a = single((ymax-ymin+1)) * single((xmax-xmin+1));
%             else
%                 a = 0;
%             end
        end
        
        function s = intersection_over_union(self, another_patch)
            intersection = self.area_of_intersection(another_patch);
            union = self.area() + another_patch.area() - intersection;
            s = double(intersection) / double(union);
            assert(s >= 0 && s <= 1);
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
        
        function [I] = indicator_on_image(self, image)
            % Creates a logical matrix with size equal to the dimensions of the input image, 
            % where ones are only set in pixel locations covered by the patch.            
            I = zeros(image.height, image.width);
            [xmin, ymin, xmax, ymax] = self.get_corners();            
            I(ymin:ymax, xmin:xmax) = 1;
            I = logical(I);
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
        
        function [P] = tightest_surrrounding_non_zeros(in_mask)
            % Given a two-dimensional array find the tightest box/patch that contains all its non-zero values.
            % The edges of the tighest box include the extremal non-zero values.
            [y, x] = find(in_mask);
            if isempty(x) || isempty(y)
                P = Patch();
            else                
                P = Patch([min(x), min(y), max(x), max(y)]);
            end
        end
        
%         function area_percent_covered_by_other_patches(frame_patch, other_patches, corners)
%                 frame_corners = corners(frame_patch, :);
%                 mask = zeros(frame_corners(4) - frame_corners(2)+1, frame_corners(3)-frame_corners(1)+1);
%                 for p = 1:length(other_patches)
%                     sp = other_patches(p);                
%                     sp_corners = corners(sp,:);
%                     xmin = sp_corners(1) - frame_corners(1) + 1;
%                     xmax = xmin + s_corners(3) - sp_corners(1);
%                     ymin = sp_corners(2) - frame_corners(2) + 1;
%                     ymax = ymin + sp_corners(4) - sp_corners(2);
%                     mask(ymin:ymax, xmin:xmax) = 1;                                               
%                 end
% 
%             end
               
        function [s, qd] = angle_displacement(xc, yc, x, y)
                % xc, yc: x,y coordinates of center (reference) patch.
                % x, y  : x,y coordinates of other (reference) patch.                               
                hypotenuse = norm([xc-x, yc-y]);                
                qd         = 0;
                if y < yc
                    if x > xc      % 1st quadrant
                        s  = asin(norm([yc-y, 0]) / hypotenuse) * Utils.radians_to_degrees;                        
                        qd = 1;
                    elseif x < xc  % 2nd qd
                        s = 90 + (90 - (asin(norm([yc-y, 0]) / hypotenuse) * Utils.radians_to_degrees));                        
                        qd = 2;
                    else
                        s =  90;
                    end                    
                elseif y > yc
                    if x < xc      % 3rd qd
                        s = 180 + (asin(norm([0, yc-y]) / hypotenuse) * Utils.radians_to_degrees);                        
                        qd = 3;
                    elseif x > xc  % 4rth qd                                     
                        s = 270 + ((asin(norm([xc-x,0]) / hypotenuse) * Utils.radians_to_degrees));                        
                        qd = 4;
                    else
                        s = 270;
                    end
                else               % on xx'.
                    if x < xc  
                        s = 180;
                    else
                        s = 0;                        
                    end
                end    
%                 assert(s>=0 && s<=360)
        end
        
        
                 
    end
    
end