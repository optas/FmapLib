classdef Image < dynamicprops
    % A class representing an arbitrary (rectangular) Image. This class plays mostly an organisational role by allowing
    % the user to keep data related to an image in a compact way.    
    % (c) Achlioptas, Corman, Guibas  - 2015  -  http://www.fmaplib.org
    
    properties (GetAccess = public, SetAccess = private)
        % Basic properties that every instance of the Image class has.
        
        CData          %   (height x width) or (height x width x 3) matrix capturing the image's color data. This
                       %   property has an identical role as the CData property of the 'matlab.graphics.primitive.Image'
                       %   class.
        height         %   (int)     -    Number of vertical pixels.
        width          %   (int)     -    Number of horizontal pixels.
        name           %   (String)  -    (default = '') A string identifying the image, e.g., 'Mandrill'.
    end
    
    methods (Access = public)
        % Class Constructor.               
        function obj = Image(varargin)     
            if nargin == 0  
                obj.CData = [];                
            elseif ischar(varargin{1}) % TODO-P Specialize on diff types. && strcmp(varargin{1}(end-3:end), '.obj')                
                if size(varargin{1}, 1) < 1 || size(varargin{1}, 2) < 1
                    error('Image constructor, stores image matrices and expects at least a 2D matrix as input.')
                end
                obj.CData = imread(varargin{1});
            else % Directly provide the matrix with the pixel conntent.
                obj.CData = varargin{1};                
            end
            [obj.height, obj.width, ~] = size(obj.CData);
            if nargin > 1 && ischar(varargin{end}) % Add potential name of picture.
                obj.name = varargin{end};
            else
                obj.name = '';
            end                           
        end
        
        function [s] = size(obj)
            s = size(obj.CData);
        end
        
        function c = color(obj)
            if isa(obj.CData, 'uint8') || isa(obj.CData, 'uint16')
                c = im2double(obj.CData);
            else 
                obj.CData;    % More on this case.
            end
        end
        
        function [h] = plot(obj)
            h = image(obj.CData);
        end
        
        function [h] = plot_patch(obj, patch, varargin)
            xmin = patch(1);
            ymin = patch(2);
            xmax = patch(3);
            ymax = patch(4);            
            h = obj.plot();
            hold on;
            plot([xmin xmax xmax xmin xmin],[ymin ymin ymax ymax ymin], ...
                 'Color', 'r', 'LineWidth', 3); %, 'LineStyle',lineStyle);
        end
        
        
        function [b] = is_rgb(obj)
            b = ndims(obj.CData) == 3;
        end        
        
        function [obj] = set_gt_segmentation(obj, segmentation)
            % Adds dynamic property 'gt_segmentation' corresponding to a groundtruth segmentation of the image.
            
            propname = 'gt_segmentation';
            if isprop(obj, propname)
                obj.(propname) = segmentation;
            else
                obj.addprop(propname);
                obj.(propname) = segmentation;
            end
        end
        
        function [F] = content_in_rectangle(obj, xmin, ymin, xmax, ymax )
            if nargin == 2                
                ymin = xmin(2);
                xmax = xmin(3);
                ymax = xmin(4);
                xmin = xmin(1);
            end
            % Add checks
            F = obj.CData(ymin:ymax, xmin:xmax, :);
        end
        
        function set_resized_image(obj, new_height, new_width)
            % to do change to varargin
            propname = 'resized';
            imres = imresize(obj.CData , [new_height, new_width], 'bilinear');
            if isprop(obj, propname)
                obj.(propname) = Image(imres);
            else
                obj.addprop(propname);
                obj.(propname) = Image(imres);
            end
        end
        
        function [resized] = get_resized_image(obj)
            resized = obj.resized;
        end
       
        function [obj] = set_patches(obj, patches)                        
            propname = 'patches';
            if isprop(obj, propname)
                obj.(propname) = patches;
            else
                obj.addprop(propname);
                obj.(propname) = patches;
            end
        end
        
        
        
        
        
 
    end
    
end