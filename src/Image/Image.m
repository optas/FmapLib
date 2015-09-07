classdef Image < dynamicprops
    % A class representing an arbitrary (rectangular) Image. This class plays mostly an organisational role by allowing
    % the user to keep data related to an image in a compact way.    
    % (c) Achlioptas, Corman, Guibas  - 2015  -  http://www.fmaplib.org
    
    properties (GetAccess = public, SetAccess = private)
        % Basic properties that every instance of the Image class has.
        
        CData          %   (height x weight) or (height x weight x 3) matrix capturing the image's color data. This
                       %   property has an identical role as the CData property of the 'matlab.graphics.primitive.Image'
                       %   class.
        height         %   (int)     -    Number of vertical pixels.
        weight         %   (int)     -    Number of horizontal pixels.
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
            [obj.height, obj.weight, ~] = size(obj.CData);
            if nargin > 1 && ischar(varargin{end}) % Add potential name of picture.
                obj.name = varargin{end};
            else
                obj.name = '';
            end                           
        end
        
        function [h] = plot(obj)
            h = image(obj.CData);
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
        
    end
    
end