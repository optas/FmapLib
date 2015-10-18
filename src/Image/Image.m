classdef Image < dynamicprops
    % A class representing an arbitrary (rectangular) Image. This class plays mostly an organisational role by allowing
    % the user to keep data related to an image in a compact way.    
    % (c) Achlioptas, Corman, Guibas  - 2015  -  http://www.fmaplib.org
    
    properties (GetAccess = public, SetAccess = private)
        % Basic properties that every instance of the Image class has.
        
        CData        %   (height x width) or (height x width x 3) matrix capturing the image's color data. This
                     %   property has an identical role as the CData property of the 'matlab.graphics.primitive.Image' 
                     %   class.
                     %   
        height       %   (int)     -    Number of vertical pixels.
        width        %   (int)     -    Number of horizontal pixels.
        name         %   (String)  -    (default = '') A string identifying the image, e.g., 'Mandrill'.
    end
    
    methods (Access = public)
        % Class Constructor.
        % Input:
        %       First argument:
        %           Case 1. Filepath (string) of an image, which will be opened by imread.
        %           Case 2. 2D or 3D Matrix containing the image color data.
        %
        %       Second argument:
        %           (optional, string) describing the name of the image.        
        
        function obj = Image(varargin)
            if nargin == 0  
                obj.CData = [];                
            elseif ischar(varargin{1}) % TODO-P Specialize on diff types. && strcmp(varargin{1}(end-3:end), '.obj')
                if size(varargin{1}, 1) < 1 || size(varargin{1}, 2) < 1
                    error('Image constructor: stores image matrices and expects at least a 2D matrix as input.')
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
        
        function [varargout] = size(obj)                        
            if nargout == 2
                varargout = cell(nargout);
                varargout{1} = obj.width;
                varargout{2} = obj.height;
            else
                varargout{1} = [obj.width, obj.height];
            end
        end
        
        function [h] = plot(obj)
            % Plots the image and returns the graphic's handle to it.            
            h = imshow(obj.CData);
            if ~ isempty(obj.name)
                title(['Image name = ' obj.name]);
            end
        end
        
        function c = color(obj)
            % Returns the every pixel with its content (i.e., color).
            if isa(obj.CData, 'uint8') || isa(obj.CData, 'uint16')
                c = im2double(obj.CData);                  % Makes each chanel have values in [0,1].
            else 
                obj.CData;    % TODO-P See what other data types are expected on an image.
            end
        end
           
        function a = area(obj)                       
            a = obj.width * obj.height;
        end
        
        function [b] = is_rgb(obj)
            b = ndims(obj.CData) == 3;
        end          
                         
        function obj = set_patches(obj, patches)
            propname = 'patches';
            if isprop(obj, propname)
                warning('Updating image to a new set of patches.');                
            else
                obj.addprop(propname);
            end            
            obj.(propname) = Patch_Collection(patches, obj);
        end
                
        function [h] = plot_patch(obj, patch)
            h = obj.plot();
            hold on;
            [xmin, ymin, xmax, ymax] = patch.get_corners();            
            plot([xmin xmax xmax xmin xmin],[ymin ymin ymax ymax ymin], 'Color', 'r', 'LineWidth', 2);
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
        
        function [F] = content_in_patch(obj, patch)
            [xmin, ymin, xmax, ymax] = patch.get_corners();            
            F = obj.CData(ymin:ymax, xmin:xmax, :);
        end
        
        function [new_im] = resize(obj, new_height, new_width)
            imres = imresize(obj.CData , [new_height, new_width], 'bilinear');
            new_im = Image(imres, [obj.name '-resized']);
        end
        
        function set_resized_image(obj, new_height, new_width)
            % to do change to varargin
            propname = 'resized';            
            if ~ isprop(obj, propname)
                obj.addprop(propname);            
            end
            obj.(propname) = obj.resize(new_height, new_width);
        end
        
        function [resized] = get_resized_image(obj)
            resized = obj.resized;
        end             
    end % End of (public) instance methods.
    
end