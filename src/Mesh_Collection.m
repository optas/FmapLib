classdef Mesh_Collection
    % A class offering a variety of utilities in order to collect,
    % maintain and experiment with collecitons of Triangular Meshes.
    
    properties
        name         %   (string)    -    (default = '') a string identifying the collection i.e., 'Tosca'.    
        meshes            
    end
    

    methods

        function obj = Mesh_Collection(name, top_directory, mesh_classes)     
            
            
            if nargin == 0                
                % Construct an empty Mesh.            
                obj.name = '';
                obj.meshes = [];                
            elseif nargin == 1   % Only a directory was given: All potential objects will be loaded (.off, .obj)
            
                imgPath = '../data/jaffe/';
imgType = '*.tiff';
imageFiles = dir([imgPath imgType]);
images     = cell(length(imageFiles), 1);
for idx = 1:length(imageFiles)
    [img, cmap] = imread([imgPath imageFiles(idx).name]);
    images{idx} = struct('raw_image', img, 'name', imageFiles(idx).name);     
end
                
                
                
            end
                
%                 ischar(varargin{1}) && strcmp(varargin{1}(end-3:end), '.off')
%                 % Constructor from input .off file.
%                 off_file = varargin{1};                
%                 [obj.vertices , obj.triangles] = Mesh_IO.read_off(off_file); 
%             elseif ischar(varargin{1}) && strcmp(varargin{1}(end-3:end), '.obj')
%                 % Constructor from input .obj file.
%                 obj_file = varargin{1};                
%                 [obj.vertices , obj.triangles] = Mesh_IO.read_obj(obj_file); 
%             else
%                 % Construct with explicitly given vertices/triangles.
%                 obj.vertices  = varargin{1};
%                 obj.triangles = varargin{2};
%             end
%             
%             % Take care of sizes and name.
%             obj.num_vertices  = size(obj.vertices, 1);
%             obj.num_triangles = size(obj.triangles, 1);                           
%             
%             if nargin > 1 && ischar(varargin{end}) % The last argument is a string (and this string is not the first vararg which is reserved for filenames).
%                 obj.name = varargin{end};
%             else
%                 obj.name = '';
%             end
%         end
     
        
        
    end

    
    
    
    
end



fprintf('Computing features\n');

if strcmp(exp_type,'train') %training features
	%Read classes and models
	dir_class = dir(opts3DSCB.train_models_dir);
	models_dir = opts3DSCB.train_models_dir;
elseif strcmp(exp_type,'test')
	%Read classes and models
	dir_class = dir(opts3DSCB.test_models_dir);
	models_dir = opts3DSCB.test_models_dir;
end

dir_class = dir_class(3:end,1);
%dir_class contains the name of all the classes

%compute number of classes
numc = length(dir_class);

%read models per class
for n = 1:numc
    %build class dir name
    class_dir_name = [models_dir dir_class(n).name];
    
    name_class = dir_class(n).name;
    
    %read models (.obj)
    models = dir([class_dir_name '/*.obj']);
    nmodels = length(models);
