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
                % Construct an empty Mesh_Collection.            
                obj.name = '';
                obj.meshes = containers.Map;                
            elseif nargin > 1   % At least a collection name and a top_directory were given.
                
                % Start by finding all potential objects (.off, .obj
                % files). % TODO-add optional -.mat format.
                all_subfiles = rdir([top_directory, '/**/'], 'regexp(name, ''\.obj$|\.off$'')');                                                
                num_meshes   = length(all_subfiles);
                if num_meshes == 0
                    warning(['The given top directory does not contain any .off or .obj files' ...
                             'in it or in any of its sub directories.']);
                    return
                end                
                obj.meshes = containers.Map;
                                
                for i=1:num_meshes
                    try                          
                        full_path = all_subfiles(i).name;
                        path_substrings = strsplit(full_path, '/');   % TODO: use separator of current system.
                        last_word = path_substrings{end};
                        mesh_name = last_word(1:end-4);    % Relying on the fact that length('.off') == length('.obj') == 4.                                                               
                        if obj.meshes.isKey(mesh_name)
                            error('Todo: meshes with same meshes are not supported at the moment.')
                        end                        
                        obj.meshes(mesh_name) = Mesh(full_path, mesh_name);
                    catch                         
                        warning([full_path, ' was not loaded.']);
                        continue
                    end
                end                            
            
%                 if exist('mesh_classes', 'var')
%                     % If two columns: ID - Class_Type (both strings or
%                     % string-int).
%                     
%                     % If more than two columms, then ID - Class_type_1, Class_type_k 
% 
%                     
%                 end
                
                obj.name = name;                
            end                
        end
        
%         function set_mesh_attributes(attribute_file)
%             
%         end
    


    end
    
    
    
    
    
end
