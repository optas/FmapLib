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
            elseif nargin == 2   % A collection name and a top_directory were given.                 
                % All potential objects will be loaded (.off, .obj).
                all_subfiles = rdir([top_directory, '/**/'], 'regexp(name, ''\.obj$|\.off$'')');                                
                
                num_meshes = length(all_subfiles);
                if num_meshes == 0
                    warning(['The given top directory does not contain any .off or .obj files' ...
                             'in it or in any of its sub directories.']);
                    return
                end                
                obj.meshes = containers.Map;                
                for i=1:num_meshes
                    try                          
                        full_path = all_subfiles(i).name;
                        path_substrings = strsplit(full_path, '/');
                        last_word = path_substrings{end};
                        mesh_name = last_word(1:end-4);    % Relying on the fact that length('.off') == length('.obj') == 4.                                                               
                        if obj.meshes.isKey(mesh_name)
                            error('todo;explain.')
                        end                        
                        obj.meshes(mesh_name) = Mesh(full_path, mesh_name);                        
                    catch                         
                        warning([full_path, ' was not loaded.']);
                        continue
                    end
                end
                            
                obj.name = name;
            end
                
                
                
        end

        
    end

    
end
