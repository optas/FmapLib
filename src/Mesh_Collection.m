classdef Mesh_Collection
    % A class offering a variety of utilities in order to collect,
    % maintain and experiment with collecitons of Triangular Meshes.
    
    properties
        name         %   (string)    -    (default = '') a string identifying the collection i.e., 'Tosca'.    
        meshes            
    end
    
    methods
        function obj = Mesh_Collection(collection_name, top_directory, attributes_file)                            
            if nargin == 0                
                % Construct an empty Mesh_Collection.            
                obj.name = '';
                obj.meshes = containers.Map;                
            elseif nargin > 1   % At least a collection name and a top_directory were given.
                
                % Start by finding all potential objects (.off, .obj files). % TODO-add optional -.mat format.                
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
                        mesh_name = extract_mesh_name(all_subfiles(i).name);
                        if obj.meshes.isKey(mesh_name)
                            error('Sorry, but currently every mesh is assumed to have a unique filename.')
                        end                        
                        obj.meshes(mesh_name) = Mesh(full_path, mesh_name);
                    catch                         
                        warning([full_path, ' was not loaded.']);
                        continue
                    end
                end                            
            
                if exist('attributes_file', 'var')
                    % ID - Attribute_1, Attribute_i

  
                end                
                obj.name = collection_name;                                
            end
           
            
            function [mesh_name] = extract_mesh_name(full_path)                
                path_substrings = strsplit(full_path, '/');   % TODO: use separator of current system.
                last_word = path_substrings{end};
                mesh_name = last_word(1:end-4);               % Relying on the fact that length('.off') == length('.obj') == 4.                                                               
            end
            
            
            
                       
            function set_attributes_of_mesh_collection(mesh_collection, attribute_file)            % Move to Mesh_IO
                
                C               = textread(attribute_file, '%s', 'delimiter', '\n');                
%                 C               = cellfun(Mesh_IO.string_or_num, C, 'uni', false);
                
                attribute_types = strsplit(C{1});
                
                for line=2:length(C)
                    mesh_name = C{line}{1};
                    if mesh_collection.meshes.isKey(mesh_name)                        
                        mesh     = mesh_collection.meshes(mesh_name);
                        content  = strsplit(C{line});                   % Attributes of mesh corresponding to line.                        
                        content  = content(2:end);     
                                                
                        mesh.set_semantics(attribute_types, )
                    end
                    
                end

            
            end
            
        end

    end
    
    
    
    
    
end
