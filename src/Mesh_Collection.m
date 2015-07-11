classdef Mesh_Collection < dynamicprops
    % A class offering a variety of utilities for collecting, maintaining and experimenting with collections of 
    % Triangular Meshes.
        
    properties               % Each Mesh_Collection object has at least the following properties.        
        name                 % (String)            -   (default = '') A string identifying the collection i.e., 'Tosca'.    
        meshes               % (Containers.Map)    -   A dictionary storing as values the meshes of the collection. The 
                             %                         keys are strings corresponding to mesh names.
    end
    
    methods
        
        function obj = Mesh_Collection(collection_name, top_directory, attributes)                                        
            % Class Constructor.
            % Input:
            %             
            %           collection_name  -  (String) Name given to colleciton of meshes (e.g., TOSCA)
            %             
            %           top_directory    -  (String) Filepath of the top-directory that includes (e.g., in sub-directories)
            %                               the meshes of the collection.
            %              
            %           attributes       -  (String, optional) File containing attribute-like information of the meshes, see Notes.
            %
            % Notes: TODO-P describe format of attributes file
            
            if nargin == 0                
                % Construct an empty Mesh_Collection.
                obj.name = '';
                obj.meshes = containers.Map;                            
            elseif nargin == 2 || nargin == 3                
                % Find all potential mesh files (.off, .obj files).
                all_subfiles = rdir([top_directory, '/**/'], 'regexp(name, ''\.obj$|\.off$'')');                                                
                num_meshes   = length(all_subfiles);
                if num_meshes == 0
                    warning(['The given top directory does not contain any .off or .obj files' ...
                             'in it or in any of its sub directories.']);
                    return
                end
                
                obj.meshes = containers.Map;                                
                for i=1:num_meshes   % Load and store the meshes.
                    full_path = all_subfiles(i).name;
                    try                          
                        mesh_name = extract_mesh_name(full_path);
                        if obj.meshes.isKey(mesh_name)
                            error('Sorry, but currently every mesh is assumed to have a unique filename.')
                        end                        
                        obj.meshes(mesh_name) = Mesh(full_path, mesh_name);
                        disp([full_path, ' was loaded.']);
                    catch                         
                        warning([full_path, ' was not loaded.']);
                        continue
                    end
                end                                        
                if exist('attributes', 'var')
                    error('Not implemented yet.')                    
                end                
                obj.name = collection_name;                                
            end
           
            
            function [mesh_name] = extract_mesh_name(full_path)                
                path_substrings = strsplit(full_path, '/');   % TODO: use separator of current system.
                last_word = path_substrings{end};
                mesh_name = last_word(1:end-4);               % Relying on the fact that length('.off') == length('.obj') == 4.                                                               
            end

                     
%             function set_attributes_of_mesh_collection(mesh_collection, attributes)            % TODO-P Move to Mesh_IO.
%                 C               = textread(attributes, '%s', 'delimiter', '\n');                
% %               C               = cellfun(Mesh_IO.string_or_num, C, 'uni', false);
%                 
%                 attribute_types = strsplit(C{1});
%                 
%                 for line=2:length(C)
%                     mesh_name = C{line}{1};
%                     if mesh_collection.meshes.isKey(mesh_name)                        
%                         mesh     = mesh_collection.meshes(mesh_name);
%                         content  = strsplit(C{line});                   % Attributes of mesh corresponding to line.                        
%                         content  = content(2:end);                                                     
%                         mesh.set_semantics(attribute_types, )
%                     end
%                     
%                 end
% 
%             
%             end
            
        end
        
        %         function B = subsref(obj, S)
        %             if strcmp(S(1).type, '()')                
        %                 B = obj.meshes(S.subs{:});
        %             end
        %         end
        
        function [C] = get_property_of_meshes(obj, property_name, mesh_list)                
                if ~ isprop(obj, property_name)                    
                    error('Property requested does not exist.');
                end
                
                C = cell(length(mesh_list), 1);
                
                for i = 1:length(mesh_list)
                    meshname = mesh_list{i};
                    C{i} = obj.(property_name)(meshname);
                end
       end
                
        function compute_laplace_beltrami_basis(obj, num_eigs, area_type)                        
            if ~ isprop(obj, 'lb_basis')
                obj.addprop('lb_basis');
                obj.lb_basis = containers.Map;            
            end
                                    
            for key = obj.meshes.keys               
                meshname = key{:};
                m    = obj.meshes(meshname);
                if exist('area_type', 'var')
                    obj.lb_basis(meshname) = Laplace_Beltrami(m, Mesh.area_of_vertices(m.vertices, m.triangles, area_type));
                else
                    obj.lb_basis(meshname) = Laplace_Beltrami(m);
                end
                obj.lb_basis(meshname).get_spectra(num_eigs);                
                disp(['Computing Laplace Beltrami basis for: ', meshname, ' done.']);
            end
        end
                
        function compute_default_feautures(obj)                        
            if ~ isprop(obj, 'raw_features')
                obj.addprop('raw_features');
                obj.raw_features = containers.Map;            
            end
                                             
            for key = obj.meshes.keys               
                meshname = key{:};               
                m     = obj.meshes(meshname);
                lb    = obj.lb_basis(meshname);
                neigs = length(lb.spectra.evals);                
                F     = Mesh_Features.default_mesh_feauture(m, lb, neigs);
                obj.raw_features(meshname) = F;
            end
        end
                
        function [C] = project_features(obj, mesh_list, eigs_num, features)
                C = cell(length(mesh_list), 1);
                
                for i = 1:length(mesh_list)
                    meshname = mesh_list{i};                    
                    
                    if strcmp(eigs_num, 'all')
                        eigs_num = length(obj.lb_basis(meshname).spectra.evals);
                    end

                    C{i} = obj.lb_basis(meshname).project_functions(eigs_num, features{i});
                    normalize = 0; % TODO-P
                    if normalize
                        C{i} = divide_columns(C{i}, l2_norm(C{i}'));
                    end
                end

        end
        

        function [fmaps] = compute_ground_truth_fmaps(obj, pairs, groundtruth)                        
            num_pairs = size(pairs, 1);
            fmaps     = cell(num_pairs, 1);
           
            for i = 1:num_pairs
                meshname_source = pairs{i,1};
                meshname_target = pairs{i,2};            
                
                lb_src   = obj.lb_basis(meshname_source);
                lb_tar   = obj.lb_basis(meshname_target);
                
                fmaps{i} = Functional_Map.groundtruth_functional_map(lb_src.spectra.evecs, lb_tar.spectra.evecs, groundtruth{i}, lb_tar.A);  
            end
        end
        
       
        
    end
 
end
