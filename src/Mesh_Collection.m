classdef Mesh_Collection < dynamicprops
    % A class offering a variety of utilities for collecting, maintaining and experimenting with collections of 
    % Triangular Meshes.
        
    properties               % Each Mesh_Collection object has at least the following properties.        
        name                 % (String)            -   (default = '') A string identifying the collection i.e., 'Tosca'.    
        meshes               % (Containers.Map)    -   A dictionary storing as values the meshes of the collection. The 
                             %                         keys are strings corresponding to mesh names.
    end
    
    methods
        
        function obj = Mesh_Collection(collection_name, top_directory, attributes_file)                                        
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
            elseif nargin == 1
                error('Wrong input arguements.')
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
                if exist('attributes_file', 'var')
                    error('Not implemented yet.')                    
                end                
                obj.name = collection_name;                                
            end
           
            
            
            function [mesh_name] = extract_mesh_name(full_path)                
                path_substrings = strsplit(full_path, '/');   % TODO: use separator of current system.
                last_word = path_substrings{end};
                mesh_name = last_word(1:end-4);               % Relying on the fact that length('.off') == length('.obj') == 4.                                                               
            end

                                 
        end
        
        function [bool] = contains(obj, mesh_name)
            % Returns true iff the collection contains a mesh with the given mesh_name.
            bool = obj.meshes.isKey(mesh_name);
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
        
       
    
        function [D, i] = set_semantics_of_meshes(obj, semantics)
            % TODO-P Deal with missing data
            fid = fopen(semantics);
            C   = textscan(fid, '%s', 'delimiter', '\n');
            C   = C{1};
            
            i = 1; % Skip comments and empty lines.
            while(all(isstrprop(C{i}, 'wspace')) || C{i}(1)=='#')
                i = i+1;
            end
                       
            A = cellfun(@strsplit, C(i:end, :), 'UniformOutput', false);
            B = cell(1);      
            m = 1;
            for i=1:length(A)                
                if i==1  || (obj.contains(A{i}{1}))   % Keep info for stored meshes (i=1 case to get semantics header).                
                    for j=1:length(A{i})
                        B{m,j} = Mesh_IO.string_or_num(A{i}{j});
                    end
                    m = m+1;
                end
            end
        
            D  = cell2dataset(B, 'readObsName', true);

            if ~ isprop(obj, 'mesh_semantics')
                obj.addprop('mesh_semantics');
                obj.mesh_semantics = D;                
            end

  
%             
%             
%             semantic_types = strsplit(C{i});
%             
%             
%             
%             
%             semantic_types = semantic_types (2:end); % First column correspond to name, is ignored.
% 
%             obj.mesh_semantics.types = semantic_types;
%             
%             for line=i:length(C)
%                 content = strsplit(C{line});
%                 if length(content) ~= length(attributes) + 1
%                     error('Semantics file format. Use ''-'' for missing values.')
%                 end
%                 mesh_name = content{1};
%                 if obj.contains(mesh_name)                        
%                    content(2:end)
% 
%                    obj.mesh_semantics.C
%                    (mesh_name, attribute_types, ); 
%                 end
% 
%             end
        end
        
        
    end
    
   
 
end
