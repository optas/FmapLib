classdef Mesh_Collection < dynamicprops
    % A class offering a variety of utilities for collecting, maintaining and experimenting with collections of 
    % Triangular Meshes.
    %
    % (c) Achlioptas, Corman, Guibas  - 2015  -  http://www.fmaplib.org
        
    properties               % Each Mesh_Collection object has at least the following properties.        
        name                 % (String)            -   (default = '') A string identifying the collection i.e., 'Tosca'.    

        meshes               % (Containers.Map)    -   A dictionary storing as values the meshes of the collection. The 
                             %                         keys are strings corresponding to mesh names.
    end
    
    methods
        
        function obj = Mesh_Collection(collection_name, top_directory, semantics_file)                                        
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
                error('Wrong input arguments.')
            elseif nargin == 2 || nargin == 3                
                % Find all potential mesh files (.off, .obj files).                
                all_subfiles = rdir([top_directory, [filesep '**' filesep]], 'regexp(name, ''\.obj$|\.off$'')');                
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
                if exist('semantics_file', 'var')
                    obj.set_semantics_of_meshes(semantics_file);
                end                
                obj.name = collection_name;                                
            end
                     
            function [mesh_name] = extract_mesh_name(full_path)                
                path_substrings = strsplit(full_path, filesep); 
                last_word = path_substrings{end};
                mesh_name = last_word(1:end-4);               % Relying on the fact that length('.off') == length('.obj') == 4.                                                               
            end                                 
        end
        
        function compute_laplace_beltrami_basis(obj, num_eigs, area_type, mesh_list)
            % Computes and stores the laplace beltrami (LB) basis of the meshes kept by the Mesh_Collection.
            % Input:
            %        num_eigs   -  (int)            
            %        
            %        mesh_list  -  (optional, cell array of strings) specifying the names of the meshes for which
            %                       the LB basis is computed. If ommited, the basis of all the meshes kept by the
            %                       Mesh_Collection are computed.
            
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
        
        function [fmaps] = compute_fmaps(obj, pairs, features, method, varargin)
            
            % features  - containers.map (String to Mesh_Features), String corresponds to the name of a mesh.
            
            options = struct('eigs', 'all', 'lambda', 0);
            options = load_key_value_input_pairs(options, varargin{:});
            
            num_pairs = size(pairs, 1);
            fmaps     = containers.Map;
           
            for i = 1:num_pairs
                src_name = pairs{i,1}; 
                trg_name = pairs{i,2};                
                lb_src   = obj.lb_basis(src_name);
                lb_tar   = obj.lb_basis(trg_name);                
                                
                if ~ fmaps.isKey(src_name)                          % Initialize a dictionary per source mesh.
                    fmaps(src_name) = containers.Map;                   
                end                
                fmap_src_dic           = fmaps(src_name);                                            
                fmap_src_dic(trg_name) = Functional_Map(lb_src, lb_tar);                                
                if strcmp(options.eigs, 'all') 
                    neigs_source = length(lb_src.spectra.evals);        
                    neigs_target = length(lb_tar.spectra.evals);        % TODO-P, we re-project the raw features.
                else
                    neigs_source = options.eigs;
                    neigs_target = options.eigs;
                end
                fmap_src_dic(trg_name).compute_f_map(method, neigs_source, neigs_target, features(src_name), features(trg_name), 'lambda', options.lambda);
            end
        end
        

        function [bool] = contains(obj, mesh_name)
            % Returns true iff the collection contains a mesh with the given mesh_name.
            bool = obj.meshes.isKey(mesh_name);
        end
        
        function s = size(obj)
            s = size(obj.meshes, 1);
        end
        
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
                   
        function compute_default_feautures(obj, wks_samples, hks_samples, mc_samples, gc_samples)                        
            if ~ isprop(obj, 'raw_features')
                obj.addprop('raw_features');
                obj.raw_features = containers.Map;            
            end
            neigs = 'all';
            for key = obj.meshes.keys               
                meshname = key{:};                                               
                feats = Mesh_Features(obj.meshes(meshname), obj.lb_basis(meshname));                                            
                feats.compute_default_feautures(neigs, wks_samples, hks_samples, mc_samples, gc_samples);
                obj.raw_features(meshname) = feats;
                disp(['Computing Default raw Features for: ', meshname, ' done.']);
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
                    normalize = 1; % TODO-P
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
               
        function [mesh_names] = meshes_with_semantic_condition(obj, semantic_var, condition)            
            D = obj.mesh_semantics;
            if ischar(condition) % This imples that the semantic is a categorical variable.
                index = find(strcmp([D.(semantic_var)], condition));                                
            else
                index = D.(semantic_var) == condition;
            end
            mesh_names = D.Properties.ObsNames(index);
        end
                        
        function [D] = set_semantics_of_meshes(obj, semantics_file)
            % TODO-P Deal with missing data
            fid = fopen(semantics_file);            
            
            if fid < 0
                warning('Cannot open semtantic file. Collection will not save any semantic information.')                 
                return
            end
            
            C   = textscan(fid, '%s', 'delimiter', '\n');
            fclose(fid);
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
        end
        
        
    end
    
   
 
end
