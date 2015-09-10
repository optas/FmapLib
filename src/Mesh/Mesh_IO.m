classdef Mesh_IO
    % Implementing a set of functions facilitating Input and Output operations
    % on Triangular Meshes.
    %
    % (c) Achlioptas, Corman, Guibas  - 2015  -  http://www.fmaplib.org
       
    methods (Static)
        function [vertices, faces] = read_off(filename)
            %   Reads mesh data from a .off file.
            %
            %   Input:
            %           filename    -   (string) file name of the .off file.
            %
            %   Output:
            %           vertices    -   (number of vertices x 3) array specifying the 3D embedding of the vertices.
            %           faces       -   (number of faces x 3) array specifying the 3D embedding of the faces.
            %
            %   Copyright (c) 2003 Gabriel Peyre            
            %   Copyright (c) 2015 Panos Achlioptas http://www.stanford.edu/~optas

            fid = fopen(filename, 'r');
            if (fid == -1)
                error(['Cannot open the file: ', filename, '.']);        
            end

            str = fgets(fid);   % -1 if eof.
            if ~strcmp(str(1:3), 'OFF')
                error('This file does not compy with the .off prototype.');    
            end

            str = fgets(fid);
            sizes = sscanf(str, '%d %d', 2);
            while length(sizes) ~= 2
                str = fgets(fid);
                sizes = sscanf(str, '%d %d', 2);
            end

            nv = sizes(1);      % Number of Vertices.
            nf = sizes(2);      % Number of Faces.

            % Read vertices.
            [vertices, cnt] = fscanf(fid, '%lf %lf %lf\n', [3,nv]);
            if cnt ~= 3*nv
                warning('Unexpected error while reading vertices.');
            end
            vertices = vertices';

            % Read faces.
            [faces, cnt] = fscanf(fid,'3 %ld %ld %ld\n', [3,inf]);    
            if cnt ~= 3*nf
                warning('Unexpected error while reading faces.');
            end        
            faces = double(faces' + 1);

            fclose(fid);
        end
        
        
        function [vertices, faces, normal] = read_obj(filename)
            %   Reads mesh data from a .obj file.
            %
            %   Input:
            %           filename    -   (string) file name of the .obj file.
            %
            %   Output:
            %           vertices    -   (number of vertices x 3) array specifying the 3D embedding of the vertices.
            %           faces       -   (number of faces x 3) array specifying the 3D embedding of the faces.
            %           normal      -   TODO
            %
            %   Copyright (c) 2008 Gabriel Peyre
            %   Copyright (c) 2015 Panos Achlioptas http://www.stanford.edu/~optas

            fid = fopen(filename, 'r');
            if (fid == -1)
                error(['Cannot open the file: ', filename, '.'])
            end
            
            try % This is a fast way of reading the file but will not work on all .obj content. TODO: expand.
                C = textscan(fid, '%c %f %f %f');                
                vertex_index  = C{1} == 'v';
                vertices      = zeros(sum(vertex_index), 3);
                vertices(:,1) = C{2}(vertex_index);
                vertices(:,2) = C{3}(vertex_index);
                vertices(:,3) = C{4}(vertex_index);
                
                face_index    = C{1} == 'f';
                faces         = zeros(sum(face_index), 3);
                faces(:,1)    = int64(C{2}(face_index));        % Using 64 bit integers in ecessive given the standard sizes of meshes.
                faces(:,2)    = int64(C{3}(face_index));
                faces(:,3)    = int64(C{4}(face_index));
                normal        = [];
                fclose(fid);                                
            catch
                frewind(fid);
                a = fscanf(fid, '%c', 1);            
                if strcmp(a, 'P')
                % This is the montreal neurological institute (MNI) specific ASCII facesangular mesh data structure.
                % For FreeSurfer software, a slightly different data input coding is
                % needed. It will be provided upon request.
                    fscanf(fid,'%f',5);
                    n_points=fscanf(fid,'%i',1);
                    vertices=fscanf(fid,'%f',[3,n_points]);
                    normal=fscanf(fid,'%f',[3,n_points]);
                    n_faces=fscanf(fid,'%i',1);
                    fscanf(fid,'%i',5+n_faces);
                    faces=fscanf(fid,'%i',[3,n_faces])'+1;
                    fclose(fid);
                    return;
                end
                frewind(fid);
                vertices = [];
                faces = [];
                normal= [];
                while 1
                    s = fgetl(fid);
                    if ~ischar(s),
                        break;
                    end
                    if ~isempty(s) && strcmp(s(1), 'f')
                        % face
                        faces(:,end+1) = sscanf(s(3:end), '%d %d %d');
                    end
                    if ~isempty(s) && strcmp(s(1), 'v')
                        % vertex
                        vertices(:,end+1) = sscanf(s(3:end), '%f %f %f');
                    end
                end
                vertices = vertices';
                faces    = faces';
                fclose(fid);
            end            
        end
        
        function write_off(filename, vertex, face, renormalize)
            % write_off - write a mesh to a OFF file
            %
            %   write_off(filename, vertex, face);
            %
            %   vertex must be of size [n,3]
            %   face must be of size [p,3]
            %
            %   Copyright (c) 2003 Gabriel Peyré

            if nargin<4
                renormalize = 0;
            end
            if size(vertex,2)~=3
                vertex=vertex';
            end
            if size(vertex,2)~=3
                error('vertex does not have the correct format.');
            end

            if renormalize==1
                m = mean(vertex);
                s = std(vertex);
                for i=1:3
                    vertex(:,i) = (vertex(:,i)-m(i))/s(i);
                end
            end

            if size(face,2)~=3
                face=face';
            end
            if size(face,2)~=3
                error('face does not have the correct format.');
            end

            fid = fopen(filename,'wt');
            if( fid==-1 )
                error('Can''t open the file.');
                return;
            end

            % header
            fprintf(fid, 'OFF\n');
            fprintf(fid, '%d %d 0\n', size(vertex,1), size(face,1));

            % write the points & faces
            fprintf(fid, '%f %f %f\n', vertex');
            fprintf(fid, '3 %d %d %d\n', face'-1);

            fclose(fid);
        end
    
        function [converted] = string_or_num(x)
            [converted, status] = str2num(x);
            if ~ status
                converted = x;
            end
        end
        
        function mat_to_off(top_directory)
            % Converts .mat files containing meshes to .off.             
            % Assumes topdir and its subirectories include .mat files that correspond to triangular meshes. 
            % Also the .mat loads a structure with 4 fields: TRIV (triangles), X,Y,Z
            % (vertices). TODO-V parameterize expected .mat structure.
            % Input:
            %        top_directory  - (String) Filename of a directory which contains the .mat files directly under it
            %                         or in subdirectories.
            
            all_subfiles = rdir([top_directory, '/**/'], 'regexp(name, ''\.mat$'')');                                                
            num_meshes   = length(all_subfiles);
            if num_meshes == 0
                warning(['The given top directory does not contain any .mat files' ...
                             'in it or in any of its sub directories.']);
                return
            end
            
            for i=1:num_meshes   % Load meshes.
                full_path = all_subfiles(i).name;
                try                          
                    S    = load(full_path);  % The loaded stucture.
                    name = fieldnames(S);  % The name of the top-field of S.
                    
                    vertices = [S.(name{1}).X, S.(name{1}).Y, S.(name{1}).Z];                                                                               
                    file_name = strcat(full_path(1:end-3), 'off');
                    Mesh_IO.write_off(file_name, vertices, S.(name{1}).TRIV);
                    disp([full_path, ' was not loaded.']);
                catch                         
                    warning([full_path, ' was not loaded.']);
                    continue
                end
           end                                        

        end
                                
    end
    
    
    
end

