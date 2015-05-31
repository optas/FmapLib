classdef Mesh_IO
    % Implementing a set of functions facilitating Input - Output operations
    % on Triangular Meshes.
       
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

    end
end

