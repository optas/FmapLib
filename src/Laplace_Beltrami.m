classdef Laplace_Beltrami < dynamicprops
    
    properties (GetAccess = public, SetAccess = private)
        W           = [];               % Weight matrix of cotangent Laplacian.
        M           = [];               % Associated Mesh of LB.   
        spectra     = containers.Map;   % A dictionary carrying various types of eigenvalues and eigenvectors
    end
    
    methods
    
        function obj = Laplace_Beltrami(in_mesh)
            % Class constructor.
            if nargin == 0                                                
                obj.M = [];          
                obj.W = [];
                obj.spectra = containers.Map;                
            else
                obj.M       = in_mesh;
                obj.spectra = containers.Map;
                if isprop(in_mesh, 'angles')
                    obj.W  = Laplace_Beltrami.cotangent_laplacian(in_mesh.vertices, in_mesh.triangles, in_mesh.angles);
                else
                    obj.W  = Laplace_Beltrami.cotangent_laplacian(in_mesh.vertices, in_mesh.triangles);
                end                       
            end                                                                            
        end
                
        function [evals, evecs] = get_spectra(obj, eigs_num, area_type)
            if ~ Mesh.is_supported_area_type(area_type)
                error('You specidied an area_type which is not supported by the Mesh Library.')
            end
            
            
            if ~ obj.spectra.isKey(area_type) || ...                  % This type of spectra has not been computed before, or,
                size(obj.spectra(area_type).evals, 1) < eigs_num      % the requested number of eigenvalues is larger than what has been previously calculated.                
                
                try % Retrieve the vertex areas or compute them.
                    A = obj.M.get_vertex_areas(area_type);
                catch                                       
                    obj.M.set_vertex_areas(area_type);
                    A = obj.M.get_vertex_areas(area_type);
                end                
                A              = spdiags(A, 0, length(A), length(A));
                [evecs, evals] = Laplace_Beltrami.compute_spectra(obj.W, A, eigs_num);                
                obj.spectra(area_type) = struct('evals', evals, 'evecs', evecs); % Store computed spectra.
           
            else                                        % Use previously computed spectra.
                evals = obj.spectra(area_type).evals;
                evals = evals(1:eigs_num);                
                evecs = obj.spectra(area_type).evecs;
                evecs = evecs(:, 1:eigs_num);
            end
        end
        
        function [Proj] = project_functions(obj, area_type, eigs_num, varargin)
            n_varargin = nargin -3; % Number of arguments passed through varargin.
            
            if n_varargin < 1
                error ('Please provide some functions to be projected on the LB basis.');
            end
            
            num_vertices = size(obj.W, 1);
            
            functions_total = 0;            
            for i=1:n_varargin
                if size(varargin{i}, 1) ~= num_vertices                    
                    error ('Wrong dimensionality. The functions must be defined over a Mesh with a number of vertices equal to that of this LB.')
                end
                functions_total = functions_total + size(varargin{i}, 2);
            end

            [~, evecs] = obj.get_spectra(eigs_num, area_type);
            assert(num_vertices  == size(evecs, 1));
                            
            % Project feauture vectors into reduced LB basis.            
            Proj = zeros(functions_total, eigs_num);            % Pre-allocate space.
            right = 0;         
            for i = 1:n_varargin
                left  = right + 1;
                right = right + size(varargin{i}, 2);
                Proj(left:right, :) = evecs(:, 1:eigs_num)' * varargin{i};
            end            
        end

    end
    
    methods (Static)
        
        function [W] = cotangent_laplacian(V, T, varargin)
                % Computes teh cotangent laplacian weights.
                % W is symmetric.
                % optional third argument is the angles of the triangles of
                % the mesh, if not provided it will be calculated.
                I = [T(:,1); T(:,2); T(:,3)];
                J = [T(:,2); T(:,3); T(:,1)];        
                                
                if nargin == 3
                    A = varargin{1};
                elseif nargin == 4
                    A = Mesh.angles_of_triangles(V, T);
                else
                    error('Too many arguments were given.')
                end
                
                S = 0.5 * cot([A(:,3); A(:,1); A(:,2)]);
                In = [I; J; I; J];
                Jn = [J; I; I; J];
                Sn = [-S; -S; S; S];
                
                nv = size(V, 1);
                W  = sparse(In, Jn, Sn, nv, nv);
                assert(isequal(W, W'))
        end
        
        function [Phi, lambda] = compute_spectra(W, vertex_areas, eigs_num)
            % Returns the sorted ..add comments..
            % TODO-P: complex spectra
            if eigs_num < 1 || eigs_num > size(W, 1)-1;
                error('Eigenvalues must be in range of [1, num_of_vertices-1].')
            end
            
            [Phi, lambda] = eigs(W, vertex_areas, eigs_num, -1e-5);
            lambda        = diag(lambda);
            lambda        = abs(real(lambda));
            [lambda, idx] = sort(lambda);
            Phi           = Phi(:,idx);
            Phi           = real(Phi);                 % W is symmetric. diag(Vertex_Areas) is PSD. Thus, the Generalized Eigen-Prob should return real.
        end
            
    end
    
end

