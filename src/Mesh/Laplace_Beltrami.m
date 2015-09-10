classdef Laplace_Beltrami < Basis
    % A class representing the cotangent discretization of the Laplace Beltrami operator, associated with a given
    % object of the class Mesh.
    %
    % (c) Achlioptas, Corman, Guibas - 2015  -  http://www.fmaplib.org
    
    properties (GetAccess = public, SetAccess = private)
        M;               % (Mesh) Associated Mesh of LB.        
        W;               % (M.num_vertices x M.num_vertices) Sparse weight matrix of cotangent Laplacian.                
        A;               % (M.num_vertices, 1) The areas associated with the vertices of the M.        
    end
    
    methods (Access = public)          
        function obj = Laplace_Beltrami(in_mesh, vertex_areas)
            % Class constructor.
            % Input:            
            %           in_mesh         -  (Mesh) A triangular Mesh for which the Laplace Beltrami will be constructed.
            %
            %           vertex_areas    -  (optional, in_mesh.num_vertices x 1) A vector carrying the areas associated                         
            %                              with the vertices. If not provided the Mesh is assumed to carry a set of 
            %                              default vertex areas.          
            % Output:
            %           obj             -  (Laplace_Beltrami) object.
            
            obj@Basis();
            if nargin == 0                                               
                obj.M = Mesh();          
                obj.W = [];
                obj.A = [];                
            else
                obj.M       = in_mesh;
                
                if isprop(in_mesh, 'angles')
                    obj.W  = Laplace_Beltrami.cotangent_laplacian(in_mesh.vertices, in_mesh.triangles, in_mesh.angles);
                else
                    obj.W  = Laplace_Beltrami.cotangent_laplacian(in_mesh.vertices, in_mesh.triangles);
                end
                
                if ~ exist('vertex_areas', 'var')
                    vertex_areas = obj.M.get_vertex_areas();
                end
            
                if any(vertex_areas <= 0 )
                    error ('The areas of the vertices provided are not strictly positive.');
                elseif obj.M.num_vertices ~= length(vertex_areas)
                    error ('The number of vertex areas provided is not identical to the number of the mesh vertices.');
                end                
                obj.A = spdiags(vertex_areas, 0, length(vertex_areas), length(vertex_areas));                 
            end
        end
        
        
%         function [evals, evecs] = get_spectra(obj, eigs_num)
%             % Computes the eigenvectors corresponding to the smallest eigenvalues of the Laplace Beltremi operator. 
%             % It stores the results in the spectra property of the LB object (obj.spectra). It also, automatically reuses 
%             % the previously computed ones when this is doable, instead of computing them from the scratch.
%             %
%             % Input:
%             %           eigs_num    -   (int) number of eigenvector-eigenvalue pairs to be computed. This must be
%             %                           smaller than the number of vertices of the associated Mesh.
%             %
%             % Output:
%             %           evals       -   ()
%             %           evecs       -   ()
%             
%             % The requested number of eigenvalues is larger than what has been previously calculated.
%             if length(obj.spectra.evals) < eigs_num      
%                 [evecs, evals] = Laplace_Beltrami.compute_spectra(obj.W, obj.A, eigs_num);                
%                 obj.spectra = struct('evals', evals, 'evecs', evecs);  % Store new spectra.
%            
%             else                                         % Use previously computed spectra.
%                 evals = obj.spectra.evals;
%                 evals = evals(1:eigs_num);                
%                 evecs = obj.spectra.evecs;
%                 evecs = evecs(:, 1:eigs_num);
%             end
%         end
%         
%         function [E] = evecs(obj, eigs_num)
%             % 'Convenience function'.
%             %  Returns only the eigenvectors of the LB. (see get_spectra()).
%             [~, E] = obj.get_spectra(eigs_num);
%         end
%         
%         function [E] = evals(obj, eigs_num)
%             % 'Convenience function'.
%             %  Returns only the eigenvectors of the LB. (see get_spectra()).
%             [E, ~] = obj.get_spectra(eigs_num);
%         end

%         function [R] = synthesize_functions(obj, coeffs)
%             eigs_num = size(coeffs, 1);            
%             R = obj.evecs(eigs_num) * coeffs;            
%         end
                
        function [Proj] = project_functions(obj, eigs_num, varargin)
            %   Projects a set of given functions on the corresponding eigenfunctions of the Laplace Beltrami
            %   operator. Each LB eigenfunction has num_vertices dimensions. I.e., as many as the
            %   vertices of its corresponding Mesh.
            %
            %   Input:
            %           eigs_num    -  (int) number of LB basis functions to be used in the projection.
            %                                      
            %           varargin{i} -  (num_vertices, k{i}) A matrix capturing k{i} functions that will be
            %                          projected on the LB. Each function is a column this matrix.
            %                          
            %   Output: 
            %           Proj        - [sum(k{i}), eigs_num] Matrix carrying the projections of all the functions
            %                         given in matrices in the varargin (thus sum(k{i}) where i=1:num_varargin ) such 
            %                         functions will be outputted. Each one has eigs_num dimensions.
            
            n_varargin = nargin -2; % Number of arguments passed through varargin.
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
            
            [~, evecs] = obj.get_spectra(eigs_num);
            assert(num_vertices  == size(evecs, 1));
                            
            % Project feauture vectors into reduced LB basis.            
            Proj = zeros(eigs_num, functions_total);            % Pre-allocate space.
            right = 0;                     
            projector = evecs' * obj.A;
            
            for i = 1:n_varargin
                left  = right + 1;
                right = right + size(varargin{i}, 2);                                
%                 Proj(:, left:right)  = evecs \ varargin{i};                
                  Proj(:, left:right)  = projector * varargin{i};                            
            end            
        end
        
    end % End of public object-tied functions.
    
    
    methods (Access = public, Hidden = true)
        
        function [Phi, lambda] = compute_spectra(obj, eigs_num)
            % Returns the eigenvalues and the eigenvectors associated with the mesh and its cotangent-LB scheme.

            if eigs_num < 1 || eigs_num > size(obj.W, 1) - 1;
                error('Eigenvalues must be in range of [1, num_of_vertices-1].');
            end            
                        
            sigma = -1e-5; % TODO-V: sigma=0 or 'SM'?
            [Phi, lambda] = eigs(obj.W, obj.A, eigs_num, sigma);
            lambda        = diag(lambda);
            
            if ~isreal(Phi) || ~isreal(lambda)
                error ('Input Mesh leads to an LB operator that has not real spectra.')
            end            
            if sum(lambda < 0) > 1                
                warning ('More than one *negative* eigenvalue were produced. LB is PSD and only the 1st eigenvalue is expected to potentially be smaller than zero (instead of exactly zero).')
            end            
            atol = 1e-08; rtol = +Inf;            
            if ~all_close(Phi' * obj.A * Phi, eye(eigs_num), atol, rtol) 
                error ('The produced eigenvectors are not orthogonal wrt. the Inner Product defined by the vertex areas.')
            end
  
            lambda        = abs(lambda);
            [lambda, idx] = sort(lambda);
            Phi           = Phi(:, idx);                        
        end
    end
    
    
    
    methods (Static)
        
        function [W] = cotangent_laplacian(V, T, varargin)
                % Computes the cotangent laplacian weight matrix. Also known as the stiffness matrix.
                % Input: 
                %                 varargin  - (optional) Angles of the mesh triangles.                
                % Output: 
                %                 W is symmetric. 
                                
                I = [T(:,1); T(:,2); T(:,3)];
                J = [T(:,2); T(:,3); T(:,1)];        
                              
                if nargin == 2
                    angles = Mesh.angles_of_triangles(V, T);                    
                elseif nargin == 3
                    angles = varargin{1};
                else
                    error('Too many arguments were given.')
                end
                
                S = 0.5 * cot([angles(:,3); angles(:,1); angles(:,2)]);
                In = [I; J; I; J];
                Jn = [J; I; I; J];
                Sn = [-S; -S; S; S];
                
                nv = size(V, 1);
                W  = sparse(In, Jn, Sn, nv, nv);
                assert(isequal(W, W'))                
        end
        
    end
    
end
