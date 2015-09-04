classdef Laplacian < Basis
    % All the goodies around the Laplacian of a graph.
    %
    % (c) Achlioptas, Corman, Guibas  - 2015  -  http://www.fmaplib.org
    
    properties (SetAccess = private)       
        L;          % (n x n)  The Laplacian matrix.
        type;       % (String) Describes the type of Laplacian (See constructor for valid values).
        G;          % (Graph)  The corresponding graph from which the Laplacian is derived.        
    end
            
    methods (Access = public)
        % Class Constructor.
        function obj = Laplacian(in_graph, laplacian_type)    
            % Input:
            %           in_graph         - (Graph)  Graph over which the Laplacian is defined.
            %
            %           laplacian_type   - (String) Describes the desired type of laplacian.
            %                              Valid values:   
            %                               'comb' - Combinatorial (unormalized) Laplacian.
            %                               'norm' - Symmetric Normalized Laplacian.
            %                               'sign' - Signless Laplacian.
            % Output:
            %           obj             -  (Laplacian) object.
            obj@Basis();
            if  nargin == 0
                obj.L = [];
                obj.type = '';
                obj.G = Graph();                
            else
                obj.G     = in_graph;
                obj.L     = Laplacian.adjacency_to_laplacian(obj.G.A, laplacian_type);
                obj.type  = laplacian_type;
            end
        end
        
        function [Proj] = project_functions(obj, eigs_num, varargin)
            %   Projects a set of given functions on the corresponding (small) eigenfunctions of the Laplacian.
            %   Each Laplacian eigenfunction has obj.G.num_vertices dimensions. I.e., as many as the vertices of its 
            %   corresponding Graph.
            %   
            %   Input:
            %           eigs_num    -  (int) number of basis functions to be used in the projection.
            %                                      
            %           varargin{i} -  (num_vertices, k{i}) A matrix capturing k{i} functions that will be projected 
            %                          on the LB. Each function is a column this matrix.
            %                          
            %   Output: 
            %           Proj        - [sum(k{i}), eigs_num] Matrix carrying the projections of all the functions
            %                         given in matrices in the varargin (thus sum(k{i}) where i=1:num_varargin ) such 
            %                         functions will be outputted. Each one has eigs_num dimensions.
            
            n_varargin = nargin -2; % Number of arguments passed through varargin.
            if n_varargin < 1
                error ('Please provide some functions to be projected on the Laplacian basis.');
            end
            
            num_vertices = size(obj.L, 1);            
            functions_total = 0;
            for i=1:n_varargin
                if size(varargin{i}, 1) ~= num_vertices                    
                    error ('Wrong dimensionality. The functions must be defined over the nodes of the corresponding graph.')
                end
                functions_total = functions_total + size(varargin{i}, 2);
            end
            
            [~, evecs] = obj.get_spectra(eigs_num);
            assert(num_vertices  == size(evecs, 1));
            assert(num_vertices  == obj.G.num_vertices);
            % Project feature vectors into Laplacian basis.            
            Proj = zeros(eigs_num, functions_total);            % Pre-allocate space.
            right = 0;                                 
            for i = 1:n_varargin
                left  = right + 1;
                right = right + size(varargin{i}, 2);                                
%                 Proj(:, left:right)  = evecs \ varargin{i};                
                  Proj(:, left:right)  = evecs' * varargin{i};                            
            end            
        end

    end % End of public object-tied functions.
    
    methods (Access = public, Hidden = true)
        function [Phi, lambda] = compute_spectra(obj, eigs_num)
            % Returns the smallest eigenvalues and their corresponding eigenvectors of the Laplacian matrix.
            % The user should use the inherited Basis.get_spectra() function instead since that functions wrap around
            % this one, and avoids unnessary computations.
                       
            if eigs_num < 1 || eigs_num > size(obj.L, 1) - 1;
                error('Eigenvalues must be in range of [1, num_of_vertices-1].');
            end            
                        
            sigma = -1e-5; % TODO-V: sigma=0 or 'SM'?
            [Phi, lambda] = eigs(obj.L, eigs_num, sigma);
            lambda        = diag(lambda);
            
            if ~isreal(Phi) || ~isreal(lambda)
                error ('Laplacian has not real spectra.')
            end            
            if sum(lambda < 0) > 1                
                warning ('More than one *negative* eigenvalue were produced. Laplacian is PSD and only the 1st eigenvalue is expected to potentially be smaller than zero (instead of exactly zero).')
            end            
            atol = 1e-08; rtol = +Inf;            
            if ~all_close(Phi' * Phi, eye(eigs_num), atol, rtol) 
                error ('The produced eigenvectors are not orthogonal.')
            end
  
            lambda        = abs(lambda);
            [lambda, idx] = sort(lambda);
            Phi           = Phi(:,idx);                        
        end
    end % End of object-tied functions.
        
    methods (Static)
        function [L] = adjacency_to_laplacian(A, laplacian_type)    
            % Computes the laplacian matrix for a graph described by its adjacency matrix.
            %  
            % Input:    A                - (n x n) Adjacency matrix of a graph with n nodes.
            %           laplacian_type   - (String, optional) Describes the desired type of laplacian.
            %
            %                           Valid values:                           
            %                               'comb' - Combinatorial (unormalized) laplacian (Default value).
            %                               'norm' - Symmetric Normalized Laplacian.
            %                               'sign' - Signless Laplacian.
            %                                       
            % Output:   L               - (n x n) sparse matrix of the corresponding laplacian.
            %                            
            % Notes:  
            %       DOI: "A Tutorial on Spectral Clustering, U von Luxburg".
            %
            % (c) Panos Achlioptas 2015  -  http://www.stanford.edu/~optas/FmapLib

            if(size(A, 1) ~= size(A, 2))
                error('Adjacency matrices are square, but a non square matrix was given.')
            end

            if nargin < 2
                laplacian_type = 'comb';  % Default output is the combinatorial Laplacian.
            end

            n = size(A, 1);
            total_weight = sum(A, 2);            
            D = spdiags(total_weight, 0, n, n);

            if strcmp(laplacian_type, 'comb')         
                L = -A + D;      
            elseif strcmp(laplacian_type, 'norm')
                Dn = spdiags (1 ./ sqrt(total_weight), 0, n, n);
                L  = Dn *( -A + D) * Dn;
            elseif strcmp(laplacian_type, 'sign')         
                L = A + D;      
            else
                error('Please provide a valid argument for the type of laplacian you want.')
            end
        end        
    end
       
end

