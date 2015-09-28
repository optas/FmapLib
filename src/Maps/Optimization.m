classdef Optimization
    % Some high level routines for doing numerical optimization. Work in progress.
  
    methods (Static)
  
        function [X] = linear_sdp_admm(X_0, C, W, ids_A, b, X_init)
            % This function solves the following optimization problem
            % min <C, X> + <W,E> 
            % subject to
            %     X is semidefinite
            %     X - X_0 <= E
            %     X_0 - X <= E
            %     X(ids_A) == b

            dim = size(X_0, 1);

            y = zeros(length(ids_A), 1);
            Z1 = zeros(dim, dim);
            Z2 = zeros(dim, dim);
            E = zeros(dim, dim);
            S = zeros(dim, dim);
            X = zeros(dim, dim);
            if size(X_init, 1) > 0
                X = X_init;
            end

            mu = sqrt(dim)/20;
            rho = 0.99;
            t = cputime;
            for iter = 1:160
                TP = Z1 + C - S - Z2 - mu*X;
                y = b*mu + TP(ids_A);
                A_star_y = zeros(dim, dim);
                A_star_y(ids_A) = y;
                A_star_y = A_star_y + A_star_y';
                for i=1:dim
                    A_star_y(i,i) = A_star_y(i,i)/2;
                end
                Z1 = (W + (A_star_y + S - C) - mu*(X_0 + E - X))/2;
                Z2 = (W - (A_star_y + S - C) - mu*(-X_0 + E + X))/2;
                Z1 = max(Z1, 0);
                Z2 = max(Z2, 0);
                E = E + (Z1 + Z2 - W)/mu;
                V = Z1 - Z2 + C - A_star_y - X*mu;
               [Q,Sigma] = eig(V);
            %     V1 = gpuArray(V);
            %     [Q1,Sigma1] = eig(V1);
            %     Q = gather(Q1);
            %     Sigma = gather(Sigma1);
                d = diag(Sigma);
                ids = find(d > 0);
                S = Q(:, ids)*Sigma(ids, ids)*Q(:,ids)';
                S = (S+S')/2;
                X_prev = X;
                X = (S - V)/mu;
                Dif = X_prev - X;
                disp = sqrt(sum(sum(Dif.*Dif)))/sqrt(sum(sum(X_0.*X_0)));
                if mod(iter, 20) == 0
                    fprintf(' inner iteration %d\n', iter);
                end
                if disp < 1e-6
                    break;
                end
                mu = mu*rho;
            end
        end

        function [X_all] = rank_min(X_0, W, n, m)            
            dim = n*m;
            C = eye(2*dim)/2;
            W = [zeros(dim, dim), W;
                W', zeros(dim, dim)];

            X_0 = [zeros(dim, dim), X_0;
                   X_0', zeros(dim, dim)];

            b = reshape(kron(ones(1,n), eye(m)), [n*m^2, 1]);
            ids = kron(ones(1,m), 1:m) + dim*2*kron(0:(m-1), ones(1,m));
            ids = (kron(ones(1,n), ids) + (2*dim+1)*m*kron(0:(n-1), ones(1,m^2)))' + 2*dim^2;

            [X] = Optimization.linear_sdp_admm(X_0, C, W, ids, b, []);

            ids_y = 1:dim;
            ids_z = (dim+1):(2*dim);
            X_all = zeros(dim,dim, 11);
            X_all(:,:,1) = X(ids_y,ids_z);

            delta = 1;
            for i = 1:10
                ids_y = 1:dim;
                ids_z = (dim+1):(2*dim);
                Y = X(ids_y, ids_y);
                Z = X(ids_z, ids_z);
                Y = inv(Y + delta*eye(dim));
                Z = inv(Z + delta*eye(dim));
                C(ids_y, ids_y) = Y;
                C(ids_z, ids_z) = Z;
                [X] = Optimization.linear_sdp_admm(X_0, C, W/8, ids, b, X);
                X_all(:,:,i+1) = X(ids_y, ids_z);
            end
            X = X(ids_y, ids_z);   % seems unecessary 
        end
        
        
        function [Y, V] = latent_basis_given_functional_maps(in_maps, weigths, latent_size)
            % Returns a set of orthonormal vectors that correspond to latent commonalities among objects of a
            % collection. These commonalities are found by exploiting a set of functional maps between pairs of related 
            % objects.
            % 
            % Input:
            %           in_maps  -  ( N x N cell array ) carrying the functional map from object i to object j, it its  
            %                       (i,j) position. N is the total number of objects in the collection.
            %
            %           weights  -  ( N x N matrix) weights(i,j) is a double reflecting how closely related object i to 
            %                       object j.
            %
            %           latent_size - (int) Specifies how many latent functions will be derived in each object.
            %
            %
            % Output: _add_
            %
            %
            % Reference: 'Image Co-Segmentation via Consistent Functional Maps, F. Wan et al. in _add_''
            
            % _add_content
            
        end
        
        
        function [X, opt] = least_edges_given_connectivity_first_try(init_graph, dist_matrix, alpha)
            % Computes the adjacency matrix of a graph for which the algrebraic connectivity is larger than alpha and
            % the sum of all its edge weigts is minimized. _add content_
            %            
            % Input:
            %           init_graph   -  
            %                       
            %           dist_matrix  -  (N x N matrix) weights(i,j) is the weight an edge between (i and j).            
            %
            %           alpha        - (int) Lower bound on algebraic connectivity of derived graph.
            %
            %
            % Output: _add_
            %
            %
            % Reference: 'Graph Weight Design for Laplacian Eigenvalue Constraints with Multi-Agent Systems Applications.'
            
            n = init_graph.num_vertices;
            L = Laplacian(init_graph, 'comb');            
            
            [~, U    ] = L.get_spectra(n-1);
            [Ulast, ~] = eigs(L.L, 1);         % Compute the largest eigen-pair.
            
            size(U)
            size(Ulast)
            U = [U(:, 2:end) Ulast];           % Remove the constant eigenvector and instert Ulast.
            size(U)
            
            alpha = alpha * eye(n);
            cvx_begin                        
                variable X(n, n) symmetric                
                minimize trace(X * dist_matrix)
                subject to
                    U.' * ( diag(diag(X * dist_matrix)) - X - alpha) * U >= 0
                    vec(X)  >= 0
                    -vec(X) >= -1                
            cvx_end      
            opt = cvx_optval;                
        end
        
        function [Q] = get_orthogonal_vectors(in_lap, alpha)
                nodes = size(in_lap, 1);
%                 A = in_lap;
                A     = in_lap - diag(repmat(alpha, nodes, 1));                
                [Q, ~] = eigs(A, nodes - 1);                

%                 sigma = 10e-6;
%                 [Q, ~] = eigs(A, 3, sigma);                
% %                 Q = fliplr(Q);
%                 Q = Q(:,1:2)
                
                                
        end
        
        function [X, opt] = least_edges_given_connectivity(dist_matrix, alpha)
            % Computes the adjacency matrix of a graph for which the algrebraic connectivity is larger than alpha and
            % the sum of all its edge weigts is minimized. _add content_
            %            
            % Input:                       
            %           dist_matrix  -  (N x N matrix) weights(i,j) is the weight an edge between (i and j).            
            %
            %           alpha        - (int) Lower bound on algebraic connectivity of derived graph.
            %
            %
            % Output: _add_
            %
            %
            % Reference: 'Graph Weight Design for Laplacian Eigenvalue Constraints with Multi-Agent Systems Applications.'
            tic
            weighted_clique  = Graph(dist_matrix, false);
            n                = weighted_clique.num_vertices;
            e                = weighted_clique.num_edges;
            [node_f, node_t, weights]  = weighted_clique.all_edges();
            W                = spdiags(weights, 0, e, e);
                       
            clique           = Graph.generate('clique', n, n);
                        
            if clique.num_edges ~= e
                error('"dist_matrix" must containt zeros only on its diagonal.')
            end
            
            I = clique.incidence_matrix();            
            toc
                                                
            alpha_i = eye(n)*alpha;
            num_iter = 2;
            
            D = eye(size(I, 2));
            
            for i =1:num_iter                
                Q = Optimization.get_orthogonal_vectors(I*D*I', alpha);
                K = []
                cvx_begin sdp quiet
                    variable K(e, e) diagonal
                    minimize trace(K * W)                
                    subject to
                        Q.' * ( I * K * I' - alpha_i ) * Q == semidefinite(size(Q,2))
                        K(:)  >= 0
                        K(:)  <= 1                    
                cvx_end        
                opt = cvx_optval;                            
%                 K
                D = K;
           end
            
            if opt == +Inf
                warning('Program was not feasible.')
                X = NaN;
                return 
            else                                                
                X = sparse(node_f, node_t, diag(K), n, n);  % Put edges on a adjacency matrix.
                X = X + X';              
                L = Laplacian(Graph(X, false), 'comb');
                lambda = L.evals(2);
                assert(lambda(2) >= alpha)
            end
            
            toc
            
        end
        
        
        minimum_effective_resistance
            

    end
end
