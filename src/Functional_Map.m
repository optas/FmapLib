classdef Functional_Map
    % A class implementing a variety of utilities related to the Functional
    % Maps framework.
        
    properties
    end
    
    methods (Static)
        
        
        function [centers, from_radius, to_radius] = ball_distortion_of_map(inmap, from_mesh, to_mesh, ngeoballs)
            % Computes the distortion of a geodesic ball by the given map. 
            % The centers of the balls are uniformally distributed random
            % points. For growing radius of the initial geodesic ball correspond the
            % radius of the mapped geodesic balls.
            % 
            % Input:
            %           inmap       -   (num_vertices x 1) The i-th vertex
            %                           of from_mesh correspond to the
            %                           inmap(i)-th vertex of to_mesh.
            %           from_mesh   -   (Mesh) Input mesh. 
            %           to_mesh     -   (Mesh) Input mesh. 
            %           ngeoballs   -   (int)  Number of geodesic balls to
            %                            be produced.                              
            %
            % Output:   centers     -   (ngeoballs x 1) List of indeices on from_mesh
            %                            corresponding to the center of a
            %                            geodesic balls.
            %           from_radius -   (num_vertices x ngeoballs) Growing
            %                           radius of a geodesic ball on
            %                           from_mesh.
            %           to_radius   -   (num_vertices x ngeoballs)
            %                           Corresponging radius on to_mesh.
            
            % TODO-E Test
            centers = Functional_Map.random_delta_functions(from_mesh, ngeoballs, 1);
            
            from_radius = zeros(from_mesh.num_vertices, ngeoballs);
            to_radius   = zeros(from_mesh.num_vertices, ngeoballs);
            
            for i = 1:ngeoballs
                from_dists                  = comp_geodesics_to_all(from_mesh.vertices(:,1), from_mesh.vertices(:,2), from_mesh.vertices(:,3), from_mesh.triangles', centers(i));
                [from_radius(:,i), id_from] = sort(from_dists);
                
                to_dists       = comp_geodesics_to_all(to_mesh.vertices(:,1), to_mesh.vertices(:,2), to_mesh.vertices(:,3), to_mesh.triangles', inmap(centers(i)));
                to_dists       = to_dists( inmap(id_from) );
                to_radius(:,i) = cummax(to_dists); % WARNING exists only in MATLAB15a (redefined below)
            end
        end
        

       function [dists] = pairwise_distortion_of_map(inmap, from_mesh, to_mesh, from_basis, to_basis, nsamples, fast)                                
           %% Document.
           %  Symmetries? e.g., decode Rodola's file           
            
%             proj_deltas = (from_basis \ deltas);                           % Delta functions projected in the from_basis.                                                
%             deltas_transfered = inmap * proj_deltas;                        
%             all_to_deltas = (to_basis \ speye(size(to_basis, 1)));                     
%             [ids, dist]  = knn(all_to_deltas'  ,deltas_transfered');         %TODO-P, 'IncludeTies', 'True');                       
            
            deltas            = Functional_Map.random_delta_functions(from_mesh, nsamples, 1);
            proj_deltas       = from_basis(deltas, :)';                     % Delta functions expressed in from_basis.
            deltas_transfered = inmap * proj_deltas;                        % Use inmap to transfer them in to_mesh.
            [ids, dist]       = knnsearch(to_basis', deltas_transfered');   % Find closest function for its tranfered on (Euclidean dist is used).
                                                                            % TODO-P,E solve 'Ties' in knn.                                              

            pairs = [ids, groundtruth(deltas)];                                                           
            if exist('fast' ,'var') && fast ~= 0  % compute true geodesics or use approx. by Dijkstra.
                dists = comp_geodesics_pairs(to_mesh.vertices(:,1), to_mesh.vertices(:,2), to_mesh.vertices(:,3), to_mesh.triangles', pairs, 1);
            else
                dists = comp_geodesics_pairs(to_mesh.vertices(:,1), to_mesh.vertices(:,2), to_mesh.vertices(:,3), to_mesh.triangles', pairs);
            end                            
        end

        function [S] = random_delta_functions(inmesh, nsamples, index_only)
            % Computes uniformly i.i.d. random delta functions of the given mesh
            % vertices. The delta function of vertex -i- is a vector 
            % which has a single non-zero entry at its i-th dimension. 
            % The total number of dimensions of such a vector is equal 
            % to the number of vertices of the given mesh. The sampling is
            % done without replacement.
            % 
            % Input:
            %           inmesh      -   (Mesh) Input mesh.             
            %           nsamples    -   (int)  Number of delta functions to
            %                            be produced.            
            %           index_only  -   (optional, Boolean) If 1, only the indices 
            %                            of the vertices are returned.
            %                                                                                     
            %
            % Output:   S           -     (num_vertices x nsamples) Sparse matrix.            
            %                       - or, (nsamples, 1) if index_only == 1.
            %                          
            % Precondition: nsamples must be at most as large as
            %               inmesh.num_vertices.
                       
            vertices = randsample(inmesh.num_vertices, nsamples);           
            if exist(index_only, 'var') && index_only == 1
                S = vertices;
                return
            end
            S  = sparse(vertices, 1:nsamples, ones(1,nsamples), inmesh.num_vertices, nsamples);                       
        end
        
        
        function X = sum_of_squared_frobenius_norms(D1, D2, L1, L2, lambda)
            % This code uses plain least squares techniques to 
            % solve the objective function using Frobenius norm squared
            % terms.
            % We aim to minimize ||X*D1-D2||^2 + lambda*||X*L1-L2*X||^2. 
            % This is done by computing the gradient wrt X, which is given
            % by (X*D1-D2)*D1' for term 1, and by 
            % (((lambda1_i-lambda2_j)^2)*X_{ij})_{ij} for term 2. This is
            % solved using linear techniques as given below.
            
            N1 = size(D1, 1);
            N2 = size(D2, 1);
            
            A_fixed = D1 * D1' ;
            B = D1 * D2' ;
            X = zeros(N2, N1);
            
            if lambda == 0   %  Un-regularized                
                X = (A_fixed\B)';
            else
                for i = 1 : N2
                    A = diag(lambda * (L1 - L2(i)) .^ 2) + A_fixed;
                    X(i, :) = A \ B(:, i);
                end
            end
        end
          
        function X = sum_of_frobenius_norms(D1, D2, L1, L2, lambda)
            % Copyright (C) 2014 Fan Wang.            
            % TODO-P: Add documentation.
            % This code uses Sedumi to write the objective function as a
            % Semi-definite program in the dual form. In this problem we
            % aim to minimize ||X*D1-D2|| + lambda*||X*L1-L2*X||. This is
            % done by solving for min t1 + lambda*t2, s.t t1 >=
            % ||X*D1-D2||, t2 >= ||X*L1-L2*X||.
            
            N1 = size(D1, 1);
            N2 = size(D2, 1);
            N1N2 = N1 * N2;
            % cvx_setspath('sedumi');
            z = sparse(N1N2, 1);
                       
            dim = numel(D2);
            K.q = [1 + dim];
            b   = [z; -1];
            D2t = sparse(D2');
            c   = [0; D2t(:)];
            At = [z' -1;
                  kron(speye(N2), D1') sparse(dim, length(K.q))];

            % L1 L2
            if lambda > 0
                dim = N1N2;
                K.q = [K.q 1 + dim];
                b = [b; -lambda];
                c = [c; 0; z];
                At = [At sparse(size(At, 1), 1);
                    sparse(1, size(At, 2)) -1;
                    diag(kron(ones(N2, 1), L1) - kron(L2, ones(N1, 1))) sparse(dim, length(K.q))];
            end

            % Solve
            % The objective for dual is max b'*y which in this case is
            % -(y_{N1N2+1} + lambda*y_{N1N2+2}). t1 = y_{N1N2+1}, t2 =
            % y_{N1N2+2}. 
            
            %The dual y satisfies s = c-A'*y, s \in Dual cone K'. The
            % existence of s in the dual gives the conditions relating
            % t1,t2 to the Frobenius norms in the objective.
            [~, y, ~] = sedumi(sparse(At), sparse(b), sparse(c), K, struct('fid', 0));
            X = reshape(y(1:N1N2), N1, N2)';
        end

        
        function X = sum_of_squared_frobenius_norms_non_diagonal(D1, D2, L1, L2, lambda)
            N1 = size(D1,1);
            N2 = size(D2,1);
            A_fixed = D1 * D1' ;
            B = D2 * D1' ;
            % The first N2 rows of A_large correspond to the first row of
            % AX', or in other words (a_1)*X', the corresponding B values
            % being the first column of B.
            A_large = kron(A_fixed,eye(N2));
            B_large = B(:);
            for n=1:N1
                for m=1:N2
                    L = zeros(N2,N1);
                    for i=1:N2
                        for j=1:N1
                            switch j
                                case m
                                    switch i
                                        case n
                                            L(i,j) = 2*(sum_square(L1(n,:))-sum_square(L2(:,m))-2*L1(n,n)*L2(m,m)-(L1(n,n)-L2(m,m))^2);
                                        otherwise,
                                            L(i,j) = 2*(L1(j,:)*L1(n,:)'-L2(m,m)*L1(n,j)-L1(n,n)*L1(j,n));
                                    end
                                otherwise,
                                    switch j
                                        case n
                                            L(i,j) = 2*(L2(:,i)'*L2(:,m)-L1(n,n)*L2(i,m)-L2(m,i)*L2(m,m));
                                        otherwise,
                                            L(i,j) = 2*(-L2(m,i)*L1(n,j)-L2(i,m)*L1(j,n));
                                    end
                            end
                        end
                    end
                    
                    A_large((n-1)*N2+m,:) = A_large((n-1)*N2+m,:) + lambda*L(:)';
                end
            end
            
            X_vec = A_large \ B_large;
            X = reshape(X_vec,N2,N1);
        end
    
    
        function [maps_new, X_all] = low_rank_filtering(maps, W)
            % Panos figuring this out:

            % maps      -   mn x mn matrix with all initial pairwise maps (m x m)
            %               between n objects.

            % W         -   n x n: global similarity between each pair of objects. When
            %               two objects are very similar, being inconcsisten will be
            %               penalized more.

            % X_all     -   must be the collection of the sequential mn x mn matrices that correspond
            % to the updates of the low rank optimization towards convergence.
            % maps_new  -   must be the last update, i.e., the X_all(:,:,end)

            n = size(W, 1);
            m = size(maps{1,2}, 1);

            % Initialization of X0 to a block matrix carrying the initial maps.
            X_0 = kron(ones(n,n), eye(m));
            for rowIndex = 1:n
                rowIndices = ((rowIndex-1)*m+1):(rowIndex*m);
                for colIndex = 1:n
                    if rowIndex == colIndex
                        continue;
                    end
                    colIndices = ((colIndex-1)*m+1):(colIndex*m);
                    X_0(rowIndices, colIndices) = maps{rowIndex, colIndex};                
                end
            end

            lambda = kron(W+eye(n), ones(m,m))/sqrt(m*n);

            % low rank optimization
            X_all = Optimization.rank_min(X_0, lambda, n, m);

            X = X_all(:,:,end);
            % write it back
            maps_new = maps;
            for rowIndex = 1:n
                rowIndices = ((rowIndex-1)*m+1):(rowIndex*m);
                for colIndex = 1:n
                    if rowIndex == colIndex
                        continue;
                    end
                    colIndices = ((colIndex-1)*m+1):(colIndex*m);
                    maps_new{rowIndex, colIndex} = X(rowIndices, colIndices);
                end
            end
        end

    end % Static.
end % ClassDef.