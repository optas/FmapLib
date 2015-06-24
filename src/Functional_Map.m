classdef Functional_Map < dynamicprops
    % A class implementing a variety of utilities related to the Functional
    % Maps framework.

    properties (GetAccess = public, SetAccess = private)
        % Basic properties that every instance of the Functional_Map class has.
        source_basis = [];        
        target_basis = [];                        
        fmap         = [];
        source_neigs = 0;
        target_neigs = 0;
    end
    
    methods (Access = public)
        % Class Constructor.        
        function obj = Functional_Map(varargin)     
            if nargin == 0                
                % Construct an empty Mesh.            
                obj.source_basis = [];
                obj.target_basis = [];                               
            else 
                obj.source_basis = varargin{1};
                obj.target_basis = varargin{2};
            end
            obj.fmap = [];
            obj.source_neigs = 0;
            obj.target_neigs = 0;
        end
        
        function [F] = compute_f_map(obj, method ,neigs_source, neigs_target, source_feat, target_feat, varargin)
            ns = obj.source_basis.M.num_vertices;  % Number of vertices on source.
            nt = obj.target_basis.M.num_vertices;  % Number of vertices on target.
            if size(source_feat, 1) ~= ns || size(target_feat, 1) ~= nt
                error('The feature functions are supposed to be given as vectors and their dimensions should equal the number of corresponding mesh vertices.')
            end
            
            options = struct('normalize', 1, 'area_normalize', 0);
            option_names = fieldnames(options); % Read the acceptable names.
            nvargs = length(varargin);
            if round(nvargs/2) ~= nvargs/2
                error('Expecting property_name/property_value pairs.')
            end
            for pair = reshape(varargin, 2, []) % Pair is {propName;propValue}.
                inp_name = lower(pair{1}); % Make case insensitive.
                if any(strcmp(inp_name, option_names))
                    options.(inp_name) = pair{2};
                else
                    error('%s is not a recognized parameter name.',inp_name)
                end
            end
            
            switch method
                case 'functions_only'                    
                    source_feat = obj.source_basis.project_functions(neigs_source, source_feat);
                    target_feat = obj.target_basis.project_functions(neigs_target, target_feat);                                        
                    
                    if options.normalize == 1 
                        source_feat = normc(source_feat);
                        target_feat = normc(target_feat);
                    else options.area_normalize == 1 
                        error('TODO-E normalize them to have unit norm wrt to source-areas.')                        
                    end
                    
                    [F] = Functional_Map.sum_of_squared_frobenius_norms(source_feat, target_feat, [], [], 0);
            end
            obj.fmap         = F;
            obj.source_neigs = neigs_source;
            obj.target_neigs = neigs_target;
        end
        
        function [D] = area_difference(obj)
            if isempty(obj.fmap)
                error('It appears that this object currently is not carrying a matrix corresponding to a functional map.')
            end
            
            transfered_basis = obj.target_basis.evecs(obj.target_neigs) * obj.fmap;
            target_inner_prod = transfered_basis' * obj.target_basis.A * transfered_basis;            
            source_evecs = obj.source_basis.evecs(obj.source_neigs);            
            source_inner_prod = source_evecs' * obj.source_basis.A * source_evecs;                                   
            D = pinv(source_inner_prod) * target_inner_prod;        
        end
        
        function [D] = area_difference2(obj)
            if isempty(obj.fmap)
                error('It appears that this object currently is not carrying a matrix corresponding to a functional map.')
            end            
            D = obj.fmap'* obj.fmap;       
        end
        
%         
%         function [D] = conformal_difference(obj)
%             if isempty(obj.fmap)
%                 error('It appears that this object currently is not carrying a matrix corresponding to a functional map.')
%             end
%             obj.fmap
%             D = obj.source_basis.evecs(:, )
%         end
%         
        
        function [dists, indices] = pairwise_distortion(obj, groundtruth, varargin)
            [dists, indices] = Functional_Map.pairwise_distortion_of_map(obj.fmap, obj.source_basis, obj.target_basis, groundtruth, varargin{:});
        end
        
        
%         function [F] = groundtruth_fmap(obj)
%             F = Functional_Map(obj.source_basis, obj.target_basis)
%             
%             F.fmap = Functional_Map.groundtruth_functional_map(source_basis.evecs, basis_to, gt_from_to, to_areas)                        
%             
%             
%             function [X] = 
%             
%         end
        
    end
    
    methods (Static)                
        
        function [dists, indices] = pairwise_distortion_of_map(inmap, source_basis, target_basis, groundtruth, varargin)
            %% Document.            
            % inmap - matrix corresponding to a functional map
            % source_basis - LB basis of source
            % target_basis - LB basis of source
            % groundtruth  - 
            
            options = struct('fast', 1, 'nsamples', 100, 'indices', [], 'symmetries' , []);        
            option_names = fieldnames(options); % Read the acceptable names.
            nvargs = length(varargin);
            if round(nvargs/2) ~= nvargs/2
                error('Expecting property_name/property_value pairs.')
            end
            for pair = reshape(varargin, 2, []) % Pair is {propName;propValue}.
                inp_name = lower(pair{1});      % Make case insensitive.
                if any(strcmp(inp_name, option_names))
                    options.(inp_name) = pair{2};
                else
                    error('%s is not a recognized parameter name.',inp_name)
                end
            end

            if ~isempty (options.indices)  % TODO-P add type_checking.
                [deltas, ~]       = Functional_Map.random_delta_functions(diag(source_basis.A), 0, options.indices);
                indices           = options.indices;
            else                
                [deltas, indices] = Functional_Map.random_delta_functions(diag(source_basis.A), options.nsamples);                           
            end
            
            s_neigs = size(inmap, 2);   % Number of source basis used to constuct the inmap.
            t_neigs = size(inmap, 1);
            
            proj_deltas = source_basis.project_functions(s_neigs, deltas); % Project random delta function on source basis.           
            deltas_transfered = inmap * proj_deltas;                       % Use inmap to transfer them in target_mesh.            
            
            target_deltas     = target_basis.evecs(t_neigs)' * sqrt(target_basis.A); % This is the set of all delta
                                                                                     % functions defined on target
                                                                                     % and expressed in the
                                                                                     % target_basis.
                                                                                     
            [ids, ~]          = knnsearch(target_deltas' , deltas_transfered');      % Find closest function for its tranfered on (Euclidean dist is used).
                                                                                                      % TODO-P,E solve 'Ties' in knn.                                              
                                                                                                                       
            pairs = [ids, groundtruth(indices)]';                                                           
            dists = comp_geos(pairs);

            if ~ isempty(options.symmetries)
                for i=1:size(options.symmetries, 2)
                    pairs = [ids, options.symmetries(indices, i)]';                                                           
                    dists = min(comp_geos(pairs), dists);
                end                
            end
            
            function [dists] = comp_geos(pairs)
                vertices  = target_basis.M.vertices;
                triangles = target_basis.M.triangles;
                if options.fast == 1                                            % Compute true geodesics or use approx. by Dijkstra.                            
                    dists = comp_geodesics_pairs(vertices(:,1), vertices(:,2), vertices(:,3), triangles', pairs, 1);
                else 
                    dists = comp_geodesics_pairs(vertices(:,1), vertices(:,2), vertices(:,3), triangles', pairs);
                end                
            end
            
        end
        
                
        function [centers, from_radius, to_radius] = ball_distortion_of_map(inmap, source_mesh, target_mesh, ngeoballs)
            % Computes the distortion of a geodesic ball by the given map. 
            % The centers of the balls are uniformally distributed random
            % points. 
            % TODO-E,P: Improve: For growing radius of the initial geodesic ball correspond the
            % radius of the mapped geodesic balls.
            % 
            % Input:
            %           inmap       -   (num_vertices x 1) The i-th vertex
            %                           of source_mesh correspond to the
            %                           inmap(i)-th vertex of target_mesh.
            %           source_mesh   -   (Mesh) Input mesh. 
            %           target_mesh     -   (Mesh) Input mesh. 
            %           ngeoballs   -   (int)  Number of geodesic balls to
            %                            be produced.                              
            %
            % Output:   centers     -   (ngeoballs x 1) List of indeices on source_mesh
            %                            corresponding to the center of a
            %                            geodesic balls.
            %           from_radius -   (num_vertices x ngeoballs) Growing
            %                           radius of a geodesic ball on
            %                           source_mesh.
            %           to_radius   -   (num_vertices x ngeoballs)
            %                           Corresponging radius on target_mesh.

            % TODO-E Test
            centers = Functional_Map.random_delta_functions(source_mesh, ngeoballs, 1);

            from_radius = zeros(source_mesh.num_vertices, ngeoballs);
            to_radius   = zeros(source_mesh.num_vertices, ngeoballs);

            for i = 1:ngeoballs
                from_dists                  = comp_geodesics_to_all(source_mesh.vertices(:,1), source_mesh.vertices(:,2), source_mesh.vertices(:,3), source_mesh.triangles', centers(i));
                [from_radius(:,i), id_from] = sort(from_dists);

                to_dists       = comp_geodesics_to_all(target_mesh.vertices(:,1), target_mesh.vertices(:,2), target_mesh.vertices(:,3), target_mesh.triangles', inmap(centers(i)), 1);
                to_dists       = to_dists( inmap(id_from) );
                    to_radius(:,i) = cummax(to_dists);
            end
        end
        
        
        function [S, indices] = random_delta_functions(vertex_areas, nsamples, indices)
            % Computes uniformly i.i.d. random delta functions of the given mesh vertices. The delta function of 
            % vertex -i- is a vector which has a single non-zero entry at its i-th dimension. The actual value at
            % the i-th dimension is equal to 1/sqrt(area_of_vertex(i)). This value is chosen to sasisfy the equation:
            %             delta_function(i)' * diag(Vertex_Areas) * delta_function(i) = 1. 
            % I.e., it makes a delta function to have a unit of mass wrt. the area inner product of a mesh.
            % The total number of dimensions of a delta function/vector is equal to the number of vertices of the given 
            % mesh. Finally, the sampling of the vertices is done done without replacement.
            % 
            % Input:
            %           vertex_areas  -   (num_vertices x 1) The areas associated with each vertex of a Mesh.
            %              
            %           nsamples      -   (int)  Number of delta functions to be produced.
            %           
            %           indices       -   (optional, k x 1) k positive integers describing the vertices over which
            %                             the delta functions will be produced. When used this function is deterministic.
            %                            
            % Output:   
            %           S             -   (num_vertices x nsamples) Sparse matrix with the delta functions as column 
            %                             vectors.
            %                                       
            %           indices       -   (nsamples, 1) The vertices of the deltas (equals to 'indices' if passed).
            %
            % Precondition: nsamples must be at most as large as length(vertex_areas).
 
            num_vertices = length(vertex_areas);
            
            if ~exist('indices' ,'var')
                indices = randsample(num_vertices, nsamples);
            else
                nsamples = length(indices);
            end                        
            areas   = 1 ./ sqrt(vertex_areas(indices));
            S       = sparse(indices, 1:nsamples, areas, num_vertices, nsamples);
        end
        
        
        
        function [X] = groundtruth_functional_map(basis_from, basis_to, gt_from_to, to_areas)                        
%             nodes_from = size(basis_from, 1);
%             nodes_to   = size(basis_to, 1);              
%             non_zero   = length(correspondences_from_to(:, 2));
%             P          = sparse(correspondences_from_to(:, 2), correspondences_from_to(:, 1), ones(non_zero,1), nodes_to, nodes_from);            
%             X          = basis_to' * P * basis_from;
            
            basis_from = basis_from(gt_from_to ~= 0, :) ;              % Remove dimensions for which you do not know the groundtruth (i.e., map -ith- vertex to 0).
            basis_from = basis_from(gt_from_to(gt_from_to ~= 0), :);   % Permute the basis to reflect the corresponding ground_truth.            
            A          = spdiags(to_areas, 0, length(to_areas), length(to_areas));      %TODO-P: we added the areas to have the real pinv of LB1.            
            X          = basis_to' * A * basis_from;                       
%             X          = basis_to' * A * basis_from * A_to;  %TODO-E,P                       
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