classdef Optimization
    % Some high level routines for doing numerical optimization. Still
    % figuring this out. TODO-V, Do your best:)
  
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
            X = X(ids_y, ids_z);
        end

    end
end
