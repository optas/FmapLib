classdef IS
    % A class offering a variety of static functions of the form  B = IS.x(A), where x is some binary property testable 
    % on A.      
    
    methods (Static)        
        
        function B = dyadic(n)
            % Checks if a number is a positive power of 2, i.e. 
            % if n = 2^k for some positive integer k.
            %
            % Input:        n - (m x n) matrix, usually just a real. 
            %
            % Output:       B - (m x n) equal to 1 in all entries that are             
            %                   dyadic and 0 elsewhere.
            %
            % Precondition: Every entry of n is less than 2^30.
            %
            % Note: for n > 2^30 the function has not being tested and
            % numerical stability should be taken into account (try
            % n = 2^50 + 1)            
            assert(n < 2^30)
            if n <= 0
                B = 0;
            else
                B = floor(log2(n)) == log2(n);
            end
        end
        
        function B = integer(n)
            % Checks if a number is a positive or negative integer.            
            %
            % Input:        n - (m x n) matrix, usually just a real. 
            %
            % Output:       B - (m x n) equal to 1 in all entries that are
            %                   integers and 0 elsewhere.
            %                              
            % Precondition: n is real.
            B = ~mod(n, 1);
        end
        
        function B = permutationaly_supported(A)
            % True IFF input A is a square matrix that has exactly
            % a single non-zero element in each of its rows and a single non
            % zero element in each of its columns. In this case the non
            % zero elements are supported in a permutation matrix of same
            % size. 
            % 
            % More info: http://en.wikipedia.org/wiki/Generalized_permutation_matrix
            %             
            % Input:    A - (n x n) real or comlex matrix.
            % 
            % Output:   1 if A is permutationaly supported, else 0.
            %             
            % Precondition: A is a square matrix.             
            n = size(A, 1);
            assert(n == size(A,2));
            B = A ~= 0;            % Non-zero elements binary mask
            col_sum = sum(B);
            row_sum = sum(B, 2);
            
            if all(col_sum == ones(1,n)) && all(row_sum == ones(n,1))
                B = true;
            else
                B = false;
            end
        end
        
        function B = full_rank(A)
            % Checks if the input matrix is full rank.
            %
            % Input:    A - (m x n) real or comlex matrix.
            % 
            % Output:   1 if the rank of A is m or n, else 0.            
            [~, D, ~] = svd(full(A));
            D = abs(diag(D));
            B = min(D(:)) > .01 .* mean(D(:));
        end
        
        function B = aprox_within_e(X1, X2, norm_type, aprox)
            n = size(X1,1);
            assert(n == size(X1,2))    % X1,X2 are square matrices.
            diff = norm(X1-X2, norm_type);
            if diff < n*aprox
                B = true;
            else
                B = false;
            end
        end
        
        function B = binary(A)
            ones  = sum(A(:) == 1);
            zeros = sum(A(:) == 0);            
            B = ones + zeros == numel(A);
        end
            
        function B = percent(A)
            B = all(A(:) <=1) && all(A(:) >= 0);
            
        end
        
        function B = positive_definite(A)
             [~, B] = chol(A);              % Computationally more efficient that eig(A).
             B = ~B;
        end
        
        function B = psd(A)
            B = min(eig(A)) >= 0;
        end
        
        function B = non_increasing(A)            
            B = all(diff(A) <= 0);
        end
        
        function B = non_decreasing(A)
            B = all(diff(A) >= 0);
        end
        
        function B = decreasing(A)            
            B = all(diff(A) < 0);
        end
        
        function B = increasing(A)
            B = all(diff(A) > 0);
        end
        
        function B = single_number(A)            
            B = (numel(A) == 1) && isnumeric(A);
        end
        
        function B = vector(A)            
            B = ismatrix(A) && min(size(A)) == 1;
        end
        
    end
    
end

