function [N] = l2_norm(A)
%   Computes the euclidean norm of each row of A.
%
%   Input:
%           A - (m x n) matrix.
%
%   Output:
%           N - (m x 1) vector, N(i) = ||A(i,:)||_2

    if size(A, 2) == 1   % A is a column vector.
        A = A';
    end
    N = sqrt(sum(A.^2, 2));
end