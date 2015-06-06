function n = normv(A)
% Compute the euclidean norm of each row of A.
%
% Input    A - (m x n) matrix
%
% Output   n - (m x 1) vector, n(i) = ||A(i,:)||_2

    n = sqrt(sum(A.^2, 2));
end