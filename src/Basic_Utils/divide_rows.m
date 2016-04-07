function B = divide_rows(A, v)
% Divides each row of A by the corresping entry of the vector v.
%
% Input    A - (m x n) matrix
%          v - (1 x n) or (n x 1) vector with non zero, real entries.
%
% Output   B - (m x n) matrix, B(i,:) = A(i,:) / v(i)
%
% (c) Panos Achlioptas 2014    http://www.stanford.edu/~optas

    B = divide_columns(A', v);
    B = B';
end