function B = divide_columns(A, v)
% Divides each column of A by the corresping entry of the vector v.
%
% Input    A - (m x n) matrix
%          v - (1 x n) or (n x 1) vector with non zero, real entries.
%
% Output   B - (m x n) matrix, B(:,i) = A(:,i) / v(i)
%
% (c) Panos Achlioptas 2014    http://www.stanford.edu/~optas

    assert (size(A,2) == size(v,1) || size(A,2) == size(v,2))
    assert (sum(v == 0) == 0);   
    if any(abs(v) < 1e-7)
        warning('Diving with elements that are smaller than 1e-7.')
    end
        
    N = size(A, 2);        
    B = zeros(size(A));
    for i=1:N               % TODO - parfor?
        B(:,i) = A(:,i) / v(i);
    end  
end
