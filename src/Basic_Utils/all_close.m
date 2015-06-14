function [bool] = all_close(A, B, atol, rtol)
    % Checks if matrices are element-wise equal within a tolerance.
    % The tolerance values are positive, typically very small numbers. 
    % The relative difference (rtol * abs(B)) and the absolute difference atol are added together to compare against the absolute difference between A and B.    
    % TODO: If either array contains one or more NaNs, 0 is returned. Infs are treated as equal if they are in the same place and of the same sign in both arrays.
    % 
    % Input:	
    %           A, B  -  (m,n) matrices to be compared.
    %
    %           rtol  -  (float) The relative tolerance parameter (see Notes).
    % 
    %           atol  -  (float) The absolute tolerance parameter (see Notes).
    % 
    % Output:	
    %           bool -  1 iff A is close to B (see Notes).
    % Notes:
    %
    %       Iff the following equationa are both element-wise True, then allclose returns 1. 
    %               absolute(A - B) <= atol
    %               absolute(A - B) <= rtol * absolute(B))        
    %       The above equation is not symmetric in A and B, so that allclose(A, B) might be different from allclose(B, A) in some rare cases.
    
    if ~exist('atol', 'var')
        atol = 1e-08;
    end
    if ~exist('rtol', 'var')
        rtol = 1e-05;
    end
    
    assert(atol > 0 & rtol >0)
       
    abs_error =  A(:) - B(:);    

    bool = all( (abs(abs_error) < atol) & (abs(abs_error) < rtol * abs(B(:))));

end

