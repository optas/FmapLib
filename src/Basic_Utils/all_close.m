function [bool] = all_close(A, B, atol, rtol)
    % Checks if matrices are element-wise equal within a tolerance. The tolerance values are positive and typically very 
    % small numbers. The relative difference and the absolute difference are both accounted to compare against the 
    % absolute difference between A and B.    
    %     
    % TODO-V: If an array contains a NaN or Inf then the other array must have at the same location also a
    % NaN or an Inf (of the same sign) to return 1. 
    %     
    % Input:	
    %           A, B  -  (m,n) matrices to be compared.
    %
    %           rtol  -  (float, or +inf) The relative tolerance parameter (see Notes).
    % 
    %           atol  -  (float, or +inf) The absolute tolerance parameter (see Notes).
    % 
    % Output:	
    %           bool -  1 iff A is close to B (see Notes).
    % Notes:
    %
    %       Iff the following equations hold both element-wise, then allclose returns 1. 
    %     
    %               1. absolute(A - B) <= atol
    %               2. absolute(A - B) <= rtol * absolute(B))        
    %
    %       If atol = +inf, then only the second equation is taken into account.
    %       If rtol = +inf, then only the fisrt equation is taken into account.    
    %       Also, observe that equation  2, is not symmetric in A and B, so that all_close(A, B) might be different from 
    %       allclose(B, A) in some rare cases.
    
    if ~exist('atol', 'var')
        atol = 1e-08;
    end
    if ~exist('rtol', 'var')
        rtol = 1e-05;
    end
    
    if atol < 0 || rtol <0
        error('Tolerance parameters must be positive real numbers.')
    end
    
    abs_difs = abs(A(:) - B(:));
    
    if atol == +Inf
        bool = all(abs_difs < (rtol * abs(B(:)))) ;
    elseif rtol == + Inf
        bool = all(abs_difs < atol);
    else
        bool = all(abs_difs < atol) && all(abs_difs < (rtol * abs(B(:))));
    end    

end

