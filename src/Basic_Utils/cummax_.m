function y = cummax_(x)
    % Computes the cummulative maximum element at every dimnension of
    % a vector.
    %
    % Input:
    %           x      -  (n x 1) Vector.                                                                                      
    %
    % Output:   y      -  (n x 1) Vector such that y(i) = max(x(1:i))
    %                          
    % Note: In versions of MATLAB15a and onwards this function is natively implemented.    
    y = x;
    for i = 2:length(x)
        y(i) = max(y(i), y(i-1));
    end
end