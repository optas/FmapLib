function [P] = sample_random_pairs(pairs, max_elem, order, replacement)
    % Computes random pairs of integers in the range [1, max_elem]. A pair (a1, b1) always consists of different
    % integers (a1 ~= b1). The sampling can be with or without replacement and the ordering betwen two elements 
    % is taken into account only if order == 1. In all other cases (a1, b1) is considered equal to (b1, a1).
    %
    % Input:
    %        pairs       - (int) number of pairs requested
    %     
    %        max_elem    - (int) samplin will happen for pairs in range [1, max_elem].
    %     
    %        order       - (binary) determines if the order between the elements of a pair matters wrt. to the uniquness of 
    %                               a pair: if 0 (a1, b1) is the same as (b1, a1).
    %        replacement - (binary) if 1 sampling is with replacement.
    %
    % Output:
    %        P           - (pairs, 2) matrix carrying produced pairs on its lines. 
 
    
    if order == 1 || replacement == 1
        error('Not implemented yet.')
    end
        
    % Case with no replacement, and no order.
    if 2*pairs > max_elem * (max_elem-1)
        error('Without replacement you cant request more pairs than all feasible pairs.')
    end
    k = randperm(max_elem/2*(max_elem-1), pairs);       % k random integers    
    q = floor(sqrt(8*(k-1) + 1)/2 + 3/2);               % Trick that maps integers to set of all possible pairs.
    p = k - (q-1).*(q-2)/2;                                                         
    P = [p;q]';
    assert(all(P(:,1) ~= P(:,2)))
    
end