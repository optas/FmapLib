function index = string_keys_to_ints(dict)
    %   Assigns every key of a dictionary to an integer in the range [1, length(dict.keys)].
    %   Input: 
    %               dict    -   (containers.Map) a dictionary with string type keys.
    %   Output:
    %               index   -   (containers.Map) a dictionary with the same keys as the dict, but with integer values.
    %
    %   Precondtition: the keys of the dict are strings.

    index = containers.Map;
    i = 1;            
    for key = dict.keys
        index(key{:}) = i;
        i = i + 1;                
    end
    assert(i-1 == length(dict.keys))
end