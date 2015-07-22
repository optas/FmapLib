function options = load_key_value_input_pairs(options, varargin)
    option_names = fieldnames(options); % Read the acceptable key names.
    nvargs = length(varargin);
    if round(nvargs/2) ~= nvargs/2
        error('Expecting property_name/property_value pairs.')
    end
    for pair = reshape(varargin, 2, []) % Pair is a cell {propName;propValue}.
        inp_name = lower(pair{1});      % Make case insensitive.
        if any(strcmp(inp_name, option_names))
            options.(inp_name) = pair{2};
        else
            error('%s is not a recognized parameter name.', inp_name)
        end
    end
 end


