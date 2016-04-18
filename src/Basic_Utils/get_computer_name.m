function name = get_computer_name()
    if ismac
        [~, name] = system('scutil --get ComputerName');
    else
        [~, name] = system('hostname');
    end
    name = name(1:end-1); % Remove trailing new line character.
end