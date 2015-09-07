function [data_path, code_path] = get_project_paths(project_name)
% Computes the personalized file paths for the data and external code libraries, assosiated with a specific project.
% Input:
%           project_name - (str) A string specifying the project's name.
% Output:   
%           data_path    - (str) File path to project's data repository.
%           code_path    - (str) File path to project's (external) code repository.

    if ismac
        [~, name] = system('scutil --get ComputerName');
    else
        [~, name] = system('hostname');
    end
    
    data_path = '';
    code_path = '';
    
    if strcmp(name(1:end-1), 'optasMacPro')        
        if strcmp(project_name, 'FmapLib')
            data_path = '/Users/optas/Documents/Git_Repos/FmapLib/data/';
            code_path = '/Users/optas/Documents/Git_Repos/FmapLib/src/Extrernal_Code';
        elseif strcmp(project_name, 'ImageJointUnderstanding')
            data_path = '/Users/optas/Dropbox/with_others/Zimo - Peter - Panos/Joint_Image_Understanding/Data/';
            code_path = '/Users/optas/Dropbox/matlab_projects/External_Packages/';            
        end
        
    elseif strcmp(name(1:end-1), 'Ettiene_pc_name')
        data_path = 'Etienne add here';
        code_path = 'Etienne add here';
    
    else
        error('Unknown computer.');
    end
        
    if isempty(data_path)
        error('Unknown Project');
    end
    
end