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
    
    if strcmp(name(1:end-1), 'optasMacPro')                                                 % Panos's Space.
        if strcmp(project_name, 'FmapLib')
            data_path = '/Users/optas/Documents/Git_Repos/FmapLib/data/';
            code_path = '/Users/optas/Documents/Git_Repos/FmapLib/src/External_Code/';
        
        elseif strcmp(project_name, 'ImageJointUnderstanding')
            data_path = '/Users/optas/Dropbox/With_others/Zimo_Peter_Panos/Joint_Image_Understanding/Data/';
            code_path = '/Users/optas/Dropbox/matlab_projects/External_Packages/';
                
        else strcmp(project_name, 'Shape_Classification')
            data_path = '/Users/optas/Dropbox/Matlab_Projects/Shape_Classification/Data/';
            code_path = '/Users/optas/Dropbox/Matlab_Projects/Shape_Classification/src/';
        end
        
    elseif strcmp(name(1:end-1), 'orionp.stanford.edu')
        if strcmp(project_name, 'ImageJointUnderstanding')
            data_path = '/orions3-zfs/projects/optas/Matlab_projects/Data/ImageJointUnderstanding/';
            code_path = '';
        end
        
    elseif strcmp(name(1:end-1), 'Etienne-HP')                                              % Etienne's space.        
        data_path = 'C:\Users\Etienne\Desktop\GitHubProj\FmapLib\data\';
        code_path = 'C:\Users\Etienne\Desktop\GitHubProj\FmapLib\src\External_Code\';
    
    elseif strcmp(name(1:end-1), 'Vignesh_Comp_Name')                                       % Vignesh's space.
        data_path = '';
        code_path = '';        
    
    else
        error('Unknown computer.');
    end
        
    if isempty(data_path)
        error('Unknown Project.');
    end
    
end