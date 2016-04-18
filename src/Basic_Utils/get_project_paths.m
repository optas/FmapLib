function [data_path, code_path] = get_project_paths(project_name)
% Computes the personalized file paths for the data and external code libraries, assosiated with a specific project.
% Input:
%           project_name - (str) A string specifying the project's name.
% Output:   
%           data_path    - (str) File path to project's data repository.
%           code_path    - (str) File path to project's (external) code repository.


    name = get_computer_name();
        
    data_path = '';
    code_path = '';
    
    if strcmp(name, 'optasMacPro')                                                          % Panos's Space.
        if strcmp(project_name, 'FmapLib')
            data_path = '/Users/optas/Documents/Git_Repos/FmapLib/data/';
            code_path = '/Users/optas/Documents/Git_Repos/FmapLib/src/External_Code/';
        
        elseif strcmp(project_name, 'Joint_Image_Analysis')
            data_path = '/Users/optas/Dropbox/Matlab_projects/Joint_Image_Analysis/Data/';
            code_path = '/Users/optas/Dropbox/Matlab_projects/External_Packages/';
                
        elseif strcmp(project_name, 'Shape_Classification')
            data_path = '/Users/optas/Dropbox/Matlab_Projects/Shape_Classification/Data/';
            code_path = '/Users/optas/Dropbox/Matlab_Projects/Shape_Classification/src/';
        end
        
    elseif strcmp(name(1:5), 'orion') && strcmp(name(end-13:end), '.stanford.edu')          % Orion server farm on Stanford.        
        if strcmp(project_name, 'FmapLib')
            data_path = '/orions3-zfs/projects/optas/Matlab_Projects/FmapLib/data/';
            code_path = 'orions3-zfs/projects/optas/Matlab_Projects/FmapLib/src/';
        
        elseif strcmp(project_name, 'Joint_Image_Analysis')
            data_path = '/orions3-zfs/projects/optas/Data/Joint_Image_Analysis/';
            code_path = '/orions3-zfs/projects/optas/Matlab_Projects/External_Matlab_Code/';             
        
        end
        
    elseif strcmp(name, 'Etienne-HP')                                                       % Etienne's space.        
        data_path = 'C:\Users\Etienne\Desktop\GitHubProj\FmapLib\data\';
        code_path = 'C:\Users\Etienne\Desktop\GitHubProj\FmapLib\src\External_Code\';
        
    else
        
        error('Unknown computer.');
        
    end
        
    if isempty(data_path)
        error('Unknown Project.');
    end
    
end