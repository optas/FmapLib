classdef doc_testing < handle & my.super.Class
    % General 'idea' about Class's functionality.
    %
    % Summary 
    % .. autosummary::
    % get_x
    %
    % Contents
                                                                      
   
    properties (SetAccess = private)        
        spectra    %:attr: A struct carrying the eigenvalues and eigenvectors of the basis.                        
    end
    
    methods
        function h = doc_testing(x)
            % lala
            %
            h.x = x
        end
        
        function x = get_x(obj)
        % how is this displayed?
            x = obj.x
        end
    end
    
    methods (Static)
        function w = my_static_function(z)

        % A static function in 
        %
        % Args:
        %       z (int): The path of the file to wrap        
        %
        % Returns:
        %       w (Laplce_Beltrami): A buffered writable file descriptor

            w = z
        end
    end
end