classdef doc_testing < handle & my.super.Class
    % General 'idea' about Class's functionality.
    %
    % Summary 
    % .. autosummary::
    %
    % Contents
                                                                      
   
    properties (SetAccess = private)        
        spectra    %:attr: A struct carrying the eigenvalues and eigenvectors of the basis.                        
    end
    
    methods
        function h = doc_testing(x)
            h.x = x
        end
        function x = get.x(obj)
        % how is this displayed?
            x = obj.x
        end
    end
    
    methods (Static)
        function w = my_static_function(z)
        % A static function in 
        %
        % :param int z: input z
        % :return: w
        % :rtype: Laplace_Beltrami

            w = z
        end
    end
end