classdef Mesh_Features < dynamicprops
    
    properties (GetAccess = public, SetAccess = private)
        m = [];
    end
    
    methods (Access = public)
        
        function obj = Mesh_Features(varargin)     
            if nargin == 0                
                obj.m = [];
            else
                obj.m = varargin{1};
            end
        end        
    end
    
    
        
end

