classdef Image_Graph < Graph
    % A class representing a Graph object associated with an Image object.
    %
    % (c) Achlioptas, Corman, Guibas  - 2015  -  http://www.fmaplib.org
    
    properties (GetAccess = public, SetAccess = private)
        I;          %  (Image) - Underlying Image object.                
    end
        
    methods (Access = public)
        % Class Constructor.               
        function obj = Image_Graph(varargin)
            % Set up super-class (Graph) arguements.
            if nargin == 0
                super_args = cell();
            else
                [h, w, ~] = size(varargin{1}.CData);
                G         = Graph.generate(varargin{2}, h, w);
                super_args{1} = G.A;       % Adjacency matrix.               % TODO-P Add construct from graph in Graph.
                super_args{2} = false;     % is_directed attribute.
                super_args{3} = G.name;    % Graph's name.
            end            
            obj@Graph(super_args{:})
            
            if nargin == 0                            % Set Image.
                obj.I = Image();
            else
                obj.I = varargin{1};             
            end            
        end
    end
     

end