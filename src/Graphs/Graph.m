classdef Graph < dynamicprops
    % A class representing an arbitrary Graph (dyadic relation). A variety of graph related algorithms are 
    % implemented here.
        
    % TODO-P: Extend to allow self-loops.
    % (c) Achlioptas, Corman, Guibas  - 2015  -  http://www.fmaplib.org
    
    properties (GetAccess = public, SetAccess = private)
        % Basic properties that every instance of the Graph class has.        
        num_vertices   %   (int)                            -    Number of vertices.
        num_edges      %   (int)                            -    Number of edges.
        A              %   (num_vertices x num_vertices)    -    Adjacency matrix capturing the edge connectivity.
        is_directed    %   (logical)                        -    True iff the graph's eges are directed.
        name = ''      %   (string)                         -    (default = '') A string identifying the graph, 
                       %                                         e.g., 'Price_Network'.
    end
    
    methods (Access = public)
        % Class Constructor.               
        function obj = Graph(varargin)     
            if nargin == 0                
                % Construct an empty Graph.
                obj.A             = sparse([]);
                obj.is_directed   = false;
            elseif ischar(varargin{1}) && strcmp(varargin{1}(end-3:end), '.edge_list')
                % Constructor from input .edge_list file.
                edge_list_file  = varargin{1};                
                obj.is_directed = varargin{2};
                obj.A           = Graph.read_adjacency(edge_list_file, 'edge_list', obj.is_directed);
            else
                % Construct a graph from explicitly given adgjacecny.
                obj.A = varargin{1};
                obj.is_directed = varargin{2};
            end            
            % Take care of node/edge size and graph's name.
            obj.num_vertices  = size(obj.A, 1);            
            obj.num_edges     = nnz(obj.A);
            
            if obj.is_directed() == false
                if mod(obj.num_edges, 2)
                    error('An undirected graph without self-loops, must have an even number of edges.')
                end
                obj.num_edges    = obj.num_edges   ./ 2;
            end
            
            if nargin > 1 && ischar(varargin{end}) % The last argument is a string (and this string is not the first 
                obj.name = varargin{end};          % vararg which is reserved for filenames).
            else
                obj.name = '';
            end
            
            if size(obj.A, 1) ~= size(obj.A, 2)                
                error('The provided adjacency matrix is not a square matrix.')                
            end
            assert(size(obj.A, 1) == obj.num_vertices)
        end
    end
           
    methods (Static)
        
        function [A] = read_adjacency(fname, format)
            % It reads and returns the adjacency matrix as it is encoded in a file with predefined format.
            % 
            % Input:    
            %           fname  - (String) The filename of the file holding the adjacency.
            %                    
            %           format - (String - optional) The type of encoding. 
            %                    
            %                    Valid values:  'edge_list' - the simplest possible.             
            %                                 
            % 
            % Notes: 'edjge_list' is almost identical to the 'snap' graph format:
            %        Thus the 1st line is #nodes, the 2nd #edges, and an (a,b) in the next ones means a is connected 
            %        with b.            
            % TODO-P: skip when see # (i.e., comments)
            if strcmp(format, 'edge_list')                
                links = load(fname);            
                n     = links(1,1);  % Number of nodes is encoded in the 1st line.
                e     = links(2,1);  % Number of edjes is encoded in the 2st line.                
                A     = sparse(links(3:end, 1), links(3:end, 2), links(3:end, 3), n, n);
                
                entries = size(links, 1) - 2;
                assert (e == entries || 2*e == entries) % Corresponding to directed and undirected graph respectively.
                
                if sum(diag(A)) ~= 0
                    disp('self-loops exist')
                end
            end            
        end
        
        function [G] = generate(graph_type, varargin)            
            % Generator of various commonly used graphs.
            % Input:
            %           graph_type  - (String) specyfing the type of graph model to be produced.
            %                         Valid Strings:
            %                                       - 'lattice'        -  a planar lattice graph.
            %                                       - 'checkerboard'   -  a lattice graph but with diagonal node 
            %                                                             edges included.
            %
            %           varargin    - Extra parameters specifying properties of each graph model.
            %                           If graph_type is 'lattice' or 'checkerboard' then:
            %                            varargin{1}   -  (int) Numer of horizontal nodes.
            %                            varargin{2}   -  (int) Numer of vertical nodes.
            %
            % Output:
            %           G           - (Graph) resuting graph model.
            
            if strcmp(graph_type, 'lattice') || strcmp(graph_type, 'checkerboard')  % Type-checking
                if nargin ~= 3
                    error('For a lattice graph, please specify the (x,y) number of nodes.')
                end
                m = varargin{1};                                    % X-axis number of nodes.
                n = varargin{2};                                    % Y-axis number of nodes.
                if m < 1 || n < 1
                    error('The number of nodes has to be at least one.')
                end
            end
                
            if strcmp(graph_type, 'lattice')
                total_edges = 2 * (((m-1) * n) + ((n-1) * m));      % Result of from  graph Theory.
                diag_vec_1  = repmat([0; ones(n-1, 1)], m, 1);      % Horizontal connections.
                diag_vec_1  = spdiags(diag_vec_1, 1, m * n, m * n);
                diag_vec_2  = repmat([1; ones(n-1, 1)], m, 1);      % Vertical connections.
                diag_vec_2  = spdiags(diag_vec_2 , n, m * n, m * n);
                adg = diag_vec_1 + diag_vec_2;
                adg = adg + adg.';                                  % Edges are symmetric.                        
                assert(nnz(adg) == total_edges);                
                default_name = sprintf('%d_%d_lattice', m, n);
                is_directed  = false;                
                G = Graph(adg, is_directed, default_name);
            elseif strcmp(graph_type, 'checkerboard')
                 % TODO implement.
            else
                error('Not valid graph type was requested.');
            end
        end
        
        
    end
    
end

