classdef Graph < dynamicprops
    % A class representing an arbitrary Graph (dyadic relation). A variety of graph related algorithms are 
    % implemented here.
        
    % TODO-P: Extend to allow self-loops.
    % (c) Achlioptas, Corman, Guibas  - 2015  -  http://www.fmaplib.org
    
    properties (GetAccess = public, SetAccess = protected)
        % Basic properties that every instance of the Graph class has.        
        num_vertices   %   (int)                            -    Number of vertices.
        num_edges      %   (int)                            -    Number of edges.
        A              %   (num_vertices x num_vertices)    -    Adjacency matrix capturing the edge connectivity.
        is_directed    %   (logical)                        -    True iff the graph's eges are directed.
        name           %   (string)                         -    (default = '') A string identifying the graph, 
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
            
            if nargin > 1 && ischar(varargin{end})  % The last argument is a string (and this string is not the first 
                obj.name = varargin{end};           % vararg which is reserved for filenames).
            else
                obj.name = '';
            end
            
            if size(obj.A, 1) ~= size(obj.A, 2)                
                error('The provided adjacency matrix is not a square matrix.')                
            end
            if any(obj.A < 0)
                error('The provided adjacency matrix has negative entries but the edges must be non negative.')                
            end
            if ~ obj.is_directed && ~issymmetric(obj.A)
                error('The provided adjacency matrix is not symmetric but an undirected graph was expected.')
            end            
            assert(size(obj.A, 1) == obj.num_vertices)            
        end
        
        function obj = set_name(obj, new_name)
            obj.name = new_name;
        end
        
        function obj = copy(this)
            % Defines what is copied when a deep copy is performed.                        
            obj = feval(class(this));  % Instantiate new object of the same class.
                        
            % Copy all non-hidden properties (including dynamic ones) % TODO: Hidden properties
            p = properties(this);
            for i = 1:length(p)
                if ~isprop(obj, p{i})   % Adds all the dynamic properties.
                    obj.addprop(p{i});
                end                
                obj.(p{i}) = this.(p{i});
            end           
        end
        
        function [node_from, node_to, weight] = all_edges(obj)
            % Computes all the edges between every pair of nodes. If the graph is undirected it returns each edge once 
            % and the. TODO-P add documentation
            %
            if obj.is_directed
                [node_from, node_to, weight] = find(obj.A);                   
            else
                [node_from, node_to, weight] = find(triu(obj.A));     % Each edge is returned once since there is not direction.                  
            end
        end
        
        function [W] = edge_weight(obj, from, to)
            if ~ all(size(from) == size(to))
                error('"from" and "to" must have the same size, since they define vertices connected by an edge.')
            end
            if length(from) > 1
                ind = sub2ind(size(obj.A), from, to);
                W   = obj.A(ind);
            else
                W = obj.A(from, to);    
            end
        end

        function [N] = neighbors(obj, vertices)            
            % Finds the neighbor vertices (edge-connected) of all input vertices.
            %
            % Parameters
            %
            % Returns
            %
            N = cell(length(vertices), 1);              % TODO-P take care of directionality. 
            for i = 1:length(vertices)
                    N{i} = find(obj.A(i,:));                               
            end                       
        end
        
        function [I] = incidence_matrix(obj)
            % Computes the oriented incidence matrix of the underlying graph. Each column of the incidence matrix corresponds to 
            % and edge of the graph and has exactly two non zero entries. If the oriented... add_content
            %
            % I - (num_vertices x num_edges) matrix. A column corresponds to an edge 
                        
            I = sparse([], [], [], obj.num_vertices, obj.num_edges, 2*obj.num_edges);            
            [node_from, node_to] = obj.all_edges();                        
            for i=1:obj.num_edges                       % TODO- speed up.
                I(node_from(i), i) = 1; 
                I(node_to(i),   i) = -1;
            end
            
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
            %                                       - 'clique'         -  a fully connected graph.
            %
            %           varargin    - Extra parameters specifying properties of each graph model.
            %                           If graph_type is 'lattice' or 'checkerboard' then:
            %                            varargin{1}   -  (int) Numer of horizontal nodes.
            %                            varargin{2}   -  (int) Numer of vertical nodes.
            %
            % Output:
            %           G           - (Graph) resuting graph model.
            
            
            if strcmp(graph_type, 'lattice') || strcmp(graph_type, 'checkerboard') || strcmp(graph_type, 'r_radius_connected')  % Type-checking
                if nargin < 3
                    error('For a lattice graph, please specify the (x,y) number of nodes.')
                end
                m = varargin{1};                                    % X-axis number of nodes.
                n = varargin{2};                                    % Y-axis number of nodes.
                if m < 1 || n < 1
                    error('The number of nodes has to be at least one.')
                end
            end
                
            if strcmp(graph_type, 'clique')
                n = varargin{1};
                directed  = false;                
                adj = ones(n,n) - diag(ones(n,1)) ; 
                G = Graph(adj, directed, sprintf('%d_clique', n));
            elseif strcmp(graph_type, 'lattice')
                total_edges = 2 * (((m-1) * n) + ((n-1) * m));      % Result of from graph Theory.
                diag_vec_1  = repmat([0; ones(n-1, 1)], m, 1);      % Horizontal connections.
                diag_vec_1  = spdiags(diag_vec_1, 1, m * n, m * n);
                diag_vec_2  = repmat([1; ones(n-1, 1)], m, 1);      % Vertical connections.
                diag_vec_2  = spdiags(diag_vec_2 , n, m * n, m * n);
                adj = diag_vec_1 + diag_vec_2;
                adj = adj + adj.';                                  % Edges are symmetric.                        
                assert(nnz(adj) == total_edges);                
                default_name = sprintf('%d_%d_lattice', m, n);
                directed  = false;                
                G = Graph(adj, directed, default_name);
            elseif strcmp(graph_type, 'checkerboard')
                adj = simple_graphs('checkerboard_dense', m, n);                              
%                 diag_vec_1  = repmat([0; ones(n-1, 1)], m, 1);      % Horizontal connections. 
%                 diag_vec_1  = spdiags(diag_vec_1, 1, m * n, m * n);
%                 diag_vec_2  = repmat([1; ones(n-1, 1)], m, 1);      % Vertical connections.
%                 diag_vec_2  = spdiags(diag_vec_2 , n, m * n, m * n);
%                 diag_vec_3  = [0; diag_vec_1(1 : (n * (m-1)))];     % Anti-diagonal connections.
%                 diag_vec_3  = spdiags(diag_vec_3, n-1, m * n, m * n);
%                 diag_vec_4  = diag_vec_3(2 : end-1);                % Diagonal connections.
%                 diag_vec_4  = spdiags(diag_vec_4, n+1, m * n, m * n);                
%                 adj = diag_vec_1 + diag_vec_2 + diag_vec_3 + diag_vec_4;
%                 adj = adj + adj.';                        
                default_name = sprintf('%d_%d_checkerboard', m, n);
                directed  = false;                
                G = Graph(adj, directed, default_name);        
            elseif strcmp(graph_type, 'r_radius_connected')
                radius = varargin{3};
                rows = m ;
                columns = n;
                all_starts = {}; all_ends = {}; all_vals = {};
                image = ones(rows, columns);
                for i1=1:2*radius + 1
                    i = i1 - radius - 1;
                    width = ceil(sqrt(radius^2 - i^2+1));
                    starts = [];
                    ends = [];
                    vals = [];
                    for j = -width:width
                        if i == 0 && j == 0
                            continue;
                        end
                        section = image(1+max(i,0):end+min(i,0), 1+max(j,0):end + min(j,0));
                        shift   = zeros(size(image,1),size(image,2));
                        shift(1 - min(i,0):end-max(i,0),1-min(j,0):end-max(j,0)) = section;
                        where = (image == shift);
                        weights = zeros(rows, columns);
                        dist=sqrt(i^2 + j^2);
                        weights(where == 1) = dist;
%                         weights = weights';
                        start_indices = find(weights)';
%                         end_indices = start_indices + i*size(image,2) + j;
                        end_indices = start_indices + j*size(image,1) + i;
                        val         = weights(weights~=0);

                        starts(end+1:end+length(start_indices)) = start_indices;
                        ends(end+1:end+length(end_indices)) = end_indices;
                        vals(end+1:end+length(val)) = val';
                    end
                    all_starts{i1} = starts;
                    all_ends{i1}   = ends;
                    all_vals{i1}   = vals;
                end
                all_starts = cell2mat(all_starts);
                all_ends   = cell2mat(all_ends);                
                all_vals   = cell2mat(all_vals);
                adj        = sparse(all_starts, all_ends, all_vals);
                assert(all_close(adj, adj', 0.0001, +Inf));
                default_name = sprintf('%d_%d_%d_radius_connected', m, n, radius);
                directed  = false;                
                G = Graph(adj, directed, default_name);                              
            else
                
                
                error('Not valid graph type was requested.');
            end
        end
        
        
    end
    
end

