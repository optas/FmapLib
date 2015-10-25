classdef Graph < dynamicprops
    % A class representing an arbitrary Graph (dyadic relation). A variety of graph related algorithms are 
    % implemented here.
    %
    % notes: Adjacency (sparse) matrix representation of a graph. For directed ones, A(i,j) is the weight of the edge
    %        from  i to j.
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
            if any(any(obj.A) < 0)
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

        function [N] = out_neighbors(obj, vertices)            
            % Finds the out-neighbor vertices (edge-connected) of all input vertices.
            %
            % Parameters
            %
            % Returns
            %                        
            if numel(vertices) == 1
                N  = find(obj.A(vertices,:));                                              
            else
                N = cell(length(vertices), 1);
                for i = 1:length(vertices)                
                        N{i} = find(obj.A(vertices(i),:));                               
                end   
            end
        end
                
        function [N, W] = in_neighbors(obj, vertices)                        
            if numel(vertices) == 1
                N = find(obj.A(:,vertices));                                              % Is column indexing in sparse as fast as row?
                W = full(obj.A(N, vertices));
            else
                N = cell(length(vertices), 1);
                W = cell(length(vertices), 1);
                for i = 1:length(vertices)                
                        N{i} = find(obj.A(:,vertices(i)));                               
                        W{i} = full(obj.A(N{i}, vertices(i)));
                end   
            end
        end
               
        function obj = add_edge(obj, from, to, weight)
            if weight <= 0
                error('Edge weight must be positive.');
            end            
            obj.A(from, to) = weight;
            if ~ obj.is_directed
                obj.A(to, from) = weight;
            end            
        end
        
        function remove_edge(obj, from, to)
            if obj.A(from, to) == 0
                warning('Request to remove a non-existing edge.')
            end
            obj.A(from, to) = 0;
            if ~ obj.is_directed
                obj.A(to, from) = 0;
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
            elseif strcmp(graph_type, 'checkerboard')       % TODO - generate directly sparse checkerboard.
                adj = simple_graphs('checkerboard_dense', m, n);
                adj = sparse(adj); 
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
                assert(all(nonzeros(adj) > 0));
                default_name = sprintf('%d_%d_%d_radius_connected', m, n, radius);
                directed  = false;                
                G = Graph(adj, directed, default_name);                              
            else                               
                error('Not valid graph type was requested.');
            end
        end
        
        function [A] = knn_to_adjacency(neighbors, weights, direction)            
            % Converts neighborhood data into the adjacency matrix of the underlying directed/weighted graph.
            % 
            % Input:                
            %        neighbors  - (N x K) neighbors(i,j) is j-th neighbor of the i-th node.
            %                             
            %        weights    - (N x K) weights(i,j) is the weight of the (directed) edge between i to j.
            %                                   
            %        direction  - (optioal, String) 'in' or 'out'. If 'in' then weights(i,j) will be given in an edge
            %                     that points towards i. Otherwise, towards j. Default = 'in'.
            %
            % Output:   
            %           A       - (N x N) sparse adjacency matrix.      
            if any(any(weights < 0))
                error('Non negative weights for an adjacency matrix are not supported.')
            end
            [n, k]      = size(neighbors);            
            temp        = repmat(1:n, k, 1)';
            i = temp(:);
            j = neighbors(:);
            v = weights(:);
            A = sparse(i, j, v, n, n)';
            assert(nnz(A) == sum(sum(weights>0)))
            if exist('direction', 'var') && strcmp(direction, 'out')
                A = A';                
            end
        end
        
        
    end
    
end

