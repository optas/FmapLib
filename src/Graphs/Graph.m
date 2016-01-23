classdef Graph < dynamicprops
    % A class representing an arbitrary Graph (dyadic relation). A variety of graph related algorithms are 
    % implemented here.
    %
    % notes: Adjacency (sparse) matrix representation of a graph. For directed ones, A(i,j) is the weight of the edge
    %        from  i to j.
    % TODO-P: Extend to allow self-loops.
    % (c) Achlioptas, Corman, Guibas  - 2015  -  http://www.fmaplib.org
    % See also: http://strategic.mit.edu/downloads.php?page=matlab_networks
    
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
            
            elseif IS.single_number(varargin{1})
                n = varargin{1};
                obj.A           = sparse([],[],[], n, n);                
                obj.is_directed = varargin{2};
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
        
        function obj = set_directed(obj, choice)
            if choice
                if ~ obj.is_directed
                    assert(nnz(obj.A) == obj.num_edges * 2)
                    obj.num_edges   = obj.num_edges * 2;
                    obj.is_directed = true;
                end                
            else
                if obj.is_directed
                    obj.A = (obj.A + obj.A') ./ 2;
                    obj.is_directed = false;
                    obj.num_edges   = nnz(obj.A);
                    self_loops      = obj.self_loops();                    
                    assert(~mod(obj.num_edges-length(self_loops), 2));
                    obj.num_edges  = obj.num_edges ./ 2;
                end                
            end
        end

        function D = in_degrees(obj)            
            D = full(sum(obj.A, 1));
        end
        
        function D = out_degrees(obj)
            D = full(sum(obj.A, 2));
        end
        
        function B = is_weighted(obj)
            B = ~IS.binary(obj.A);
        end
                
        function S = self_loops(obj)           
           S = find(diag(obj.A) ~= 0);                      
        end
                
        function obj = remove_self_loops(obj)
           N = obj.num_vertices;
           self_loop_num = sum(diag(obj.A) ~= 0);           
           obj.A(1: N+1 :N^2) = 0;
           obj.num_edges = obj.num_edges - self_loop_num;
        end
        
        function obj = drop_weights(obj)
            obj.A(obj.A > 0) = 1;
        end
        
        function b = has_self_loops(obj)
            b = true;
            if all(diag(obj.A) == 0)
                b = false;
            end
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
                W   = full(obj.A(ind));
            else
                W = full(obj.A(from, to));
            end
        end
        
        function [N, W] = out_neighbors(obj, vertices)            
            % Finds the out-neighbor vertices (edge-connected) of all input vertices.
            %
            % Parameters
            %
            % Returns
            %                        
            if numel(vertices) == 1
                N  = find(obj.A(vertices,:)');                                              
                W  = full(obj.A(vertices, N)');
            else
                N = cell(length(vertices), 1);
                W = cell(length(vertices), 1);
                for i = 1:length(vertices)                
                        N{i} = find(obj.A(vertices(i), :))';                               
                        W{i} = full(obj.A(vertices(i), N{i})');
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
                        N{i} = find(obj.A(:, vertices(i)));                               
                        W{i} = full(obj.A(N{i}, vertices(i)));
                end
            end
        end
               
        function obj = add_edge(obj, from, to, weight)
            if weight <= 0
                error('Edge weight must be positive.');
            end
            
            if obj.A(from, to) == 0;    %  This is a new edge.
                obj.num_edges = obj.num_edges + 1;
            end
            
            obj.A(from, to) = weight;
            if ~ obj.is_directed
                obj.A(to, from) = weight;
            end            
        end
        
        function remove_edge(obj, from, to)
            if obj.A(from, to) == 0
                warning('Request to remove a non-existing edge.')
                return
            end
            obj.A(from, to) = 0;
            if ~ obj.is_directed
                obj.A(to, from) = 0;
            end
            obj.num_edges = obj.num_edges - 1;
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
        
        function obj = flip_edges(obj)
            obj.A = obj.A';             
        end
                
        function [T] = triangle_edges_in_ego_network(obj, ego)
            if ~ obj.is_directed
                T = [];
                direct_neighb = obj.in_neighbors(ego);       % Direct neighbors.                
                direct_neighb = setdiff(direct_neighb, ego); % Disreguard self-loops.                
                
                for i = 1:length(direct_neighb)                    
                    i_neighbs = obj.in_neighbors(direct_neighb(i));
                    i_neighbs = setdiff(i_neighbs, direct_neighb(i)); % Disreguard self-loops.
                    common = intersect(i_neighbs, direct_neighb);
                    t = length(common);
                    if t ~= 0
                        T(end+1:end+t,:) = [repmat(direct_neighb(i), t, 1) common];
                    end
                end
                T = sort(T, 2);
                T = unique(T, 'rows');
            else
                error('NIY.')
            end
        end
        
        function t = num_triangles(obj)            
            % TODO-P does it work for undirected graphs?
            Au = obj.A;
            n  = size(Au, 1); 
            Au(1:n+1:n*n) = 0;  % Remove self-loops.
            Au(Au ~= 0) = 1; 
            t = trace(Au^3) / 6;
        end
        
        function [T] = triangles(obj)
            N = obj.num_vertices;            
            T = [];
            for i = 1:N
                t_in_i = obj.triangle_edges_in_ego_network(i);    
                t      = size(t_in_i, 1);
                if t > 1                    
                    T(end+1:end+t,:) = [t_in_i repmat(i, t, 1)];
                elseif t == 1
                    T(end+1:end+t,:) = [t_in_i repmat(i, t, 1)'];
                end
            end                    
            T = unique(sort(T, 2), 'rows');
            assert(obj.num_triangles == size(T,1))
        end
       
        function [S, C] = connected_components(obj, weak)
            % Computes connected components (cc) of the underlying graph. This is a wrapper function of the 
            % 'graphconncomp' of the bioinformatics toolbox.
            % 
            %   Input:                 
            %           weak - (logical) If true, then the directionality of the edges is being disregarded for the
            %                  computation of the connectedness.
            %
            %   Output: 
            %           S    - (int) Number of connected components.
            %           C    - (obj.num_nodes x 1) Vector encoding the cc each node belongs to.
            %   More info at: Bioinformatics Toolbox - 'graphconncomp'.
            if ~ exist('weak', 'var')
                weak = false;
            end
            [S, C] = graphconncomp(obj.A, 'Directed', obj.is_directed, 'WEAK', weak);
        end
        
        function G = giant_connected_component(obj)
            [S, C] = connected_components(obj);
            if S > 1
                [num_nodes, ids] = hist(C, unique(C));
                [num_nodes, pos] = max(num_nodes);          
                id               = ids(pos);                % CC with maximum counter (num_nodes)
                ids              = find(id == C);
                G                = obj.subgraph('nodes', ids);
                assert(G.num_vertices == num_nodes);
            else
               G = obj.copy;
            end                    
        end
        
        function G = subgraph(obj, type_of_list, list)
            switch type_of_list
                case 'nodes'
                    if length(list) ~= length(unique(list))
                        error('Repeated nodes were given.')
                    end
                    list = sort(list);
                    new_adjacency   = obj.A(list, list);
                    G = Graph(new_adjacency, obj.is_directed);                   % TODO- More fancy book-keeping + dyn_properties.
                otherwise
                    error('NIY.')
            end
                
                
                
            
        end
            
        
%         function [E] = ego_network(obj, node)
% 
%             if ~ obj.is_directed
%                 [direct_neighb, weights] = obj.in_neighbors(node);   % Direct neighbors.
%                 direct_size = length(direct_neighb);
%                 center = repmat(node, direct_size, 1);
%                 E = [center, direct_neighb];
%                 
%                 for i = 1:direct_size
%                     n = direct_neighb(node);
%                     other_nodes = setdiff(direct_neighb, n);
%                     weights = obj.edges(repmat(n, length(other_nodes)), other_nodes);
%                     to_add = find(weights);
%                     if ~isempty(to_add)
%                         adjacency(n, other_nodes(to_add)) = weights(to_add);
%                     end
%                 end
%                 
% 
%             else
%                 error('Not implemented yet.');
%             end
%            
%             
%         end
        
        
%         function [T] = triangles(obj)
%             
%         end
    
        function [PR, converged, iter] = page_rank(obj, dumping, max_iter)
            % TODO - check diagonal is zero. 
            % Graph is connected etc.            
            if ~exist('max_iter', 'var')
                max_iter = 500;
            end
            if ~ exist('dumping', 'var')
                dumping = 0.85;
            end
            
            N  = obj.num_vertices;
            PR = (1/N) * ones(N, 1);           
            
            if any(abs(sum(obj.A, 2) - 1) > 1e-6)    % The out-edges of some node do not sum to 1.
                A  = divide_columns(obj.A' , sum(obj.A, 2));                            
            else
                display('Already, out-degree normalized.')
                A = obj.A';
            end
            iter = 1;
            delta_PR = true;                 
            converged = true;
            tele_transport = ((1-dumping) / N) .* ones(N,1);
            % Iterate until norm of PR changes less than 1e-6                
            while delta_PR
                prev_PR = PR;               
                PR       = (dumping .* A * PR) + tele_transport;
                delta_PR = any(PR - prev_PR > 1e-6);                                                           
                iter = iter + 1;
                if iter > max_iter
                    converged = false;
                    break
                end                    
            end            
            
        end
%        power_iter_rank = pinv((eye(length(A)) - d*A'))*(((1-d)/N)*ones(length(A),1));  % Alternative 1.
%        if dumping == 1: [~, PR] = eigs(obj.A').

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
                n        = varargin{1};
                directed = false;                
                adj      = ones(n,n) - diag(ones(n,1)) ; 
                G        = Graph(adj, directed, sprintf('%d_clique', n));
            
            elseif strcmp(graph_type, 'lattice')
                adj = lattice(m, n);
                default_name = sprintf('%d_%d_lattice', m, n);
                directed  = false;                
                G = Graph(adj, directed, default_name);
            
            elseif strcmp(graph_type, 'checkerboard')
                adj = checker_board(m, n);
                adj = sparse(adj); 
                default_name = sprintf('%d_%d_checkerboard', m, n);
                directed  = false;                
                G = Graph(adj, directed, default_name);        
            
            elseif strcmp(graph_type, 'star')
                center       = varargin{1};
                neighb       = varargin{2};
                direction    = varargin{3};
                if ~any(strcmpi(direction, {'in', 'out', 'no'}))
                    error('Direction input must be "in" or "out" or "no".')
                end
                adj          = star(center, neighb, direction);
                default_name = sprintf('star_centered_at_%d', center);                
                G = Graph(adj, ~strcmp(direction, 'no'), default_name);                
            
            elseif strcmp(graph_type, 'bipartite_fc')
                left_nodes  = varargin{1};
                right_nodes = varargin{2};                                  
                num_nodes   = left_nodes  + right_nodes;                
                i           = vec(repmat(1:left_nodes, right_nodes, 1));
                j           = vec(repmat(left_nodes+1:left_nodes+right_nodes, right_nodes, 1)');
                assert(length(i) == left_nodes * right_nodes );
                adj          = sparse(i, j, ones(length(i), 1), num_nodes,  num_nodes);                    
                default_name = sprintf('bipartite_fc_%d_%d', left_nodes, right_nodes);                
                directed     = false;
                G = Graph(adj+adj', directed, default_name);
                
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
            
            function adj = checker_board(r, c)               % TODO - generate directly sparse checkerboard.
                diagVec1 = repmat([ones(c-1,1); 0],r,1);
                diagVec1 = diagVec1(1:end-1);           
                diagVec2 = [0; diagVec1(1:(c*(r-1)))];                                                          
                diagVec3 = ones(c*(r-1),1);                                                                     
                diagVec4 = diagVec2(2:end-1);                                                                   
                adj = diag(diagVec1,1) + diag(diagVec2,c-1) + diag(diagVec3,c) + diag(diagVec4,c+1);
                adj = adj + adj.';                            
            end
          
            function adj = lattice(m, n)
                total_edges = 2 * (((m-1) * n) + ((n-1) * m));      % Result of from graph Theory.
                diag_vec_1  = repmat([0; ones(n-1, 1)], m, 1);      % Horizontal connections.
                diag_vec_1  = spdiags(diag_vec_1, 1, m * n, m * n);
                diag_vec_2  = repmat([1; ones(n-1, 1)], m, 1);      % Vertical connections.
                diag_vec_2  = spdiags(diag_vec_2 , n, m * n, m * n);
                adj = diag_vec_1 + diag_vec_2;
                adj = adj + adj.';                                  % Edges are symmetric.                        
                assert(nnz(adj) == total_edges);                
            end
            
            function adj = star(center, neighbors, directionality)
                if IS.single_number(neighbors)
                    num_nodes = neighbors+1;                    
                    neighbors = 1:neighbors;
                    if any(center == neighbors)
                        neighbors = setdiff(neighbors, center);
                        neighbors(end+1) = neighbors(end) + 1;
                    end
                    center = repmat(center, num_nodes-1, 1);                    
                    if strcmp(directionality, 'out')
                        adj = sparse(center, neighbors, ones(num_nodes-1,1), num_nodes,  num_nodes);                             
                    elseif strcmp(directionality, 'in')
                        adj = sparse(neighbors, center, ones(num_nodes-1,1), num_nodes,  num_nodes);                    
                    else
                        adj = sparse(center, neighbors, ones(num_nodes-1,1), num_nodes,  num_nodes);                             
                        adj = adj + adj';
                    end
                else
                    error('Not implemented yet');
                end
               
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

