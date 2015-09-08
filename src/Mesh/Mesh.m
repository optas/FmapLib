classdef Mesh < dynamicprops
    % This class represents a triangular Mesh of a 3D surface. It is equipped with a variety of related functions 
    % mostly of geometric nature.    
    %
    % (c) Achlioptas, Corman, Guibas  - 2015  -  http://www.fmaplib.org
    
    properties (GetAccess = public, SetAccess = private)
        % Basic properties that every instance of the Mesh class has.        
        num_vertices    %   (int)                -    Number of mesh vertices.
        num_triangles   %   (int)                -    Number of mesh triangles.
        vertices        %   (num_vertices  x 3)  -    Array containing the 3D coordinates of each vertex of the mesh.
        triangles       %   (num_triangles x 3)  -    Array containing for each triangle the index of its 3 vertices.
                        %                             The indices refer to positions in the array 'vertices'.
        name            %   (string)             -    (default = '') A string identifying the mesh, e.g., 'Happy_Dragon'.
    end
    
    methods (Access = public)
        % Class Constructor.
        % NOTE: we have assumed that .off and .obj give the vertices so that the
        % hand-right rule applies. 
        function obj = Mesh(varargin)     
            if nargin == 0                
                % Construct an empty Mesh.            
                obj.vertices      = [];
                obj.triangles     = [];                
            elseif ischar(varargin{1}) && strcmp(varargin{1}(end-3:end), '.off')
                % Constructor from input .off file.
                off_file = varargin{1};                
                [obj.vertices , obj.triangles] = Mesh_IO.read_off(off_file); 
            elseif ischar(varargin{1}) && strcmp(varargin{1}(end-3:end), '.obj')
                % Constructor from input .obj file.
                obj_file = varargin{1};                
                [obj.vertices , obj.triangles] = Mesh_IO.read_obj(obj_file); 
            else
                % Construct with explicitly given vertices/triangles.
                obj.vertices  = varargin{1};
                obj.triangles = varargin{2};
            end            
            % Take care of sizes and name.
            obj.num_vertices  = size(obj.vertices, 1);
            obj.num_triangles = size(obj.triangles, 1);                           
            
            if nargin > 1 && ischar(varargin{end}) % The last argument is a string (and this string is not the first vararg which is reserved for filenames).
                obj.name = varargin{end};
            else
                obj.name = '';
            end
        end
        
        function obj = copy(this)
            % Define what is copied when a deep copy is performed.
            % Instantiate new object of the same class.
            obj = feval(class(this));
                        
            % Copy all non-hidden properties (including dynamic ones)
            % TODO: Hidden properties?
            p = properties(this);
            for i = 1:length(p)
                if ~isprop(obj, p{i})   % Adds all the dynamic properties.
                    obj.addprop(p{i});
                end                
                obj.(p{i}) = this.(p{i});
            end           
        end
        
        function [h] = plot(this, vertex_function)
            figure; hold;
            if ~exist('vertex_function', 'var')                
                h = trisurf(this.triangles, this.vertices(:,1), this.vertices(:,2), this.vertices(:,3));                
            else                
                h = trisurf(this.triangles, this.vertices(:,1), this.vertices(:,2), this.vertices(:,3), vertex_function, 'EdgeColor', 'none');                                                

%                 vertex_function - mean(vertex_function)
%                 patch('Faces', this.triangles, 'Vertices', this.vertices, 'FaceColor', 'interp', 'FaceVertexCData', vertex_function, 'EdgeColor', 'none');                               
            end
            axis equal;
            

%             axis equal; 
%             cameratoolbar; cameratoolbar('SetCoordSys','none');
%             hold on
        end
               
        function [M] = normal_expanding_mesh(obj, normal_dist, normal_dirs)
            if any(size(normal_dirs) ~= [obj.num_vertices, 3])
                error('')
            end
            
            expanded_vertices = obj.vertices + normal_dist * normal_dirs;
            M = Mesh(expanded_vertices, obj.triangles);
        end
        
        function [B] = bounding_box(obj)
            % Returns the box (i.e., rectangular cuboid) that tangentially encloses the mesh.
            % Output:   
            %           B  -  (8 x 3) 3D coordinates of the 8 corners of the cuboid. 
            %           TODO - P: explain orientation            
            
            mins = min(obj.vertices);
            maxs = max(obj.vertices);
            B = [mins; ...
                 maxs(1), mins(2), mins(3);  maxs(1), maxs(2), mins(3);  mins(1), maxs(2), mins(3); ...                    
                 mins(1), mins(2), maxs(3); ...
                 maxs(1), mins(2), maxs(3);  maxs(1), maxs(2), maxs(3);  mins(1), maxs(2), maxs(3)]; 

            assert(max(pdist(B)) == obj.max_euclidean_distance())
            
            % TODO-E plot it.
            %             faces = [1,2,3,4; 1,2,5,6; 1,4,5,8; 2,3,6,7; 3,4,7,8; 5,6,7,8]'
        
        end
        
        function d = max_euclidean_distance(obj)
            % Returns the distance of the longest straight line connecting two vertices of the mesh.            
            d = sqrt(sum((min(obj.vertices)-max(obj.vertices)).^2));
        end
        
        
    end

    methods (Access = public)        
        
        % Setters and Getters.                
        % TODO-P: Logic to think about: want if not setted to get an exception?
        %       Are the viariables setted only via settter?
        
        function obj = set_mesh_name(obj, new_name)
                obj.name = new_name;
        end
        
        function obj = set_triangle_areas(obj)            
            Mesh.add_or_reset_property(obj, 'triangle_areas',  @Mesh.area_of_triangles, obj.vertices, obj.triangles);
        end
        
        function obj = set_edge_lengths(obj)
            Mesh.add_or_reset_property(obj, 'edge_lengths',    @Mesh.edge_length_of_triangles, obj.vertices, obj.triangles);
        end
        
        function obj = set_triangle_normals(obj)
            Mesh.add_or_reset_property(obj, 'triangle_normals', @Mesh.normals_of_triangles, obj.vertices, obj.triangles, 1);           
        end 
        
        function obj = set_vertex_normals(obj)
            Mesh.add_or_reset_property(obj, 'vertex_normals',   @Mesh.normals_of_vertices, obj.vertices, obj.triangles)
        end
        
        function obj = set_triangle_angles(obj)
            if isprop(obj, 'edge_lengths')
                Mesh.add_or_reset_property(obj, 'angles', @Mesh.angles_of_triangles, obj.edge_lengths);
            else
                Mesh.add_or_reset_property(obj, 'angles', @Mesh.angles_of_triangles, obj.vertices, obj.triangles);
            end
        end
                
        function obj = set_default_vertex_areas(obj, area_type)
            % Computes and stores on the mesh the vertex areas according to a rule described by 'area_type'.
            % Also it makes for this mesh the default rule for finding vertex areas to be the one given by
            % 'area_type'. This is useful for avoiding re-specifying the area_type in other function calls that
            % operate on the mesh.
            % Input:
            %        area_type - (String) 'barycentric' or 'voronoi'.
            %
            
            if ~ Mesh.is_supported_area_type(area_type)            
                error(strcat('Area must be one of these strings: ', strjoin(Mesh.valid_area_strings(), ', '), '.'))
            end
            propname = 'default_v_area';
            if isprop(obj, propname) && ~strcmp(obj.(propname), area_type)                
                warning ('Changing default vertex areas from %s ', obj.(propname), 'to %s', area_type, '.')
                obj.(propname) = area_type;
            else
                obj.addprop(propname);
                obj.(propname) = area_type;
            end            
            prop_name = strcat(area_type, '_v_area');
            Mesh.add_or_reset_property(obj, prop_name, @Mesh.area_of_vertices, obj.vertices, obj.triangles, area_type);            
        end
        
        function [area_type] = get_default_vertex_areas(obj)
            propname = 'default_v_area';
            if ~isprop(obj, propname)
                error('Default value has not been set.')                
            else
                area_type = obj.(propname);
            end            
        end
        
        function obj = set_vertex_areas(obj, area_type)
            if ~ Mesh.is_supported_area_type(area_type)                            
                error(strcat('Area must be one of these strings: ', strjoin(Mesh.valid_area_strings(), ', '), '.'))
            end            
            prop_name = strcat(area_type, '_v_area');
            Mesh.add_or_reset_property(obj, prop_name, @Mesh.area_of_vertices, obj.vertices, obj.triangles, area_type);
        end
        
        function [A] = get_vertex_areas(obj, area_type)
            if ~exist('area_type', 'var')
                area_type = obj.get_default_vertex_areas();
            end
            
            if ~ Mesh.is_supported_area_type(area_type)            
                error(strcat('Area must be one of these strings: ', strjoin(Mesh.valid_area_strings(), ', '), '.'))
            end

            prop_name = strcat(area_type, '_v_area');
            if isprop(obj,prop_name)
                A = obj.(prop_name);
            else
                ME = MException('Mesh:variable_not_initialized.', ...
                    'Variable %s has not been initialized via set_area_vertices() method.', ['vertex_area_', area_type]);   
                throw(ME)
            end           
        end
               
    end
    
    methods (Static)
        
        function [N] = normals_of_vertices(V, T, weights)
            % Computes the normalized outward normal at each vertex by adding the weighted normals of each triangle a 
            % vertex is adjacent to. The weights that are used are the actual area of the triangle a normal comes from.
            %             
            % Input:
            %           V   -   (num_of_vertices x 3) 3D coordinates of
            %                   the mesh vertices.
            %             
            %           T   -   (num_of_triangles x 3) T[i] are the 3 indices
            %                   corresponding to the 3 vertices of the i-th
            %                   triangle. The indexing is based on -V-.                
            %           
            %           weights - (optional, num_of_triangles, x 1) These are positive values that expres how much
            %           each corresponding traingle will contribute on the vertex normal. Default = . TODO-P
            %             
            % Output:   N   -   (num_of_vertices x 3) an array containing
            %                   the normalized outward normals of all the vertices.

            
            if exist('weights', 'var')
                if any(size(weights) ~= [size(T,1), 1]) || any(weights < 0)
                    error('Expecteing a positive weight per triangle.')
                end                    
                N = Mesh.normals_of_triangles(V, T, 1);
                N = repmat(weights, [1, 3]) .* N;
            else
                N = Mesh.normals_of_triangles(V, T);
            end
            
            N = [accumarray(T(:), repmat(N(:,1), [3,1])) , accumarray(T(:), repmat(N(:,2) , [3,1])), accumarray(T(:), repmat(N(:,3), [3,1]))];
            N = N ./ repmat(l2_norm(N), [1, 3]);
        end
        
        function [N] = normals_of_triangles(V, T, normalize)
            % Computes the outward normal vector of each triangle of a given mesh.            
            %
            % Input:
            %           V           -   (num_of_vertices x 3) 3D coordinates of the mesh vertices.
            %                           
            %           T           -   (num_of_triangles x 3) T[i] are the 3 indices corresponding to the 3 vertices of                         
            %                           the i-th triangle. The indexing is based on -V-.                
            %             
            %           normalize   -   (int, optional) if 1, the normals
            %                           will be normalized to have unit lenght.
            %                           
            % Output:   
            %           N           -   (num_of_triangles x 3) an array containing
            %                           the outward normals of all the triangles.
            
            % TODO-P: We do an assumption on the order of the vertices.
            N = cross( V(T(:,1),:) - V(T(:,2),:), V(T(:,1),:) - V(T(:,3),:));
            if (exist('normalize', 'var') && normalize == 1)
                N = divide_columns(N', l2_norm(N))';                    
            end
        end
        
        function [A] = area_of_triangles(V, T)
            % Computes the area of each triangle, in a triangular mesh.
            % Input:
            %           V  - (num_of_vertices x 3) 3D coordinates of
            %           the mesh vertices.
            %           T  - (num_of_triangles x 3) T[i] are the 3 indices
            %           corresponding to the 3 vertices of the i-th
            %           triangle. The indexing is based on -V-.                
            % 
            % Output:   A - (num_of_triangles x 1) an array containint
            %           the areas of all the triangles.

            A = Mesh.normals_of_triangles(V, T);
            A = l2_norm(A)/2;
        end  

        function [L] = edge_length_of_triangles(V, T)
            % Computes the length of each edge, of each triangle in the underlying triangular mesh.
            % Input:
            %           V  - (num_of_vertices x 3) 3D coordinates of
            %           the mesh vertices.
            %           T  - (num_of_triangles x 3) T[i] are the 3 indices
            %           corresponding to the 3 vertices of the i-th
            %           triangle. The indexing is based on -V-.                
            %
            % Output:                
            %           L - (num_of_triangles x 3) L[i] is a triple
            %           containing the lengths of the 3 edges
            %           corresponding to the i-th triange. The
            %           enumeration of the triangles is the same at in
            %           -T- and the order in which the edges are
            %           computes is (V2, V3), (V1, V3) (V1, V2). I.e.
            %           L[i][2] is the edge lenght between the 1st
            %           vertex and the third vertex of the i-th triangle.

            L1 = l2_norm(V(T(:,2),:) - V(T(:,3),:));                                      % TODO-E: Verify this is correct.
            L2 = l2_norm(V(T(:,1),:) - V(T(:,3),:));                                      % I would prefer edge(1,2), (1,3), (2,3)
            L3 = l2_norm(V(T(:,1),:) - V(T(:,2),:));                               
            L  = [L1 L2 L3];

        end

        function [A] = angles_of_triangles(varargin)
            % Computes for each triangle the 3 angles among its edges.
            % Input:
            %   option 1:   [A] = angles_of_triangles(V, T)
            % 
            %               V   - (num_of_vertices x 3) 3D coordinates of
            %                     the mesh vertices.
            %               T   - (num_of_triangles x 3) T[i] are the 3 indices
            %                     corresponding to the 3 vertices of the i-th
            %                     triangle. The indexing is based on -V-.                
            %
            %   option 2:   [A] = angles_of_triangles(L)
            %                               
            %               L - (num_of_triangles x 3) L[i] is a triple
            %                   containing the lengths of the 3 edges
            %                   corresponding to the i-th triange. 
            %
            % Output:
            %                 

            if nargin == 2
                V = varargin{1};   % Vertices.
                T = varargin{2};   % Triangles.                
                L  = Mesh.edge_length_of_triangles(V, T);
            else
                L = varargin{1};
            end

            L1 = L(:, 1); L2 = L(:, 2); L3 = L(:, 3);
            A1 = (L2.^2 + L3.^2 - L1.^2) ./ (2. * L2 .* L3);
            A2 = (L1.^2 + L3.^2 - L2.^2) ./ (2 .* L1 .* L3);
            A3 = (L1.^2 + L2.^2 - L3.^2) ./ (2 .* L1 .* L2);
            A  = [A1, A2, A3];
            A  = acos(A);            
        end
        
        function [Av] = area_of_vertices(V, T, area_type)
            % TODO-E: It computes a whole matrix potentially not
            % necessary!
            if ~exist('area_type', 'var') || strcmp(area_type, 'barycentric')    % Default is barycentric.
                Ar = Mesh.area_of_triangles(V, T);
                nv = size(V, 1);
                I   = [T(:,1);T(:,2);T(:,3)];
                J   = [T(:,2);T(:,3);T(:,1)];
                Mij = 1/12*[Ar; Ar; Ar];
                Mji = Mij;
                Mii = 1/6*[Ar; Ar; Ar];
                In  = [I;J;I];
                Jn  = [J;I;I];
                Mn  = [Mij;Mji;Mii];
                M   = sparse(In, Jn, Mn, nv, nv);                   
                Av  = full(sum(M,2));
            else
                % Compute based on Voronoi.
                % (following the algorithm decribed in 
                % Discrete Differential-Geometry Operators for Triangulated 2-Manifolds, 
                % Mark Meyer, Mathieu Desbrun, Peter Schroeder and Alan H. Barr, VisMath 2002).
                error('Not implemented yet.');

            end
        end

        function [bool] = is_supported_area_type(area_type)
            % Returns 1 iff the area_type (string) corresponds to a
            % method of creating vertex areas, that is supported by Mesh class.                
            valid_area_types = Mesh.valid_area_strings();            
            index = strcmpi(area_type, valid_area_types);  % Keywords are case independent.
            bool = any(index);
        end
        

        function [df] = gradient_of_function(f, V, T, N, A)            
            % Computes the gradient of a function defined on the vertices of a mesh. The function is assumed to be
            % interpolated linearly (via the barycentric basis functions) at each triangle. For more information see: 
            %           'Polygon Mesh Processing, Botsch et al., 1st edition - 2010,  page 44.'
            %             
            % Input:
            %           f  - (num_of_vertices x 1)  A vector representing a function with a value at every vertex of a 
            %                                       mesh. f[i] is the function's value on vertex -i-.
            %
            %           V  - (num_of_vertices x 3)  3D coordinates of the mesh vertices.            
            %             
            %           T  - (num_of_triangles x 3) T[i] are the 3 indices corresponding to the 3 vertices of the i-th
            %                                       triangle. The indices refer to vertices of -V-.                
            %
            %           N  - (num_of_triangles x 3) N(i,:) are the coordinates of the outward normal of the i-th
            %                                       triangle. 
            %             
            %           A  - (num_of_triangles x 1) an array containing the areas of all the triangles.
            %           
            % Output:   
            %           df - (num_of_triangles x 3) The gradient of the function i.e., one vector per triangle.
            
            idj = [2 3 1];
            idK = [3 1 2];
            df = zeros(size(T, 1), 3);
            for i = 1:3
                j = idj(i);
                k = idK(i);
                df = df + repmat(f(T(:,k)), [1,3]) .* ( V(T(:,j),:) - V(T(:,i),:) );
            end
            
            N = N./repmat(l2_norm(N), [1, 3]);
            df = cross(N, df, 2) ./ repmat(2 * A, [1, 3]);
        end
        
    end
    
    
    
    
    methods (Static, Access = private)
        % Functions used only internally from other functions of this class.
        
        function [S] = valid_area_strings()
            % We implement the following types of vertex-area constructions.            
            S = {'barycentric', 'voronoi'};
        end            
        
        function [] = add_or_reset_property(obj, propname, setter, varargin)
            % Check if the dynamic property already exists. In this case it
            % only updates it via the setter and the varargin. Otherwise, 
            % it first adds it on the object.            
            if isprop(obj, propname)
                obj.(propname) = setter(varargin{:});
            else
                obj.addprop(propname);
                obj.(propname) = setter(varargin{:});
            end
        end        
    end 
   
end
