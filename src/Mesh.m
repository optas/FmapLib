classdef Mesh < dynamicprops
    % This class represents the triangular Mesh of a 3D surface.
    % Todo: Seperate Input Class (deal with .off, .obj etc.)
    
    properties (GetAccess = public, SetAccess = private)
        % Basic properties that every instance of the Mesh class must have.
        
        num_vertices    %   (int)               -    number of vertices of mesh.
        num_triangles   %   (int)               -    number of triangles of mesh.num_
        vertices        %   (num_vertices  x 3) -    array containing the 3D embedding of each vertex of the mesh.
        triangles       %   (num_triangles x 3) -    array containing for each face, the index of the 3 vertices
                        %                            that are part of it. The index ids are wrt. the array 'vertices'.
        name = ''       %   (string)            -    (default = '') a string identifying the mesh, e.g., 'Dragon_1'.
    end
    

    methods (Access = public)
        
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


%          function s = saveobj(obj)          
%              possible_dyn_props = {'triangle_areas', 'edge_lengths', 'angles', 'cot_laplacian'};  % make it a hidden attributre of the class
%              for i=1:length(possible_dyn_props)                                
%                 prop_name = possible_dyn_props{i};
%                 if isprop(obj, prop_name) % Object has this dynamic property.                    
%                     prop_i                      = findprop(obj, prop_name);
%                     s.dynamicprops(i).name      = prop_i.Name;
%                     s.dynamicprops(i).value     = obj.(prop_name);
%                     s.dynamicprops(i).setAccess = prop_i.SetAccess;
%                     s.dynamicprops(i).getAccess = prop_i.GetAccess;
%                 end
%              end
%          end
        
    end
    
    
    methods (Access = public)        
        % Setters and Getters.
        
        function obj = set_name(obj, name)
            obj.name = name;
        end
            
        function obj = set_triangle_areas(obj)
            % TODO: Add safeguard against property already existing. 
            obj.addprop('triangle_areas');            
            obj.triangle_areas = Mesh.area_of_triangles(obj.vertices, obj.triangles);
        end
    
        function obj = set_edge_lengths(obj)
            obj.addprop('edge_lengths');
            obj.edge_lengths = Mesh.edge_length_of_triangles(obj.vertices, obj.triangles);
        end
        
        function obj = set_triangle_angles(obj)
            obj.addprop('angles');
            if isprop(obj, 'edge_lengths')
                obj.angles = Mesh.angles_of_triangles(obj.edge_lengths);
            else
                obj.angles = Mesh.angles_of_triangles(obj.vertices, obj.triangles);
            end
        end        

        function obj = set_vertex_areas(obj, area_type)
            if ~ Mesh.is_supported_area_type(area_type)                            
                error(strcat('Area must be one of these strings: ', strjoin(Mesh.valid_area_strings(), ', '), '.'))
            end
            prop_name = strcat(area_type, '_v_area');
            obj.addprop(prop_name);
            obj.(prop_name) = Mesh.area_of_vertices(obj.vertices, obj.triangles, area_type);
        end

        function [A] = get_vertex_areas(obj, area_type)
            if ~ Mesh.is_supported_area_type(area_type)            
                error(strcat('Area must be one of these strings: ', strjoin(Mesh.valid_area_strings(), ', '), '.'))
            end
            prop_name = strcat(area_type, '_v_area');
            if isprop(obj,prop_name)
                A = obj.(prop_name);
            else
                ME = MException('Mesh:variable_not_initialized.', ...
                    'Variable %s has not been initialized via set_area_vertices() method.', strjoin('vertex_area_', area_type));
                throw(ME)
            end           
        end
        
    end
    
    
    methods (Static, Access = private)        
        function [S] = valid_area_strings()
            S = {'barycentric', 'voronoi'};
        end        
    
    end
    
    methods (Static)
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
                
                A = cross( V(T(:,1),:) - V(T(:,2),:), V(T(:,1),:) - V(T(:,3),:));           % Shall we save these normals?                                        
                A = normv(A)/2;
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
                                
                L1 = normv(V(T(:,2),:) - V(T(:,3),:));                                      % Verify this is correct.
                L2 = normv(V(T(:,1),:) - V(T(:,3),:));                                      % I would prefer edge(1,2), (1,3), (2,3)
                L3 = normv(V(T(:,1),:) - V(T(:,2),:));                               
                L  = [L1 L2 L3];

            end
            
            function [A] = angles_of_triangles(varargin)
                % Computes for each triangle the 3 angles among its edges.
                % Input:
                %   option 1:   [A] = angles_of_triangles(V, T)
                % 
                %               V  - (num_of_vertices x 3) 3D coordinates of
                %               the mesh vertices.
                %               T  - (num_of_triangles x 3) T[i] are the 3 indices
                %               corresponding to the 3 vertices of the i-th
                %               triangle. The indexing is based on -V-.                
                %
                %   option 2:   [A] = angles_of_triangles(L)
                %                               
                %               L - (num_of_triangles x 3) L[i] is a triple
                %               containing the lengths of the 3 edges
                %               corresponding to the i-th triange. 
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
                % TODO: It computes a whole matrix potentially not
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
                    Av  = sum(M,2);

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
                area_type = lower(area_type);  % Keywords are case independent.
                index = strcmpi(area_type, valid_area_types);
                bool = any(index);
            end

     end
   
end
                      
%     %%%
% % Compute weighted and unweighted cot LB
% % Note that area weights are normalized to sum to 1
% %%

%
%     [W, A] = cotLaplacian(mesh);
%     mesh.origAreaWeights = 2*A;
%     mesh.areaWeights = 2*A;
%     % % Normalize mesh area to sum to 1
%     mesh.areaWeights = mesh.areaWeights / sum(mesh.areaWeights);
%     % Using **negative** cotLaplacian
%     mesh.cotLaplace = -2*W;
%
%     %% For convenience...
%     mesh.A = spdiags(mesh.areaWeights,0,mesh.nv,mesh.nv);
%     mesh.Ai = spdiags(1./mesh.areaWeights,0,mesh.nv,mesh.nv);
%     mesh.L = mesh.cotLaplace;
%     mesh.Lw = mesh.Ai*mesh.L;
% end
