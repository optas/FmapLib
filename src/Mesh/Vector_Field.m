classdef Vector_Field
    % A Class representing a Vector Field defined on a Mesh object.
    % Note: 
    %       This is stub class with most functionality expected to be added by September 2015.
    %
    % (c) Achlioptas, Corman, Guibas  - 2015  -  http://www.fmaplib.org
    properties
        M = Mesh();     % Presumably we will keep a pointer to corresponding Mesh.
        spectra = [];   % Presumably we will keep a harmonic basis for it.
    end
    
    methods 
        % Class Constructor
        function obj = Vector_Field()
            if nargin == 0
                obj.M = Mesh();
            end            
        end
        
        % Ploting a Vector Field.        
        function [F] = plot(this, vertex_function)
            F = figure; 
            if ~exist('vertex_function', 'var')                
                patch('faces', this.M.triangles, 'vertices',this.M.vertices, ...
                  'FaceColor', [1,0.936,0.859], 'EdgeColor', [0.5,0.5,0.5]);
            else
                show_func(mesh,f);  % TODO-P
            end
            hold on;
            vf_fquiver(mesh,vf,mesh.nf);
            colorbar; hold off;
        end
        
    end

    

    methods (Static)
    
        function compute_spectra()            
            % compute the harmonic VF and basis functions (explain relation to Hodge decomposion)
            % TODO-E re-write dec.nabla1 with Mesh class + maybe break it down into smaller functions (some of
            % which might static functions of the Vector Field or the Mesh class.
           
            % can we add as a separate function the Laplacian over the Edges?
            % Any more intuitive name for the omega2U?           
            [vv,dd] = eigs(mesh.dec.nabla1,k,'sm');              
            hvf = omega2U(mesh, vv(:,i) / dd(i,i));    
        end
    
        
        function flow_of_vector_field()
            %I.e., the non intersecting self-map.
        end
       
        function [div_vf] = divergence_of_vector_field(vf, V, T, N, Av)           
            % Input:
            %           vf - (num_of_vertices x 3) Vector field values at each
            %           face.
            %           V  - (num_of_vertices x 3) 3D coordinates of
            %           the mesh vertices.
            %           T  - (num_of_triangles x 3) T[i] are the 3 indices
            %           corresponding to the 3 vertices of the i-th
            %           triangle. The indexing is based on -V-.                
            %           N  - (num_of_triangles x 3) N(i,:) are the
            %           coordinates of the outward normal of the i-th
            %           triangle. The length of this normal should
            %           conrerespond to its weight in the sum.
            %           Av - (num_of_vertices x 1) an array containing
            %           the areas of all the vertices.
            %
            % Output:   div_vf - (num_of_vertices x 1) Divergence of vf: one 
            %           value per vertex.

            N = N./repmat(l2_norm(N), [1, 3]);
            vf = cross(vf, N, 2);
            
            idj = [2 3 1];
            idK = [3 1 2];
            div_vf = zeros(size(V, 1), 1);
            for i = 1:3
                j = idj(i);
                k = idK(i);
                scalar_prod =  sum( vf .* ( V(T(:,j),:) - V(T(:,i),:) ) , 2 ) ;
                div_vf = div_vf + accumarray(T(:,k), scalar_prod);
            end
            
            div_vf = div_vf ./ ( 2 * Av );
        end
        
        
        
        
    end
end

