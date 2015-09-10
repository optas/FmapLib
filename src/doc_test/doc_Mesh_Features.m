classdef doc_Mesh_Features < dynamicprops
    % A class for creating and manipulating features on triangular Meshes. Usually, such features correspond to 
    %  functions defined on the vertices of a mesh. Examples include:
    %    * the Wave Kernel Signature,
    %    * the Heat Kernel Signature,
    %    * the Multi-scale Gaussian Curvature,
    %    * the Multi-scale Mean Curvature.
    %
    % (c) Achlioptas, Corman, Guibas  - 2015  -  http://www.fmaplib.org
    
    properties (GetAccess = public, SetAccess = private)
        % Each Mesh_Feature object has at least the following properties.
        M                   % (:class:`~Mesh.Mesh`) Mesh over which the features are computed.
        LB                  % (:class:`~Mesh.Laplace_Beltrami`) Associated LB operator.                             
        F                   % (Matrix) carrying the computed features as column vectors.
        index               % Structure holding the type of feature stored at each column of F (e.g., F(:,10) is 'wks').
    end
        
    methods (Access = public)
        
        function h = plot(obj)
            h = imagesc(obj.F);
        end
        
    end % Object's methods.
    
    methods (Static)
   
        function [mean_curv] = mean_curvature(inmesh, laplace_beltrami, smoothing_time)                                        
            % Computes the mean curvature at each vertex of a given mesh.
            % This implementation utilizes the Laplace Beltrami (LB) operator of the mesh (see Notes).
            % 
            %
            % Parameters
			% ----------
            %           inmesh            :  :class:`~Mesh.Mesh`
			%                                Input mesh with inmesh.num_vertices vertices.
            %                                
            %           laplace_beltrami  :  :class:`~Mesh.Laplace_Beltrami`
			%                                The corresponding LB of the inmesh.
            %                                             
            %           smoothing         :  (k x 1, optional) Vector 
			%                                Values corresponding to time samples for the heat diffusion smoothing TODO-P: explain more.            
            %
            % Returns
			% -------
            %           mean_curv         :  (num_vertices x k+1)
			%                                The mean curvature of each vertex.
            %                                If smoothing is applied then the first column contains
            %                                the mean curvature and the k-following columns contain
            %                                the k-smoothed versions of it.
            %
            % References
			% ----------
            %       Meyer, M. Desbrun, P. Schroder, and A. H. Barr. "Discrete
            %       Differential-Geometry Operators for Triangulated 2-Manifolds."
			%
            % .. todo:: TODO-E: Let's talk for completeness this to take as optional arguement the normals of vertices.
            
            if isprop(inmesh, 'vertex_normals')
                N = inmesh.vertex_normals;
            else
                N = Mesh.normals_of_vertices(inmesh.vertices, inmesh.triangles);
            end
            
            mean_curv = 0.5 * sum(N .* (laplace_beltrami.W * inmesh.vertices), 2);
            
            if exist('smoothing_time', 'var')
                mean_curv_smooth = Mesh_Features.heat_diffusion_smoothing(laplace_beltrami.W, mean_curv, smoothing_time);
                mean_curv        = [mean_curv, mean_curv_smooth];
            end
        end
    
    end    
end


