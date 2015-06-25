classdef Mesh_Features < dynamicprops
    
    methods (Static)
   
        function [mean_curv] = mean_curvature(inmesh, laplace_beltrami, smoothing_time)                                        
            % Computes the mean curvature at each vertex of a given mesh.
            % This implementation utilizes the Laplace Beltrami (LB) operator of
            % the mesh (see Notes).            
            %
            % Input:    inmesh            -  (Mesh) Input mesh with inmesh.num_vertices 
            %                                vertices.
            %                                
            %           laplace_beltrami  -  (Laplace_Beltrami)The corresponding LB of the
            %                                inmesh.
            %             
            %           smoothing         -  (k x 1, optional) Vector with values corresponding to time samples 
            %                                for the heat diffusion smoothing TODO-P: explain more.            
            %
            % Output:   mean_curv         -  (num_vertices x k+1) The mean curvature of each vertex.
            %                                If smoothing is applied then the first column contains
            %                                the mean curvature and the k-following columns contain
            %                                the k-smoothed versions of it.
            %
            % Notes:
            %       Meyer, M. Desbrun, P. Schroder, and A. H. Barr. "Discrete
            %       Differential-Geometry Operators for Triangulated 2-Manifolds."
            % TODO-E: Let's talk for completeness this to take as optional arguement the normals of vertices.
            
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
        
        function [gauss_curv] = gaussian_curvature(inmesh, laplace_beltrami, smoothing_time)                                        
            % Computes the Gaussian curvature at each vertex of a given mesh.
            % (Optional) A smoothing using the heat diffusion can be done
            % as post-processing.
            %
            % Meyer, M. Desbrun, P. Schroder, and A. H. Barr. "Discrete
            % Differential-Geometry Operators for Triangulated 2-Manifolds."
            %
            % Input:  inmesh            -  (Mesh) 
            %         laplace_beltrami  -  (Laplace_Beltrami)The corresponding LB of the
            %                                inmesh.
            %         smoothing         -  (k x 1, Optional) vector with time for the
            %                              heat diffusion processing.
            %
            % Output: gauss_curv        -  (num_vertices x k+1) The gaussian curvature of each vertex.
            %                              If smoothing is applied then the first column contains
            %                              the mean curvature and the k-following columns contain
            %                              the k-smoothed versions of the mean.
            
            if isprop(inmesh, 'angles')
                angles = inmesh.angles;
            else
                if isprop(inmesh, 'edge_lengths')
                    L = inmesh.edge_lengths;
                else
                    L = Mesh.edge_length_of_triangles(inmesh.vertices, inmesh.triangles);
                end
                angles  = Mesh.angles_of_triangles(L);
            end
            
            areas = diag(laplace_beltrami.A);
            
            gauss_curv = ( 2 * pi - accumarray(inmesh.triangles(:), angles(:))) ./ areas;   % TODO-E: Is it OK that areas not normalized?
        
            if exist('smoothing_time', 'var')                                               
                gauss_curv_smooth = Mesh_Features.heat_diffusion_smoothing(laplace_beltrami.W, gauss_curv, smoothing_time);
                gauss_curv        = [gauss_curv, gauss_curv_smooth];
            end
        end
        
                                        
        function [signatures] = wave_kernel_signature(evecs, evals, energies, sigma)
            % Computes the wave kernel signature according to the spectrum of a graph 
            % derived operator (e.g., the Cotangent Laplacian). 
            % This signature was introduced in the paper of M. Aubry, U. Schlickewei and D. Cremers
            % "The Wave Kernel Signature: A Quantum Mechanical Approach To Shape Analysis" 
            % http://www.di.ens.fr/~aubry/texts/2011-wave-kernel-signature.pdf
            %
            % Usage:  [signatures] = wave_kernel_signature(evecs, evals, energies, sigma)
            %
            % Input:  evecs        - (n x k) Eigenvectors of a graph operator arranged
            %                        as columns. n is the number of nodes of the graph.
            %         evals        - (k x 1) Corresponding eigenvalues.
            %         energies     - (1 x e) Energy values over which the kernel is
            %                                evaluated.
            %         sigma        - (float, optional) Controls the variance of
            %                        the fitted gausian. Default = TODO_add.
            %
            % Output: signatures   - (n x e) Matrix with the values of the WKS for
            %                                different energies in its columns.            
            
            if(size(evals, 1) ~= size(evecs, 2))
                error('The number of eigenvalues given does not aggree with the number of eigenvectors.')
            end

            k            = size(evals, 1);                % Number of eigenvalues.
            e            = size(energies, 2);             % Number of energies.

            energies     = repmat(energies, k, 1);
            evals        = repmat(log(evals),1, e);

            gauss_kernel = exp(- (( energies-evals ).^2 ) / (2*sigma^2) );
            signatures   = evecs.^2 * gauss_kernel;
            scale        = sum(gauss_kernel, 1);            
            signatures   = divide_columns(signatures, scale);
            assert(all(all(signatures >= 0)));            
        end
        
        
        function [signatures] = heat_kernel_signature(evecs, evals, T)
            % Computes the heat kernel signature according to the spectrum of a graph operator (e.g., Laplacian).
            % The signature was introduced in the paper of J. Sun, M. Ovsjanikov and L. Guibas (2009)
            % "A Concise and Provably Informative Multi-Scale Signature-Based on Heat Diffusion."
            %
            % Usage:  [signatures] = heat_kernel_signature(evecs, evals, T)
            %
            % Input:  evecs         - (n x k) Eigenvectors of a graph operator arranged as columns. n denotes the number of nodes of the graph.
            %         evals         - (k x 1) corresponding eigenvalues
            %         T             - (1 x t) time values over which the kernel is evaluated
            %
            % Output: 
            %         signatures    - (n x t) matrix with the values of the HKS for different T in  its columns
    
            if(size(evals, 1) ~= size(evecs, 2))
                error('The number of eigenvalues given does not aggree with the number of eigenvectors.')
            end
            
            low_pass_filter = exp(- evals * T);
            signatures = evecs.^2 * low_pass_filter;
            scale = sum(low_pass_filter, 1);
            signatures = divide_columns(signatures, scale);
            assert(all(all(signatures >= 0)))
        end
    
        
        function [E, sigma] = energy_sample_generator(recipie, emin, emax, nsamples, variance)
            % TODO: add comments-explanation. variance of the WKS gaussian (wih respect to the 
            % difference of the two first eigenvalues). For easy or precision tasks 
            % (eg. matching with only isometric deformations) you can take
            % it smaller.  Yes, smaller => more distinctiveness.            
            default_variance = 5;
            sigma            = 0;
            if emin < 1e-6
                warning('The smallest eigenvalue (emin) is smaller than 1e-6.')
            end
            switch recipie                
                case 'log_linear'   % Version Used in original WKS.
                    E = linspace(log(emin), (log(emax) / 1.02), nsamples);
                    if ~exist('variance', 'var')               
                        sigma = (E(2) - E(1)) * default_variance;
                    else
                        sigma = (E(2) - E(1)) * variance;
                    end
                    
                case 'linear'
                    E = linspace(emin, emax / 1.02, nsamples);
                    delta = ( emax - emin ) / nsamples;
                    if ~exist('variance', 'var')               
                        sigma = delta * default_variance;
                    else
                        sigma = delta * variance;
                    end  
                
                case 'log_sampled'    % Version Used in original HKS.
                    % TODO-P: keep here - take care of variance param.
                    % When you first transform you range to log-scale and then
                    % sample linearly, you sample more small ranged values.
                    tmin = abs(4*log(10) / emax);
                    tmax = abs(4*log(10) / emin);                                        
                    E    = exp(linspace(log(tmin), log(tmax), nsamples));                    
                    
                otherwise
                    error('Given recipie is not recognised.')
            end
            assert(length(E) == nsamples)
        end
        
        
        function [signatures] = global_point_signature(evecs, evals)
            % Raif's embedding.
            % Requires a discrete approximation of area-weighted cotanget Laplacian . Not a graphical one.
            % DOI: "Laplace-Beltrami eigenfunctions for deformation invariant shape representation. R. Rustamov, SGP 2007."
            assert(all(evals~=0) & all(evals > 0))
            if any(abs(evals) < 1e-7) 
                warning('Eigenvalues with magnitude as small as 1e-7 are used.')
            end            
            signatures = evecs * diag((1 ./ sqrt(evals)));            
        end

               
        function [smoothed_fct] = heat_diffusion_smoothing(W, fct, diffusion_time) 
            % Computes the heat diffusion of a function for a given time using an implicit Euler scheme.
            % As a result the function will appear smoother.
            %
            % Input:  W                  -  (n x n) Matrix approximating the Laplace-Beltrami
            %                               operator
            %         fct                -  (n x 1) Vector containing the function values
            %         diffusion_time     -  (1 x k) Vector contining the diffusion times
            %
            % Output: smoothed_fct       -  (n x k) Matrix with the values of k smoothed function
            
            if size(diffusion_time, 1) > 1
                 diffusion_time = diffusion_time';
            end
            
            if size(diffusion_time, 1) ~= 1
                error('Variable "Time" should be vector');
            end
            if size(W, 1) ~= size(W, 2)
                error('Given LB matrix is not square');
            end
            if size(W, 2) ~= size(fct, 1)
                error('Uncompatible size between function and operator');
            end
            if size(fct, 2) ~= 1
                error('Variable "Function" should be a column vector');
            end
            
            diffusion_time = sort(diffusion_time);    % Time must be increasing.
            n = size(W, 2);
            k = size(diffusion_time, 2);
            
            smoothed_fct       = zeros(n, k);
            smoothed_fct(:,1)  = (speye(n, n) + (diffusion_time(1) * W )) \ fct;
            for i = 2:k
                smoothed_fct(:,i) = ( speye(n, n) + (diffusion_time(i) - diffusion_time(i-1) ) * W ) \ smoothed_fct(:, i-1) ;
            end
        end
        
        

        function [signatures] = D2(inmesh, laplace_beltrami, num_samples, isprob)
            signatures = 0;
            % TODO-V: Stub
%         The simplest shape descriptor is the "D2 shape descriptor." This consists of taking pairs of random points, 
%         computing their distances, and placing those distances into a histogram. The histogram can them be looked at 
%         as a probability distribution for how far apart points are from each other on the object. Note that because I normalized 
%         for scale by scaling down by the RMS distance of points, this distance could be unbounded in the case of outliers (i.e. having 
%         tons of points really close to the center and 1 point really far away). This case is rare, and after experimentation it seems 
%         that most points are within a distance 3 of the origin after the first step. So I simply chose to put everything further away from 
%         3 in the 3 bin. I then have 100 bins total spaced from a distance of 0 to a distance of 3, and I can take as many random samples as I want to construct the histogram in this interval.
%         NOTE: In order for the comparisons to make sense between shapes, the same number of random samples should be taken in each shape. If we wanted a different number of samples for some reason, but we still wanted to be as correct as possible, we should scale each bin down by the number of samples (thus, turning the histogram into an actual probability density function which sums to 1). I didn't do this in my program since I'm always comparing objects with the same number of random samples.

% Computes an approximation the geodesic distance of all vertices to a set.
            % This approximation comes from an heat diffusion process.
            %
            % Input:    inmesh            -  (Mesh) Input mesh with inmesh.num_vertices 
            %                                vertices.
            %                                
            %           laplace_beltrami  -  (Laplace_Beltrami)The corresponding LB of the
            %                                inmesh.
            %           num_samples       - (1 x 1) Scalar with value of number of sample pairs
            %                               we need to use to generate
            %                               histogram.
            %           isprob            - (1 x 1) Scalar taking value 1
            %                                if further scaling is to be
            %                                done.
            %         
            %
            % Output:   signatures         -  (300 x 1) histogram containing the histogram of computed distances.
            

%           The first part of the code is written to compute num_samples
%           number of unique pairs of vertices.

            num_vertices = imesh.num_vertices;
            S = zeros(num_vertices*(num_vertices-1)/2,2);
            pointer = 1;
            for i=1:num_vertices
                for j=i+1:num_vertices
                    S(pointer,:) = [i j];
                    pointer = pointer+1;
                end
            end
            sample_pairs = datasample(S,num_samples);
            
%           For every vertex, the below code computes the geo-desic distances
%           of the vertex from every other vertex, and then pulls out the
%           relevant distances alone. This is a heavy computation step.
%
            I = eye(num_vertices);
            samp_val = zeros(num_samples,1);
            pointer = 1;
            for i=1:num_vertices
                list_of_vertex_i = find(sample_pairs(:,1) == i);
                if ~isempty(list_of_vertex_i)
                    geo_dist = geodesic_distance_to_set(inmesh,laplace_beltrami,I(:,i));
                    samp_val(pointer:pointer+length(list_of_vertex_i)-1) = geo_dist(list_of_vertex_i,2);
                    pointer = pointer+list_of_vertex_i;
                end
            end
            
%           Once the distances have been computed, we bin it in the set of
%           bins with centers 0 to 300, in gaps of 0.03. Also if
%           normalizing by the number of samples has to be done, it is done as well. 

            xcenters = linspace(0,3,100);
            signatures = hist(samp_val,xcenters);
          
            if isprob==1
                signatures = signatures/num_samples;
            end
        end
        

        function [signatures] = hist_of_euclidean_distance(inmesh, num_bins)        
        % Input:
        %           num_bins    -  (int) Number of bins the D2 histogram will have.
                
            hyper_diameter = norm(max(inmesh.vertices)- min(inmesh.vertices));                                                              
            xcenters       = linspace(hyper_diameter/num_bins, hyper_diameter, num_bins);            
            v_num          = inmesh.num_vertices;           
            
            signatures = zeros(v_num, num_bins);
            for i = 1:v_num
                distlist = sum((repmat(inmesh.vertices(i,:), [v_num-1, 1]) - inmesh.vertices([1:i-1 i+1:end],:)).^2, 2);  % TODO-V check if pdist2 is faster.
                signatures(i,:) = hist(distlist, xcenters)/(v_num - 1);
            end    
        end

          function [geo_dist] = geodesic_distance_to_set(inmesh, laplace_beltrami, indicator_fct)
            % Computes an approximation the geodesic distance of all vertices to a set.
            % This approximation comes from an heat diffusion process.
            %
            % Input:    inmesh            -  (Mesh) Input mesh with inmesh.num_vertices 
            %                                vertices.
            %                                
            %           laplace_beltrami  -  (Laplace_Beltrami)The corresponding LB of the
            %                                inmesh.
            %           indicator_fct     -  (num_vertices x 1) Vector with values 1 when 
            %                                the vertex belongs to the set and 0 when ouside.            
            %
            % Output:   geo_dist         -  (num_vertices x 1) Geodesic distance between the 
            %                               i-th vertex and the set defined by indicator_fct.
            %
            % Notes:
            %       K. CRANE, C. WEISCHEDEL, M. WARDETZKY
            %       "Geodesics in Heat: A New Approach to Computing Distance
            %       Based on Heat Flow.", 2013
            %
            
            if isprop(inmesh, 'triangle_normals')
                N = inmesh.triangle_normals;                
            else
                N = Mesh.normals_of_triangles(inmesh.vertices, inmesh.triangles, 1);
            end
            
            if isprop(inmesh, 'barycentric_v_area')
                area_vertices = inmesh.barycentric_v_area;
            else           
                area_vertices = Mesh.area_of_vertices(inmesh.vertices, inmesh.triangles, 'barycentric');
            end
            
            if isprop(inmesh, 'triangle_areas')
                area_triangles = inmesh.triangle_areas;
            else           
                area_triangles = Mesh.area_of_triangles(inmesh.vertices, inmesh.triangles);
            end
            
            if isprop(inmesh, 'edge_lengths')
                L = inmesh.edge_lengths;
            else           
                L = Mesh.edge_length_of_triangles(inmesh.vertices, inmesh.triangles);
            end
            
            t              = mean(L(:))^2;
            heat_diff      = ( sparse(1:inmesh.num_vertices, 1:inmesh.num_vertices, area_vertices) + t * laplace_beltrami.W ) \ indicator_fct;
            grad_heat_diff = Mesh.gradient_of_function(heat_diff, inmesh.vertices, inmesh.triangles, N, area_triangles);
            grad_heat_diff = - grad_heat_diff ./ repmat(l2_norm(grad_heat_diff), [1, 3]);
            
            div_vf   = Mesh.divergence_of_vector_field(grad_heat_diff, inmesh.vertices, inmesh.triangles, N, area_vertices);
            geo_dist = full( laplace_beltrami.W \ ( area_vertices .* div_vf ) );
            geo_dist = geo_dist - min(geo_dist);
        end
        
  end    

end


