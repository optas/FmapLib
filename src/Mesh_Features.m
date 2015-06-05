classdef Mesh_Features < dynamicprops
    
%     properties (GetAccess = public, SetAccess = private)
%         m = [];
%     end
%     

%     methods (Access = public)
%         function obj = Mesh_Features(varargin)     
%             if nargin == 0                
%                 obj.m = [];
%             else
%                 obj.m = varargin{1};
%             end
%         end        
%     end


    methods (Static)
        
        function [smoothed_func] = laplacian_smoothing(W, Function, Time)
            % Computes the heat diffusion of a function for a given time.
            % As a result the function will appear smoother.
            %
            % Input:  W - (n x n) Matrix approximating the Laplace-Beltrami
            %         operator
            %         Function - (n x 1) Vector containing the function
            %         values
            %         Time     - (1 x k) Vector contining the diffusion
            %         times
            %
            % Output: smoothed_func - (n x k) Matrix with the values of k 
            %         smoothed function
            
            if size(Time, 1) > 1
                 Time = Time';
            end
            assert(size(Time, 1) == 1, 'Variable "Time" should be vector');
            assert(size(W, 1) == size(W, 2), 'Given LB matrix is not square');
            assert(size(W, 2) == size(Function, 1), 'Uncompatible size between function and operator');
            assert(size(Function, 2) == 1, 'Variable "Function" should be a column vector');
            
            n = size(W, 2);
            k = size(Time, 2);
            
            smoothed_func = zeros(n, k);
            for i = 1:k
                smoothed_func(:,i) = ( speye(n, n) + Time(i) * W) \ Function ;
            end
        end
        
        function [mean_curv] = mean_curvature(inmesh, laplace_beltrami, smoothing_time)                                        
            % Computes the mean curvature at each vertices of the mesh
            % using the Laplace-Beltrami operator end vertices positions.
            % (Optional) A smoothing using the heat diffusion can be done
            % as post-processing.
            %
            % http://en.wikipedia.org/wiki/Mean_curvature
            %
            % Input:  inmesh            - (Mesh class) 
            %         laplace_beltrami  - (Laplace_Beltrami class) contains
            %         the Laplace-Beltrami operator defined for inmesh.
            %         smoothing         - (k x 1) vector with time for the
            %         heat diffusion processing. (Optional)
            %
            % Output: mean_curv         - (nv x k) Smoothed mean curvature
            %         at each vertices. If smoothing_time is not given, k = 1
            %         and no smoothing is applied.
            
            % If not given compute LB operator
            if ~exist('laplace_beltrami', 'var')
                laplace_beltrami = Laplace_Beltrami(inmesh);
            end
            
            try % Retrieve the vertex normals or compute them.
                N = inmesh.vertex_normal;
            catch           
                try
                    tn = inmesh.triangle_normal;
                catch
                    tn = Mesh.normals_of_triangles(inmesh.vertices, inmesh.triangles);
                end
                N  = Mesh.normals_of_vertices(inmesh.triangles, tn);
            end
            
            mean_curv = 0.5 * sum(N .* (laplace_beltrami.W * inmesh.vertices), 2);
        
            if exist('smoothing_time', 'var')
                mean_curv = Mesh_Features.laplacian_smoothing(laplace_beltrami.W, mean_curv, smoothing_time);
            end
        end
        
        function [gauss_curv] = gaussian_curvature(inmesh, smoothing_time)                                        
            % Computes the gauss curvature at each vertices of the mesh.
            % (Optional) A smoothing using the heat diffusion can be done
            % as post-processing.
            %
            % Meyer, M. Desbrun, P. Schroder, and A. H. Barr. ?Discrete Differential-Geometry Operators for Triangulated 2-Manifolds.?
            %
            % Input:  inmesh            - (Mesh class) 
            %         smoothing         - (k x 1) vector with time for the
            %         heat diffusion processing. (Optional)
            %
            % Output: gauss_curv        - (nv x k) Smoothed gauss curvature
            %         at each vertices. If smoothing_time is not given, k = 1
            %         and no smoothing is applied.
            
            % If not given compute LB operator
            if ~exist('laplace_beltrami', 'var')
                laplace_beltrami = Laplace_Beltrami(inmesh);
            end
            
            try % Retrieve the angles or compute them.
                angles = inmesh.angles;
            catch           
                try
                    L = inmesh.edge_lengths;
                catch
                    L = Mesh.edge_length_of_triangles(inmesh.vertice, inmesh.triangles);
                end
                angles  = Mesh.angles_of_triangles(L);
            end
            
            try % Retrieve the areas or compute them.
                areas = inmesh.barycentric_v_area;
            catch           
                areas = Mesh.area_of_vertices(inmesh.vertices, inmesh.triangles, 'barycentric');
            end
            
            gauss_curv = ( 2 * pi - accumarray(inmesh.triangles(:), angles(:))) ./ areas;
        
            if exist('smoothing_time', 'var')
                gauss_curv = Mesh_Features.laplacian_smoothing(laplace_beltrami.W, gauss_curv, smoothing_time);
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
            %
            % (c) Panos Achlioptas 2014   http://www.stanford.edu/~optas

            assert(size(evals, 1) == size(evecs, 2));

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
            % Output: signatures    - (n x t) matrix with the values of the HKS for different T in  its columns
            %
            % (c) Panos Achlioptas 2014   http://www.stanford.edu/~optas
                
            assert(size(evals, 1) == size(evecs, 2))
            low_pass_filter = exp(- evals * T);
            signatures = evecs.^2 * low_pass_filter;
            scale = sum(low_pass_filter, 1);     %TODO-E
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

        function [signatures] = D2()
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
            

      function [WKS] = wks_Aubry(evecs, evals, energies, sigma)   

        % Added by Panos to make the function work
        N            = length(energies);
        num_vertices = size(evecs, 1);
        log_E = log(abs(evals))';        
        e = energies;
        % End of Panos addition
                
        WKS = zeros(num_vertices, N);
        C = zeros(1,N); %weights used for the normalization of f_E

        for i = 1:N
            WKS(:, i) = sum(evecs.^2.* ...
                       repmat( exp((-(e(i) - log_E).^2) ./ (2*sigma.^2)), num_vertices, 1),2);
            
            C(i) = sum(exp((-(e(i)-log_E).^2)/(2*sigma.^2)));
        end

        % normalize WKS
        WKS(:,:) = WKS(:,:)./repmat(C,num_vertices,1);        
      end
      
      function [hks] = hks_Sun(evecs, evals, A, ts, scale)

        % INPUTS
        %  evecs:  ith each column in this matrix is the ith eigenfunction of the Laplace-Beltrami operator
        %  evals:  ith element in this vector is the ith eigenvalue of the Laplace-Beltrami operator
        %  A:      ith element in this vector is the area associated with the ith vertex
        %  scale:  if scale = true, output the scaled hks
        %          o.w. ouput the hks that is not scaled
        %  ts :    time slices to evaluate the HKS

        % OUTPUTS
        %  hks: ith row in this matrix is the heat kernel signature of the ith vertex


           %area = sum(A);
           %A = (1/area) * A;
           %evals = area * evals;
           %evecs = sqrt(area) * evecs;

           if scale == true, 
              hks = abs( evecs(:, 2:end) ).^2 * exp( ( abs(evals(2)) - abs(evals(2:end)) )  * ts);
              Am = sparse([1:length(A)], [1:length(A)], A);
              colsum = sum(Am*hks);
              scale = 1.0./ colsum; 
              scalem = sparse([1:length(scale)], [1:length(scale)], scale);
              hks = hks * scalem;
           else
              hks = abs( evecs(:, 2:end) ).^2 * exp( - abs(evals(2:end)) * ts);

           end
      end

  end    

end


