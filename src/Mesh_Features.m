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
%         function [mean_curv] = mean_curvature(, )                
            % Mean curvature
%             Max = 5; step = 5;
%             meanCurv = sqrt(sum((Delta*mesh.vertices).^2, 2));
%             f5 = zeros(mesh.nv, Max);
%             f6 = zeros(mesh.nv, Max);
%             f5(:,1) = meanCurv;
%             f6(:,1) = log(abs(meanCurv) + 1e-10);
%             for i = 2:Max
%                 f5(:, i) = (speye(size(Delta)) - step*Delta)\f5(:, i-1);
%                 f6(:, i) = (speye(size(Delta)) - step*Delta)\f6(:, i-1);
%             end
% 
%             F = [F, f5, f6];
%             mask = [mask; ones(size(f5, 2), 1)*(max(mask) + 1)];
%             mask = [mask; ones(size(f6, 2), 1)*(max(mask) + 1)];            
%         end
        
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
        

        function [signatures] = D2_panos(inmesh, num_bins, xcenters)
        % TODO-V
        % Input:
        %           num_bins    -  (int) Number of bins the D2 histogram will have.
        
            xcenters  =     % Optional arguement defining for each bin its center. (makes num_bins obsolete)
                
            all_dists = pdist(inmesh.vertices);  % This is too expensive(!) + it used euclidean dist and not geodesics.     \
                                                 % Do other also propose
                                                 % to use Euclideans dist?                                                  
            
            max_dist = max(max(all_dists));        % Maximum distance between two vertices.
            xcenters = linspace(max_dist/num_bins, max_dist, num_bins);  % One way to define centers.
            v_num = inmesh.num_vertices;
            

            signatures = zeros(v_num, num_bins);            
            for i = 1:v_num
                distlist = sum((repmat(inmesh.vertices(i,:), [v_num-1, 1]) - inmesh.vertices([1:i-1 i+1:end],:)).^2, 2);
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


