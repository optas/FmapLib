classdef Mesh_Features < dynamicprops
    
%     properties (GetAccess = public, SetAccess = private)
%         m = [];
%     end
%     
%     
%     
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
            % This signature was introduced in the following paper: 
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

            gauss_kernel = exp(- (( energies-evals ).^2 ) / (2*sigma^2));
            signatures   = evecs.^2 * gauss_kernel;
            scale        = sum(signatures, 1);
            signatures   = divide_columns(signatures, scale);

            assert(all(all(signatures >= 0)));            
        end
        
        
        function [signatures] = HKS(evecs, evals, T)
            % Computes the heat kernel signature according to the spectrum of a graph operator (e.g., Laplacian).
            % The basic idea is introduced in the paper of Sun, J. and Ovsjanikov, M. and Guibas, L. (2009)
            % "A Concise and Provably Informative Multi-Scale Signature-Based on Heat Diffusion."
            %
            % Usage:  [signatures] = HKS(evecs, evals, T)
            %
            % Input:  evecs         - (n x k) Eigenvectors of a graph operator arranged as columns. n denotes the number of nodes of the graph.
            %         evals         - (k x 1) corresponding eigenvalues
            %         T             - (1 x t) time values over which the kernel is evaluated
            %
            % Output: signatures    - (n x t) matrix with the values of the HKS for different T in  its columns
            %
            % (c) Panos Achlioptas 2014   http://www.stanford.edu/~optas
                
                assert(size(evals, 1) == size(evecs, 2))
                signatures = evecs.^2 * exp( - evals * T);
                scale = sum(signatures, 1);
                signatures = divide_columns(signatures, scale);
                assert(all(all(signatures >= 0)))
        end
    
        function [E, sigma] = energy_sample_generator(recipie, emin, emax, nsamples, variance)
            % variance of the WKS gaussian (wih respect to the 
            % difference of the two first eigenvalues). For easy or precision tasks 
            % (eg. matching with only isometric deformations) you can take
            % it smaller.  Yes, smaller => more distinctiveness.            
            default_variance = 5;
            switch recipie                
                case 'log_linear'
                    E = linspace(log(emin), (log(emax) / 1.02), nsamples);
                    if ~exist('variance', 'var')               
                        sigma = (E(2) - E(1)) * default_variance;
                    else
                        sigma = (E(2) - E(1)) * variance;
                    end
                    
                case 'linear'
                    delta = ( emax - emin ) / nsamples;
                    if ~exist('variance', 'var')               
                        sigma = delta * default_variance;
                    else
                        sigma = delta * variance;
                    end            
            end
        end
    
        
      function [WKS] = wks_aubrey(evecs, evals, energies, sigma)   

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
      
      function [hks] = HKS_Sun(evecs, evals, A, scale)

        % INPUTS
        %  evecs:  ith each column in this matrix is the ith eigenfunction of the Laplace-Beltrami operator
        %  evals:  ith element in this vector is the ith eigenvalue of the Laplace-Beltrami operator
        %  A:      ith element in this vector is the area associated with the ith vertex
        %  scale:  if scale = true, output the scaled hks
        %          o.w. ouput the hks that is not scaled

        % OUTPUTS
        %  hks: ith row in this matrix is the heat kernel signature of the ith vertex


           %area = sum(A);
           %A = (1/area) * A;
           %evals = area * evals;
           %evecs = sqrt(area) * evecs;

           tmin = abs(4*log(10) / evals(end));
           tmax = abs(4*log(10) / evals(2));
           nstep = 100;

           stepsize = (log(tmax) - log(tmin)) / nstep;
           logts = log(tmin):stepsize:log(tmax);
           ts = exp(logts);

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


