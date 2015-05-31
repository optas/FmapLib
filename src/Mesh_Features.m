classdef Mesh_Features < dynamicprops
    
    properties (GetAccess = public, SetAccess = private)
        m = [];
    end
    
    
    
    
    methods (Access = public)
        function obj = Mesh_Features(varargin)     
            if nargin == 0                
                obj.m = [];
            else
                obj.m = varargin{1};
            end
        end        
    end

    
    % Mean curvature
Max = 5; step = 5;
meanCurv = sqrt(sum((Delta*mesh.vertices).^2, 2));
f5 = zeros(mesh.nv, Max);
f6 = zeros(mesh.nv, Max);
f5(:,1) = meanCurv;
f6(:,1) = log(abs(meanCurv) + 1e-10);
for i = 2:Max
    f5(:, i) = (speye(size(Delta)) - step*Delta)\f5(:, i-1);
    f6(:, i) = (speye(size(Delta)) - step*Delta)\f6(:, i-1);
end

F = [F, f5, f6];
mask = [mask; ones(size(f5, 2), 1)*(max(mask) + 1)];
mask = [mask; ones(size(f6, 2), 1)*(max(mask) + 1)];
    
    
    


    methods (Static)
        function [mean_curv] = mean_curvature(, )
        
        end
        
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

        gauss_kernel = exp(- (( energies-evals ).^2 ) / 2*sigma^2);
        signatures   = evecs.^2 * gauss_kernel;
        scale        = sum(signatures, 1);
        signatures   = divide_columns(signatures, scale);

        assert(all(all(signatures >= 0)));
        % wks_variance = 6; % variance of the WKS gaussian (wih respect to the 
        % difference of the two first eigenvalues). For easy or precision tasks 
        % (eg. matching with only isometric deformations) you can take it smaller
    end




  end
        
end

