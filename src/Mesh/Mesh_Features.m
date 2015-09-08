classdef Mesh_Features < dynamicprops
    % A class for creating and manipulating features on triangular Meshes. Usually, such features correspond to
    % functions defined on the vertices of a mesh. Examples include:
    %     the Wave Kernel Signature,
    %     the Heat Kernel Signature,
    %     the Multi-scale Gaussian Curvature,
    %     the Multi-scale Mean Curvature.
    %
    % (c) Achlioptas, Corman, Guibas  - 2015  -  http://www.fmaplib.org
    
    properties (GetAccess = public, SetAccess = private)
        % Each Mesh_Feature object has at least the following properties.
        M;                   % (Mesh) Mesh over which the feautures are computed.
        LB;                  % (Laplace_Beltrami) Associated LB operator.                             
        F;                   % (Matrix) carrying the computed feautures as column vectors.
        index;               % Strucutre holding the type of feature stored at each column of F (e.g., F(:,10) is 'wks').
    end
        
    methods (Access = public)
        % Class Constructor.        
        function obj = Mesh_Features(inmesh, laplace_beltrami)     
            if nargin == 0                
                % Construct empty Mesh_Features.            
                obj.M  = Mesh();
                obj.LB = Laplace_Beltrami();                
            else
               obj.M   = inmesh;
               obj.LB  = laplace_beltrami;               
            end
            obj.F  = [];
            obj.index  = struct();
        end
        
        function obj = copy(this)
            % Define what is copied when a 'deep copy' is performed.
            % I.e. when obj.copy() is called.
            obj = feval(class(this)); % Instantiate new object of the same class.
                        
            % Copy all non-hidden properties (including dynamic ones).            
            p = properties(this);
            for i = 1:length(p)
                if ~isprop(obj, p{i})   % Adds all the dynamic properties.
                    obj.addprop(p{i});
                end                
                obj.(p{i}) = this.(p{i});
            end           
        end

        function newobj = keep_only(obj, features)
            % Generates a new Mesh_Features object that contains only a subset of the current features.
            % Input:
            %        option1. features - (cell array of strings) contains the names of the feature types to be
            %                            kept. E.g., features = {'wks', 'mc'}.
            %             
            %        option2. features - (n x 1) vector containing integers. The integers correspond to the
            %                            columns (features) of obj.F that will be kept.
                        
            newobj = obj.copy();            
            if iscell(features)         % e.g., {'wks', 'hks'}                
                new_feats = [];
                feat_per_category = zeros(length(features), 1);
                pos       = 0;
                
                for i = 1:length(features)
                    ind          = obj.index.(features{i});
                    lb           = ind(1);
                    rb           = ind(2);
                    f_i          = rb - lb + 1;
                    new_feats(:, pos+1:pos+f_i) = obj.F(:, lb:rb);
                    pos = pos + f_i;
                    feat_per_category(i) = f_i;
                end
                newobj.set_features(new_feats, features, feat_per_category);                
            else
                feat_cols  = sort(features);                 % Todo-P: add check on integers - handle double entries.
                feat_names = fieldnames(obj.index);
                new_feat_names    = cell(1);
                feat_per_category = [];
                
                count = 1;
                for i=1:length(feat_names)
                    ind = obj.index.(feat_names{i});                    
                    feat_type_i = sum(  bitand(feat_cols >= ind(1), feat_cols <= ind(2)) );
                    if feat_type_i > 0
                        new_feat_names{count}    = feat_names{i};
                        feat_per_category(count) = feat_type_i;
                        count = count + 1;
                    end                    
                end
                new_feats = obj.F(:, feat_cols);
                newobj.set_features(new_feats, new_feat_names, feat_per_category);                
            end
              
        end
        
        function [] = set_features(obj, new_feats, feature_names, feat_per_category)
            if sum(feat_per_category) ~= size(new_feats,2)
                error('Mismatch in size of new features and sum of features per category.')
            end
            
            obj.F = new_feats; % Todo add dimension checking.
            obj.index = [];    % Reset index.
            Mesh_Features.index_features(obj, feature_names, feat_per_category);                            
        end
        
        function [] = append_features(obj, new_feats, feature_names, feat_per_category)
        %   TODO - finish implementation
            if sum(feat_per_category) ~= size(new_feats,2)
                error('Mismatch in size of new features and sum of features per category.')
            end            
            obj.F = [obj.F new_feats];
        %             obj.index = [];    % Reset index.
        %             Mesh_Features.index_features(obj, feature_names, feat_per_category);                            
        end
        
        
        function [nf] = size(obj)
            % Returns the number of features stored by the current object.
            nf = size(obj.F, 2);
        end
        
        function [F] = project_features(obj, basis, elems, varargin)
            % Projects the features into the given basis and returns the resulting coeffients.
            % TODO-D Add info.
            %
            % Example:   project_features(Laplace_Beltrami, 10, 'normalize', 1, 'store', 0)
            options = struct('normalize', 1, 'store', 1);
            options = load_key_value_input_pairs(options, varargin{:});
            
            if options.normalize == 1                 
                F = basis.project_functions(elems, obj.F);
                F = divide_columns(F, sqrt(sum(F.^2)));
            else
                F = basis.project_functions(elems, obj.F);
            end
            
            % Storing the coefficients.
            prop_name = 'projected';
            if options.store == 1 && isprop(obj, prop_name)
                obj.projected = F;
            elseif options.store == 1 
                obj.addprop(prop_name);
                obj.projected = F;    
            end            
        end
        
        function [C] = covariance_matrix(obj, feats)
            if nargin < 2
                feats = obj.projected;               
            end
                       
            feats = feats - repmat(mean(feats, 2), 1, size(feats,2)) ;  % Center data.            
            C     = feats * feats';
        end

        function obj = compute_default_feautures(obj, neigs, wks_samples, hks_samples, mc_samples, gc_samples)
            % 'Convenience' function. 
            %  Computes any of the the implemented mesh features with default parameters. 
            %  If X_samples is zero, the the feauture type X is not computed.

            if strcmp(neigs, 'all')
                neigs = length(obj.LB.spectra.evals);
            end
            
            obj.F = Mesh_Features.default_mesh_feautures(obj.M, obj.LB, neigs, wks_samples, hks_samples, mc_samples, gc_samples);
           
            feature_names = {'wks', 'hks', 'mc', 'gc'};
            features_per_categ = [wks_samples, hks_samples, mc_samples, gc_samples];
            Mesh_Features.index_features(obj, feature_names, features_per_categ);
            assert(obj.size() == sum(features_per_categ));            
        end    
        
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
            % Input:    
            %           inmesh            -  (Mesh) Input mesh with inmesh.num_vertices 
            %                                vertices.
            %                                
            %           laplace_beltrami  -  (Laplace_Beltrami) The corresponding LB of the inmesh.
            %                                             
            %           smoothing         -  (k x 1, optional) Vector with values corresponding to time samples 
            %                                for the heat diffusion smoothing TODO-P: explain more.            
            %
            % Output:   
            %           mean_curv         -  (num_vertices x k+1) The mean curvature of each vertex.
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
            % Input:  
            %         inmesh            -  (Mesh) 
            %         
            %         laplace_beltrami  -  (Laplace_Beltrami)The corresponding LB of the inmesh.
            %                              
            %         smoothing         -  (k x 1, Optional) vector with time for the
            %                              heat diffusion processing.
            %
            % Output: 
            %         gauss_curv        -  (num_vertices x k+1) The gaussian curvature of each vertex.
            %                              If smoothing is applied then the first column contains
            %                              the mean curvature and the k-following columns contain
            %                              the k-smoothed versions of the mean.
            %             
            % Notes: See Meyer, M. Desbrun, P. Schroder, and A. H. Barr. "Discrete
            %            Differential-Geometry Operators for Triangulated 2-Manifolds."
            
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
        
%     function [signatues] = mesh_saliency(vertices, neighborhoods, mean_curvature)
%     https://github.com/daeyun/saliency-for-3d-meshes
%             params.scale       = 0.003;
%             epsilon            = params.scale * sqrt(sum((min(vertices) - max(vertices)).^2));
%             params.windowSize  = 0.33 * sqrt(sum((min(Mesh.v)-max(Mesh.v)).^2));
%             suppressedLevelSaliency = {};
%             for level = 1:5
%                 sigma = round((level+1) * epsilon, 6);
%                 F1 = gaussian_weighted_average(mean_curvature, neighborhoods, sigma, cut_off)            
%                 F2 = gaussian_weighted_average(mean_curvature, neighborhoods, 2 * sigma, cut_off)            
%                 levelSaliency = abs(F1-F2);
%                 suppressedLevelSaliency{level} = nonlinearSuppression(Mesh, levelSaliency, sigma);            
%                 signatures = sum(cat(2,suppressedLevelSaliency{:}),2);
%             end                
%         end
        
%         function [suppressedLevelSaliency] = nonlinearSuppression(Mesh, levelSaliency, sigma)
%             global cache_ params_;
%             if ~isfield(cache_, 'D');
%                 cache_.D = pdist2(Mesh.v, Mesh.v);
%             end
%             D = cache_.D;
%             if ~isfield(cache_, 'avgDist');
%                 D_ = D;
%                 D_(D==0) = Inf;
%                 cache_.avgMinDist = mean(min(D_,[],2));
%             end
%             avgMinDist = cache_.avgMinDist;
% 
%             levelSaliency = normalizeRange(levelSaliency);
% 
%             levelSaliencyMat = repmat(levelSaliency, [1 size(D, 1)]);
%             levelSaliencyMat(D > params_.windowSize) = -Inf;
% 
%             [globalMax,globalMaxI] = max(levelSaliency);
% 
%             [~,maxI] = max(levelSaliencyMat,[],2);
%             isLocalMax=(maxI==(1:size(levelSaliencyMat,1))');
%             isLocalMax(globalMaxI)=0;
% 
%             meanLocMax = mean(levelSaliency(isLocalMax));
% 
%             suppressedLevelSaliency = levelSaliency*(globalMax-meanLocMax)^2;
%         end
        

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
            % Computes the heat diffusion of a function at a given time using an implicit Euler scheme.
            % The resulting function will be smoother.
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
        

        function [signature] = D2_shape_distribution(inmesh, pairs, num_bins, true_geodesic)
            % Computes the global shape descriptor "D2" which concicely describes the histogram of pairiwse geodesic
            % distances between points that live on the mesh. This signature characterizes the whole mesh at once
            % and is not suitable for vertex-wise feature extraction.
            %
            % Input:    
            %           inmesh         -  (Mesh) Input mesh.
            % 
            %           pairs          - (int) Number of pairs to be sampled
            %                                           
            %           true_geodesic  -  (binary) TODO-P
            %
            % Output:   
            %           signatures     - (num_bins x 1) vector representing the histogram of the D2 distances.
            %
            % Notes:
            % "Shape Distributions of R. Osada, T. Funkhouser, B. Chazelle, and D. Dobkin, Princeton University"
            %
            % TODO: Sample points on the shape instead of vertices, to make methods insensitive to tessalations.
            
            pairs = sample_random_pairs(pairs, inmesh.num_vertices, 0, 0);  % A list of unique pairs of vertices.
            dists = comp_geodesics_pairs(inmesh.vertices(:,1), inmesh.vertices(:,2), inmesh.vertices(:,3), inmesh.triangles', pairs, 1);            
            xcenters   = linspace(max(dists)/num_bins, max(dists), num_bins);
            signature = hist(dists, xcenters);
            signature = signature / sum(signature);                      % Converting histogram into a pmf.            
        end
        

        function [signatures] = hist_of_euclidean_distance(inmesh, num_bins)        
        % Input:
        %           num_bins    -  (int) Number of bins the D2 histogram will have.
                
            hyper_diameter = norm(max(inmesh.vertices)- min(inmesh.vertices));                                                              
            xcenters       = linspace(hyper_diameter/num_bins, hyper_diameter, num_bins);            
            v_num          = inmesh.num_vertices;           
            
            signatures = zeros(v_num, num_bins);
            for i = 1:v_num
                distlist = sum((repmat(inmesh.vertices(i,:), [v_num-1, 1]) - inmesh.vertices([1:i-1 i+1:end],:)).^2, 2);  % TODO-P check if pdist2 is faster.
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
            %       "Geodesics in Heat: A New Approach to Computing Distance Based on Heat Flow.", 2013        

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
        
        function [F] = default_mesh_feautures(inmesh, laplace_beltrami, neigs, wks_samples, hks_samples, mc_samples, gc_samples)                       
            % Computes the wks, hks, mean_curvature and gaussian curvature with default parameter values.
            
            evals = laplace_beltrami.evals(neigs);
            evecs = laplace_beltrami.evecs(neigs);                
            if wks_samples > 1
                [energies, sigma] = Mesh_Features.energy_sample_generator('log_linear', evals(2), evals(end), wks_samples);
                wks_sig           = Mesh_Features.wave_kernel_signature(evecs(:,2:end), evals(2:end), energies, sigma);                
            else
                wks_sig           = [];
            end
            if hks_samples > 1
                heat_time         = Mesh_Features.energy_sample_generator('log_sampled', evals(2), evals(end), hks_samples);
                hks_sig           = Mesh_Features.heat_kernel_signature(evecs(:,2:end), evals(2:end), heat_time);
            else
                hks_sig           = [];
            end            
            if mc_samples > 1
                heat_time         = Mesh_Features.energy_sample_generator('log_sampled', evals(2), evals(end), mc_samples-1);
                mc_sig            = Mesh_Features.mean_curvature(inmesh, laplace_beltrami, heat_time);                    
            else
                mc_sig            = [];
            end
            if gc_samples > 1
                heat_time         = Mesh_Features.energy_sample_generator('log_sampled', evals(2), evals(end), gc_samples-1);
                gc_sig            = Mesh_Features.gaussian_curvature(inmesh, laplace_beltrami, heat_time);
            else
                gc_sig            = [];
            end
            
            F                     = [hks_sig wks_sig mc_sig gc_sig];
        end
        
        function [F] = extract_subleveled_feautures(features, levels)
            % Given as input a set of features and a set of levels representing percentiles, it extracts a new set of
            % features that is non zero only for values that are within the presribed percentiles. 
            % Note: This is a non-linear transformation of the features.
            % Input:
            % Output:
            % Example:
            
            if ~all(all(levels >= 0) && all(levels <= 100))
                error('Levels must be a vectors with values in [0,100] interval.');
            end

            [basis_size, num_feats] = size(features);            
            F = zeros(basis_size, num_feats * length(levels));                        
            pertile_per_feature = prctile(features, levels);
            bound =  0;           
            for i = 1:length(levels)
                pf = pertile_per_feature(i,:);
                mask = (features >= repmat(pf, basis_size, 1));
                F(:, bound+1:bound+num_feats) = features .* mask;   % At each iteration num_features features are created.
                bound = bound + num_feats;
            end
        end
        
        function [F] = gaussian_weighted_average(vertex_function, neighborhoods, sigma, cut_off)
            % Computes the gaussian weighted average of a function defined over mesh vertices. That is, for each vertex 
            % it averages the function values within its neighborhood, with a weight inversely proportional to the
            % distance of the neighbor vertex from the center one.  TODO-re cast/ re say.
            gauss_kernel = exp(-0.5 .* (neighborhoods.distances ./ sigma).^2);
            if nargin == 4               
                gauss_kernel(neighborhoods.distances >= cut_off) = 0;
            end           
            V = vertex_function(neighborhoods.ids);            
            F = sum(V .* gauss_kernel, 2 ) ./ sum(gauss_kernel, 2);               
        end
    
    end    

    methods (Static, Access = private)
        % Functions used only internally from other functions of this class.        
        function obj = index_features(obj, feature_types, feature_samples)                        
            
            total = length(feature_types);
            if (total ~= max(size(feature_samples)))
                error('Mismatch in size of feature types and computed feature samples.')
            end
                        
            pos = 0;            
            for i = 1:total
                feat_name = lower(feature_types{i});      % Make case insensitive.
                feat_samples = feature_samples(i);
                if feat_samples < 1
                    continue                              % Feautures of this type were not calculated.
                else                
                    obj.index.(feat_name) = [pos+1, pos+feat_samples];
                    pos = pos + feat_samples;
                end
            end
        end
        
    end
end


