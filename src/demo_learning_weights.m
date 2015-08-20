%%  A Script researching Learning Weights of probe functions of Functional Maps.
    clr;
    gitdir;    
    cd FmapLib/src
%%
    Cats = Mesh_Collection('tosca_cats', '../data/input/tosca_cats/');    
    num_eigs =50; area_type = 'barycentric';
    
    Cats.compute_laplace_beltrami_basis(num_eigs, area_type);

    Cats.compute_default_feautures();    
           
    cat_names      = {'cat0', 'cat1', 'cat2', 'cat3'};
    
    pairs          = {'cat0', 'cat1'; ... 
                      'cat0', 'cat2'; ...
                      'cat0', 'cat3'; ...
                      };
                  
    groundtruth{1} = (1:Cats.meshes('cat0').num_vertices)';
    groundtruth{2} = (1:Cats.meshes('cat0').num_vertices)';
    groundtruth{3} = (1:Cats.meshes('cat0').num_vertices)';
    [gt_fmaps]     = Cats.compute_ground_truth_fmaps(pairs, groundtruth);
        
%%
    raw_feats   = Cats.get_property_of_meshes(cat_names, 'raw_features');    
    proj_feats  = Cats.project_features(cat_names, 'all', raw_feats);       %TODO-P, take of 'all' in a smoother way.
    
    norm_mesh1 = sqrt(sum(proj_feats{1}.^2));
    for i = 1:length(cat_names)
        proj_feats{i} = divide_columns(proj_feats{i}, norm_mesh1);
    end
   
    lbs = Cats.get_property_of_meshes(cat_names, 'lb_basis');
    
    % W    = the way you do the regularization, W contains all matrices (of same size as the fmaps). These matrixec
    %        give you the penalty of the the LB regularization        (this work for the laplacian basis.)
    
    W = cell(length(cat_names)-1, 1);  % Regularizer of Laplacian
    for i = 1:length(cat_names)-1
        W{i} = abs(repmat(abs(lbs{i+1}.spectra.evals) , [1, num_eigs]) - repmat(abs(lbs{1}.spectra.evals)' , [num_eigs, 1])) + 1;
    end
 
    %%
    num_features = size(raw_feats{1}, 2);
    mask         = 1:num_features;          % TODO-P make this default.
    lambda       = 1e4;                     % Weight of regularizaiton (same for all maps).
    epsilon      = 1e-3;                    % Smoothening of (Nuclear norm) to be able to differentiate it.
    init_guess = ones(num_features, 1); 
    
    % The reference shape is the first in proj_feats
    funObj       = @(x) oracle(x, gt_fmaps, proj_feats{1}, proj_feats(2:end), W, lambda, mask, 'nuclear', epsilon);
    % Try auto-grad and see difference.
    
    
    options.MaxIter     = 100;
    options.MaxFunEvals = options.MaxIter;
    %  second arg below is is initial solution.

    [x, fmin, exitflag] = minFunc(funObj, init_guess, options);
   %%
    options = optimoptions(@fminunc,'GradObj','on', 'Display', 'iter', 'Algorithm', 'quasi-newton');
    x = fminunc(funObj, init_guess, options);
    
    %%
    
    Sol = zeros(nbEigen, nbEigen, Max);
    A = zeros(nbEigen*Max, nbEigen);
    for i = 1:Max
        Sol(:,:,i) = findFunMap(eye(nbEigen), FRef, F(:,:,i), D, W(:,:,i), alpha);
        A(nbEigen*(i-1) + 1:nbEigen*i,1:end) = Sol(:,:,i) - map(:,:,i);
    end

    [~, ~, V] = svd(A);

    Ysub = fliplr(V);
    Y = Ysub(:, 1:20);