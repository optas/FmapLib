%%  A Script researching Ensemble Procedures for Functional Maps.
    clr;
    gitdir;
    cd FmapLib/src        

%% Load two Meshes and compute their LBs.
    meshfile       = '../data/input/tosca_small/michael11.off';
    mesh1          = Mesh(meshfile, 'mike11');
    mesh1.set_default_vertex_areas('barycentric');    
    LB1            = Laplace_Beltrami(mesh1);       
    feats1         = Mesh_Features(mesh1, LB1);    
    
    meshfile       = '../data/input/tosca_small/michael12.off';
    mesh2          = Mesh(meshfile, 'mike12');
    mesh2.set_default_vertex_areas('barycentric');
    LB2            = Laplace_Beltrami(mesh2);
    feats2         = Mesh_Features(mesh2, LB2);

%%  Compute Mesh features that will be used to produce Fmaps.
    neigs          = 100;                           % Eigenvectors that will be used in computing wks/hks.
    wks_samples    = 100;
    hks_samples    = 100;    
    mc_samples     = 50; 
    gc_samples     = 50;
    feats1.compute_default_feautures(neigs, wks_samples, hks_samples, mc_samples, gc_samples);
	feats2.compute_default_feautures(neigs, wks_samples, hks_samples, mc_samples, gc_samples);    
    
%% Compute (whp.) non-corresponding features and make an Fmap.
    random_samples = 150;
    deltas1        = Functional_Map.random_delta_functions(mesh1.barycentric_v_area, random_samples);
    rnd_feats_1    = Mesh_Features(mesh1, LB1);
    rnd_feats_1.set_features(deltas1, {'rand_deltas'}, random_samples);
    
    deltas2        = Functional_Map.random_delta_functions(mesh2.barycentric_v_area, random_samples);
    rnd_feats_2    = Mesh_Features(mesh2, LB2);
    rnd_feats_2.set_features(deltas2, {'rand_deltas'}, random_samples);  

%%  Compute some Fmaps.
    fmap_method    = 'frobenius_square';
    lambda         = 20;                                                    
    all_map        = Functional_Map(LB1, LB2);
    all_map.compute_f_map(fmap_method, neigs, neigs, feats1, feats2, 'lambda', lambda);
    
    all_map.plot_transferred_xyz();    
    
    err_map        = Functional_Map(LB1, LB2);
    err_map.compute_f_map(fmap_method, neigs, neigs, rnd_feats_1, rnd_feats_2, 'lambda', lambda);
        
%% NEW    
    tic
    fmap_method    = 'frobenius_with_covariance';
    lambda         = 20;                                                    
    cov_map        = Functional_Map(LB1, LB2);    
    cov_map.compute_f_map(fmap_method, neigs, neigs, feats1, feats2, 'lambda', lambda);  
    toc
        
%% Looking on mis-alignment error.
%  Observe that the optimized map on noise, has smaller misalignment on the trained noise than a orthodox map.
    Fs = all_map.projected_source_features();
    Ft = all_map.projected_target_features();
    misalignment = sum(sum(abs((all_map.fmap * Fs) - Ft), 1)) ./ size(Fs, 2) % Average alignment error (L1 dist) of probe functions.       
    
    Fs = err_map.projected_source_features();
    Ft = err_map.projected_target_features();
    
    misalignment = sum(sum(abs((all_map.fmap * Fs) - Ft), 1)) ./ size(Fs, 2) % Average alignment error (L1 dist) of probe functions.          
    misalignment = sum(sum(abs((err_map.fmap * Fs) - Ft), 1)) ./ size(Fs, 2) % Average alignment error (L1 dist) of probe functions.
    
    err_map.plot_transferred_xyz();
    
%%  Evaluate F-maps.
    fid = fopen('../data/input/tosca_symmetries/michael.sym'); % TODO-P add to IO.read_symmetries(); 
    C   = textscan(fid, '%s', 'delimiter', '\n');   % Read symmetries
    fclose(fid);
    symmetries  = str2double(C{:});    
    groundtruth = (1:mesh1.num_vertices)';          % Groundtruth node-correspondence (Tosca: within same class, i node matches i).

    [dists_a, indices] = all_map.pairwise_distortion(groundtruth, 'nsamples', 500,     'symmetries', symmetries);                
    [dists_c, ~]       = cov_map.pairwise_distortion(groundtruth, 'indices',  indices, 'symmetries', symmetries);
%     [dists_e,  ~]      = err_map.pairwise_distortion(groundtruth, 'indices',  indices, 'symmetries', symmetries);
    
    mean(dists_a)
    mean(dists_c)
%     mean(dists_e)
    
%% Ensembling.    
    num_of_maps = 30;
    lambda      = 20;                                        % Regulizer meta-parameter.
    ensembles   = cell(num_of_maps, 1);
    all_feats   = size(feats1.F, 2);
    min_feats   = ceil(0.1 * all_feats);    
    fmap_method = 'frobenius';
    
    Fs = all_map.projected_source_features();
    Ft = all_map.projected_target_features();

    for i=1:num_of_maps        

        feats_i     = floor(min_feats + rand() * (all_feats - min_feats));        
        sample_feat = randsample(all_feats, feats_i)';

        source_f = feats1.keep_only(sample_feat);
        target_f = feats2.keep_only(sample_feat);
        
%         source_f.F  = feats1.F(:, sample_feat);                             % TODO-P remove .F and make it a class.
%         target_f.F  = feats2.F(:, sample_feat);       
                
        ensembles{i}.M = Functional_Map(LB1, LB2);        
        ensembles{i}.M.compute_f_map(fmap_method, neigs, neigs, source_f, target_f, 'lambda', lambda);
        
        fmap = ensembles{i}.M;
                
%         Fs = fmap.projected_source_features();
%         Ft = fmap.projected_target_features();
%         error_i_l1    = sum(sum(abs((fmap.fmap * Fs) - Ft), 1)) ./ feats_i;               % Average alignment error (L1 dist) of probe functions.

        error_i_l1    = sum(sum(abs((fmap.fmap * Fs) - Ft), 1)) ./ all_feats;               % Average alignment error (L1 dist) of all functions.
        fprintf('%f, %d\n', error_i_l1, feats_i);
        
%         source_f.index 
        
        ensembles{i}.error       = error_i_l1;
    end
%%    
% Aggregating Ensembles
%   1. Just the average of them (i.e., bagging).
    Xmean = zeros(size(ensembles{1}.M.fmap));
    for i=1:num_of_maps            
        Xmean = Xmean + ensembles{i}.M.fmap;
    end
    Xmean = (1/(num_of_maps)) .* Xmean;
    Fmean = Functional_Map(LB1, LB2);
    Fmean.set_fmap(Xmean);

% %   2. Linear combination of fmaps based on individual error weights.    
%     Xlin = zeros(size(ensembles{1}.M.fmap));
%     for i=1:num_of_maps            
%         Xlin = Xlin + ((1/ensembles{i}.error) * ensembles{i}.M.fmap);
%     end
%     Flin = Functional_Map(LB1, LB2);
%     Flin.set_fmap(Xlin ./ num_of_maps);    
    
    
%   3. Combine the learned errors to create a map that emphasizes them.        
    all_errors  = zeros(num_of_maps, 1);
    for i=1:num_of_maps
        all_errors(i) = ensembles{i}.error;          
    end
    total_error = sum(all_errors);    
    w = zeros(num_of_maps, 1);
    
    for i=1:num_of_maps            
        w(i) = max(all_errors) - all_errors(i) + min(all_errors);
    end    
    w = w ./ sum(w);
    
    Xlin2 = zeros(size(ensembles{1}.M.fmap));
    for i=1:num_of_maps 
        Xlin2 = Xlin2  + (w(i) .* ensembles{i}.M.fmap);
    end

    Flin2 = Functional_Map(LB1, LB2);
    Flin2.set_fmap(Xlin2);
    
%   2. TEST-delete
    [pos, m]= max(w)      ;
    Flin = ensembles{2}.M ;
    
    
    
% Evaluating distortions.
    % Load symmetries.    
    fid = fopen('../data/input/tosca_symmetries/michael.sym'); % TODO-P add to IO.read_symmetries();
    C   = textscan(fid, '%s', 'delimiter', '\n');          % Read symmetries
    fclose(fid);
    symmetries   = str2double(C{:});    
    
    % Make groundtuth correspondence (for tosca, within same class i node matched i).
    groundtruth = (1:mesh1.num_vertices)';      
      
    
    %%
    [dists_m, indices]       = Fmean.pairwise_distortion(groundtruth, 'nsamples', 400, 'symmetries', symmetries);            
    
    %%
    
    [dists_l, ~]             = Flin.pairwise_distortion(groundtruth,  'indices', indices, 'symmetries', symmetries);                    
    [dists_l2, ~]            = Flin2.pairwise_distortion(groundtruth, 'indices', indices, 'symmetries', symmetries);        
    
    % Use All probe functions.
    
    fmap.compute_f_map(fmap_method, neigs, neigs, feats1, feats2, 'lambda', lambda);                      
    [dists_a, ~] = fmap.pairwise_distortion(groundtruth, 'indices', indices, 'symmetries', symmetries);
            
    % Use the optimal/groundtruth map.
    X_opt               = Functional_Map.groundtruth_functional_map(LB1.evecs(neigs), LB2.evecs(neigs), groundtruth, diag(LB2.A));       
    [dists_o,  ~]       = Functional_Map.pairwise_distortion_of_map(X_opt, LB1, LB2, groundtruth, 'indices', indices, 'symmetries', symmetries);
    
    % Use only the wks/hks features.    
    feats1_intr.F = feats1.F(:, 1: wks_samples + hks_samples);                                                      % TODO-P utilize feats.index.
    feats2_intr.F = feats2.F(:, 1: wks_samples + hks_samples);
    fmap.compute_f_map(fmap_method, neigs, neigs, feats1_intr, feats2_intr, 'lambda', lambda);               % Use hks/wks only.
    [dists_i, ~] = fmap.pairwise_distortion(groundtruth, 'indices', indices, 'symmetries', symmetries);
    
    mean(dists_o)  % optimal dist
    mean(dists_l)  % ensemble linear combination 1
    mean(dists_l2) % ensemble linear combination 2
    mean(dists_i)  % all wks/hks (intrinsic)
    mean(dists_a)  % all probes
    mean(dists_m)  % bagging/mean ensemble

   
    
