%%  A Script researching Ensemble Procedures for Functional Maps.
    clr;
    gitdir;
    cd FmapLib/src

%% Load two Meshes and make a F-map.
    num_eigs          = 300;
    wks_samples       = 150;
    hks_samples       = 100;
    curvature_samples = 200;
  
    meshfile       = '../data/input/kid_rodola/0001.isometry.1.off';
    mesh1          = Mesh(meshfile, 'rodola_1_1');        
    mesh1.set_default_vertex_areas('barycentric');    
%     LB1            = Laplace_Beltrami(mesh1);
%     [evals, evecs] = LB1.get_spectra(num_eigs);
%     save('../data/output/ensembles/LB1', 'LB1');              
    load('../data/output/ensembles/LB1');    
    [evals, evecs] = LB1.get_spectra(num_eigs);
    
    [energies, sigma] = Mesh_Features.energy_sample_generator('log_linear', evals(2), evals(end), wks_samples);
    wks_sig           = Mesh_Features.wave_kernel_signature(evecs(:,2:end), evals(2:end), energies, sigma);    
    heat_time         = Mesh_Features.energy_sample_generator('log_sampled', evals(2), evals(end), hks_samples);
    hks_sig           = Mesh_Features.heat_kernel_signature(evecs(:,2:end), evals(2:end), heat_time);
    heat_time         = Mesh_Features.energy_sample_generator('log_sampled', evals(2), evals(end), curvature_samples-1);
    mc_sig            = Mesh_Features.mean_curvature(mesh1, LB1, heat_time);    
    gc_sig            = Mesh_Features.gaussian_curvature(mesh1, LB1, heat_time);
    source_probes     = [hks_sig wks_sig mc_sig gc_sig];    

        
    meshfile          = '../data/input/kid_rodola/0002.isometry.1.off';
    mesh2             = Mesh(meshfile, 'rodola_2_1');    
    mesh2.set_default_vertex_areas('barycentric');    
%     LB2               = Laplace_Beltrami(mesh2); 
%     [evals, evecs]    = LB2.get_spectra(num_eigs);
%     save('../data/output/ensembles/LB2', 'LB2');                  
    load('../data/output/ensembles/LB2');    
    [evals, evecs] = LB2.get_spectra(num_eigs);
    
    [energies, sigma] = Mesh_Features.energy_sample_generator('log_linear', evals(2), evals(end), wks_samples);
    wks_sig           = Mesh_Features.wave_kernel_signature(evecs(:,2:end), evals(2:end), energies, sigma);    
    heat_time         = Mesh_Features.energy_sample_generator('log_sampled', evals(2), evals(end), hks_samples);
    hks_sig           = Mesh_Features.heat_kernel_signature(evecs(:,2:end), evals(2:end), heat_time);
    heat_time         = Mesh_Features.energy_sample_generator('log_sampled', evals(2), evals(end), curvature_samples-1);
    mc_sig            = Mesh_Features.mean_curvature(mesh2, LB2, heat_time);    
    gc_sig            = Mesh_Features.gaussian_curvature(mesh2, LB2, heat_time);
    target_probes     = [hks_sig wks_sig mc_sig gc_sig];    

%% Ensembling
    num_of_maps = 10;
    ensembles   = cell(1, num_of_maps);
    all_feats   = size(source_probes, 2);
    min_feats   = ceil(0.1*all_feats);
    for i=1:num_of_maps        
        feats_i      = floor(min_feats + rand() * (all_feats - min_feats));
        sample_feat  = randsample(all_feats, feats_i)';
        source_f     = source_probes(:, sample_feat)  ;
        target_f     = target_probes(:, sample_feat)  ;
        Fi           = Functional_Map(LB1, LB2);
        Fi.compute_f_map('functions_only', num_eigs, num_eigs, source_f, target_f, 'normalize', 1);    
        
        ensembles{i}.map = Fi;               
        ensembles{i}.nfeat = feats_i;
        ensembles{i}.s_feats = source_f;
        ensembles{i}.t_feats = target_f;
    end
    
%% Aggregating Ensembles
    Xmean = zeros(size(ensembles{1}.map));
    for i=1:num_of_maps            
        Xmean = Xmean + ensembles{i}.map.fmap;
    end
    Xmean = (1/(num_of_maps)) .* Xmean;
    Fmean = Functional_Map(LB1, LB2);
    Fmean.set_fmap(Xmean);
        
%% Evaluating distortions
    % Load symmetries.    
    C = textread('../data/kid_rodola/sym.txt', '%s', 'delimiter', ' ');  % Read symmetries:
    C = cell2mat(C); symmetries = str2num(C);           
    % Make groundtuth     
    groundtruth = (1:mesh1.num_vertices)';                  

    all_dists = cell(num_of_maps, 1);    
    indices   = randsample(mesh1.num_vertices, 200);
    
    for i=1:num_of_maps        
        [dists, ~]       = ensembles{i}.map.pairwise_distortion(groundtruth, 'indices', indices, 'symmetries', symmetries);        
        [dists2, ~]      = Fmean.pairwise_distortion(groundtruth, 'indices', indices, 'symmetries', symmetries);        
        all_dists{i}.ens = dists;
        all_dists{i}.agg = dists2;
    
    end        
    for i=1:num_of_maps       %TODO-P save sampling numbers
        fprintf('Ensmble with %d feautures. Mean: %f STD:%f \n', ensembles{i}.nfeat, mean(all_dists{i}.ens), std(all_dists{i}.ens))
        fprintf('A_Fmap Mean:%f STD:%f \n', mean(all_dists{i}.agg), std(all_dists{i}.agg))
        ensembles{i}.nfeat
    end

   
    
