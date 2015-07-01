%%  A Script demonstrating the basic functionalities of the FmapLib (Work in progress).
    clr;
    gitdir;
    cd FmapLib/src

%% Load a Mesh and calculate basic quantities.
    meshfile  = '../data/input/tosca_cats/cat0.off';
    inmesh    = Mesh(meshfile, 'cat0');    
    inmesh.set_triangle_angles();
    inmesh.set_vertex_areas('barycentric');              
    LB = Laplace_Beltrami(inmesh, inmesh.get_vertex_areas('barycentric'));
    eigs_num = 100;
    [evals, evecs] = LB.get_spectra(eigs_num);  % Calculate the first 100 spectra, based on barycentric vertex areas.
    save('../data/output/cat0_mesh_and_100_LB', 'inmesh', 'LB'); 
% Load Precomputed ones.
%     load('../data/output/cat0_mesh_and_100_LB', 'inmesh', 'LB');
%     [evals, evecs] = LB.get_spectra(eigs_num);

%%
    tic
    pairs = 1000; num_bins = 100;
    D21 = Mesh_Features.D2_shape_distribution(inmesh, pairs, num_bins);    
    toc

    tic
    pairs = 10000; num_bins = 100;
    D22 = Mesh_Features.D2_shape_distribution(inmesh, pairs, num_bins);
    toc
    
    tic
    pairs = 30000; num_bins = 100;
    D23 = Mesh_Features.D2_shape_distribution(inmesh, pairs, num_bins);
    toc

%% Two Meshes and an F-map.
    num_eigs       = 150;
    wks_samples    = 150;
    hks_samples    = 100;
  
    meshfile       = '../data/kid_rodola/0001.isometry.1.off';
    mesh1          = Mesh(meshfile, 'rodola_1_1');        
    mesh1.set_vertex_areas('barycentric');    
%     LB1            = Laplace_Beltrami(mesh1, mesh1.get_vertex_areas('barycentric'));
%     [evals, evecs] = LB1.get_spectra(num_eigs);
%     save('../data/output/LB1', 'LB1');              
    load('../data/output/LB1');    
    [evals, evecs] = LB1.get_spectra(num_eigs);
    
    [energies, sigma] = Mesh_Features.energy_sample_generator('log_linear', evals(2), evals(end), wks_samples);
    wks_sig           = Mesh_Features.wave_kernel_signature(evecs(:,2:end), evals(2:end), energies, sigma);    
    heat_time         = Mesh_Features.energy_sample_generator('log_sampled', evals(2), evals(end), hks_samples);
    hks_sig           = Mesh_Features.heat_kernel_signature(evecs(:,2:end), evals(2:end), heat_time);
    
%     heat_time         = Mesh_Features.energy_sample_generator('log_sampled', evals(2), evals(end), curvatures-1);
%     mean_curvature    = Mesh_Features.mean_curvature(mesh1, LB1, heat_time);    
%     gauss_curvature   = Mesh_Features.gaussian_curvature(mesh1, heat_time);
        
%     source_probes     = LB1.project_functions(num_eigs, wks_sig, hks_sig);     
    source_probes     = [hks_sig wks_sig];
    
    meshfile          = '../data/kid_rodola/0002.isometry.1.off';
    mesh2             = Mesh(meshfile, 'rodola_2_1');    
    mesh2.set_default_vertex_areas('barycentric');    
%     LB2               = Laplace_Beltrami(mesh2); 
%     [evals, evecs]    = LB2.get_spectra(num_eigs);
%     save('../data/output/LB2', 'LB2');              
    load('../data/output/LB2');    
    [evals, evecs] = LB2.get_spectra(num_eigs);
    
    [energies, sigma] = Mesh_Features.energy_sample_generator('log_linear', evals(2), evals(end), wks_samples);
    wks_sig           = Mesh_Features.wave_kernel_signature(evecs(:,2:end), evals(2:end), energies, sigma);    
    heat_time         = Mesh_Features.energy_sample_generator('log_sampled', evals(2), evals(end), hks_samples);
    hks_sig           = Mesh_Features.heat_kernel_signature(evecs(:,2:end), evals(2:end), heat_time);
    
%     target_probes     = LB2.project_functions(num_eigs, wks_sig, hks_sig);         
    target_probes     = [hks_sig wks_sig];

    % Load symmetries.    
    C = textread('../data/kid_rodola/sym.txt', '%s', 'delimiter', ' ');  % Read symmetries:
    C = cell2mat(C); symmetries = str2num(C);            
    

   
    %%    
    gt_map = (1:mesh1.num_vertices)';                   % Ground truth correspondences from Source_Mesh to Target_Mesh.        
    F = Functional_Map(LB1, LB2);
    X1 = F.compute_f_map('functions_only', 40, 40, source_probes, target_probes, 'normalize', 1);    
    [dists, random_points] = F.pairwise_distortion(gt_map, 'nsamples', 200, 'symmetries', symmetries);
    mean(dists)    
    [dists2, random_points2] = Functional_Map.pairwise_distortion_of_map(F.fmap, LB1, LB2, gt_map, 'indices', random_points, 'symmetries', symmetries);       
    assert(all(dists==dists2) && all(random_points == random_points2));   
    
    
    %%
    [~, source_basis] = LB1.get_spectra(num_eigs);
    [~, target_basis] = LB2.get_spectra(num_eigs);    
    X_opt             = Functional_Map.groundtruth_functional_map(source_basis, target_basis, gt_map, mesh2.get_vertex_areas('barycentric'));  %TODO-P, Ugly barycentric.    
    [dists3, random_points3] = Functional_Map.pairwise_distortion_of_map(X_opt, LB1, LB2, gt_map, 'indices', random_points, 'symmetries', symmetries); 
    all(random_points == random_points3)
    mean(dists3)
    %%

    
    
    %% Make the Fmap Using all eigenvectors    
%     lambda            = 20;                
%     X                 = Functional_Map.sum_of_squared_frobenius_norms(source_probes, target_probes, LB1.evals(num_eigs) , LB2.evals(num_eigs), lambda);                

    %%  Add to unit-test    
    % Make the groundtruth Functional Map
    gt_map            = (1:mesh1.num_vertices)';   % Ground truth correspondences from Source_Mesh to Target_Mesh.    
    [~, source_basis] = LB1.get_spectra(num_eigs);
    [~, target_basis] = LB2.get_spectra(num_eigs);
    X_opt             = Functional_Map.groundtruth_functional_map(source_basis, target_basis, gt_map, mesh2.get_vertex_areas('barycentric'));  %TODO-P, Ugly.
    % Evaluate the X_opt
    eval_points = 200;
%     [dists,  random_points]  = Functional_Map.pairwise_distortion_of_map(X_opt, mesh1, mesh2, source_basis, target_basis, gt_map, 'nsamples', eval_points);            
    [dists2, random_points2] = Functional_Map.pairwise_distortion_of_map(X_opt, LB1, LB2, gt_map, 'indices', random_points);
    assert(length(dists) == eval_points);
    assert(all(dists==dists2) && all(random_points == random_points2));
    length(unique(random_points)) == length(random_points)
    %%
    % Use symmetries.    
    C = textread('../data/kid_rodola/sym.txt', '%s', 'delimiter', ' ');  % Read symmetries:
    C = cell2mat(C); sym = str2num(C);            
    [dists3, random_points3] = Functional_Map.pairwise_distortion_of_map(X_opt, mesh1, mesh2, source_basis, target_basis, gt_map, 'indices', random_points, 'symmetries', sym);
    assert(all(dists3 <= dists2));

    %% Use a small number of eigenvectors to do the Fmap.
    [~, source_basis] = LB1.get_spectra(10, 'barycentric');
    [~, target_basis] = LB2.get_spectra(10, 'barycentric');    
    X_opt_small       = Functional_Map.groundtruth_functional_map(source_basis, target_basis, gt_map, mesh2.get_vertex_areas('barycentric'));
    [dists2, random_points2] = Functional_Map.pairwise_distortion_of_map(X_opt_small, mesh1, mesh2, source_basis, target_basis, gt_map, 'indices', random_points, 'fast');
    assert(mean(dists2) > mean(dists)) % TODO-P Change to something more reasonable.
    
    
    

    
    