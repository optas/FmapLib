%%  A Script demonstrating the basic functionalities of the FmapLib (Work in progress).
    clr;
    gitdir;
    cd FmapLib/src

%% Load a Mesh and calculate basic quantities.
    meshfile  = '../data/kid_rodola/0001.isometry.1.off';
    inmesh    = Mesh(meshfile, 'rodola_1_1');
    inmesh.set_triangle_angles();
    inmesh.set_vertex_areas('barycentric');
    
    % Calculate the first 100 spectra, based on barycentric vertex areas.
    LB             = Laplace_Beltrami(inmesh);
    [evals, evecs] = LB.get_spectra(100, 'barycentric');
    save('../data/output/mesh_and_LB', 'inmesh', 'LB');

    % Load Precomputed ones.
%     load('../data/output/mesh_and_LB', 'inmesh', 'LB');
%     [evals, evecs] = LB.get_spectra(100, 'barycentric');


%%  %%  Testing new way of computing geodesics. 
    id                = 1;
    set_indicator     = zeros(inmesh.num_vertices, 1);
    set_indicator(id) = 1;
       
   
    [geo_dist]        = Mesh_Features.geodesic_distance_to_set(inmesh, LB, set_indicator);
        
    figure;
    trisurf(inmesh.triangles, inmesh.vertices(:,1), inmesh.vertices(:,2), inmesh.vertices(:,3), geo_dist);
    axis equal; shading interp;
      
    geo_dist2 = comp_geodesics_to_all(inmesh.vertices(:,1), inmesh.vertices(:,2), inmesh.vertices(:,3), inmesh.triangles', id, 1);
       
    figure;
    trisurf(inmesh.triangles, inmesh.vertices(:,1), inmesh.vertices(:,2), inmesh.vertices(:,3), geo_dist2);
    axis equal; shading interp; colorbar;
    
    d1 = geo_dist ./ max(geo_dist);
    d2 = geo_dist2 ./ max(geo_dist2);    
    norm(d1-d2) / norm(d2)
    norm(geo_dist -geo_dist2) / norm(geo_dist)

       
    %% Two Meshes and a F-map.
    num_eigs       = 100;
    wks_samples    = 150;
    hks_samples    = 100;
    curvatures     = 100;

    meshfile       = '../data/kid_rodola/0001.isometry.1.off';
    mesh1          = Mesh(meshfile, 'rodola_1_1');        
    LB1            = Laplace_Beltrami(mesh1);    
    [evals, evecs] = LB1.get_spectra(num_eigs, 'barycentric');
    save('../data/output/LB1', 'LB1');          

%     load('../data/output/LB1');    
%     [evals, evecs] = LB1.get_spectra(num_eigs, 'barycentric');

    [energies, sigma] = Mesh_Features.energy_sample_generator('log_linear', evals(2), evals(end), wks_samples);
    wks_sig           = Mesh_Features.wave_kernel_signature(evecs(:,2:end), evals(2:end), energies, sigma);    
    heat_time         = Mesh_Features.energy_sample_generator('log_sampled', evals(2), evals(end), hks_samples);
    hks_sig           = Mesh_Features.heat_kernel_signature(evecs(:,2:end), evals(2:end), heat_time);
    
%     heat_time         = Mesh_Features.energy_sample_generator('log_sampled', evals(2), evals(end), curvatures-1);
%     mean_curvature    = Mesh_Features.mean_curvature(mesh1, LB1, heat_time);    
%     gauss_curvature   = Mesh_Features.gaussian_curvature(mesh1, heat_time);
    
    %TODO-P Normalize prob functions
    from_probes       = LB1.project_functions('barycentric', num_eigs, wks_sig, hks_sig);
 
    meshfile          = '../data/kid_rodola/0002.isometry.1.off';
    mesh2             = Mesh(meshfile, 'rodola_2_1');
    LB2               = Laplace_Beltrami(mesh2);            
    [evals, evecs]    = LB2.get_spectra(num_eigs, 'barycentric');
    save('../data/output/LB2', 'LB2');          
    [energies, sigma] = Mesh_Features.energy_sample_generator('log_linear', evals(2), evals(end), wks_samples);
    wks_sig           = Mesh_Features.wave_kernel_signature(evecs(:,2:end), evals(2:end), energies, sigma);    
    heat_time         = Mesh_Features.energy_sample_generator('log_sampled', evals(2), evals(end), hks_samples);
    hks_sig           = Mesh_Features.heat_kernel_signature(evecs(:,2:end), evals(2:end), heat_time);
    to_probes         = LB2.project_functions('barycentric', num_eigs, wks_sig, hks_sig);    

    
%%  
%     lambda = 20;        
%     X      = Functional_Map.sum_of_squared_frobenius_norms(from_probes, to_probes, LB1.get_spectra(num_eigs, 'barycentric'), LB2.get_spectra(num_eigs, 'barycentric'), lambda); 
%%
    correspondences = [1:mesh1.num_vertices; 1:mesh2.num_vertices]';
    [~, basis_from] = LB1.get_spectra(num_eigs, 'barycentric');
    [~, basis_to]   = LB2.get_spectra(num_eigs, 'barycentric');
    X_opt = Functional_Map.groundtruth_functional_map(basis_from, basis_to, correspondences);

%% 

    groundtruth = (1:mesh1.num_vertices)';
    evaluation_samples  = 20;    
    [dist]              = Functional_Map.pairwise_distortion_of_map(X_opt, mesh1, mesh2, basis_from, basis_to, evaluation_samples, groundtruth, 1);                                          
    
