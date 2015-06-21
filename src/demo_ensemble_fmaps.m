%%  A Script researching Ensemble Procedures for Functional Maps.
    clr;
    gitdir;
    cd FmapLib/src

%% Load a Mesh and calculate basic quantities.
    meshfile  = '../data/kid_rodola/0001.isometry.1.off';
    inmesh    = Mesh(meshfile, 'rodola_1_1');    
    inmesh.set_triangle_angles();
    inmesh.set_vertex_areas('barycentric');            
    sum(inmesh.get_vertex_areas('barycentric'))             % We don't normalize vertex areas to sum to 1.
        
    num_eigs       = 200;
    LB             = Laplace_Beltrami(inmesh);              % Calculate the first 100 spectra, based on barycentric vertex areas.
    [evals, evecs] = LB.get_spectra(num_eigs, 'barycentric'); 
    save('../data/output/mesh_and_LB', 'inmesh', 'LB');
 
    % Load Precomputed ones.
%     load('../data/output/mesh_and_LB', 'inmesh', 'LB');
%     [evals, evecs] = LB.get_spectra(num_eigs, 'barycentric');
%     

%%  Compute Mesh Features.
    wks_samples = 100; hks_samples = 100; curvatures = 100;
    [energies, sigma] = Mesh_Features.energy_sample_generator('log_linear', evals(2), evals(end), wks_samples);
    wks_sig           = Mesh_Features.wave_kernel_signature(evecs(:,2:end), evals(2:end), energies, sigma);    
    heat_time         = Mesh_Features.energy_sample_generator('log_sampled', evals(2), evals(end), hks_samples);
    hks_sig           = Mesh_Features.heat_kernel_signature(evecs(:,2:end), evals(2:end), heat_time);    
    heat_time         = Mesh_Features.energy_sample_generator('log_sampled', evals(2), evals(end), curvatures-1);
    mean_curvature    = Mesh_Features.mean_curvature(inmesh, LB, heat_time);    
    gauss_curvature   = Mesh_Features.gaussian_curvature(inmesh, LB, heat_time);
    
    %     TODO-P,E Normalize prob functions (before or after?)
    source_probes     = LB.project_functions('barycentric', num_eigs, wks_sig, hks_sig, mean_curvature, gauss_curvature);

    save('../data/output/features_1', 'source_probes')
%%
    num_feat = size(source_probes, 2);
    