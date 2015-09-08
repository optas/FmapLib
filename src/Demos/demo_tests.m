%%
%%    Script for dirty tests (to be deleted).
%%
 
[dp, sp] = get_project_paths('FmapLib');

%  Load mesh 1.
meshfile       = [dp '/input/tosca_small/michael1.off'];
mesh1          = Mesh(meshfile, 'mike1');
mesh1.set_default_vertex_areas('barycentric');              % Associate an area with each vertex via the 'barycentric' rule.

% Compute a basis for functions defined on the mesh vertices (Laplace Beltrami).
LB1            = Laplace_Beltrami(mesh1);                   % Uses the cotangent scheme for the laplacian discretisation.

% Generate functions over the mesh vertices.
feats1         = Mesh_Features(mesh1, LB1);                 % Mesh node features.

% Parameters for the function generation.
hks_samples    = 0;                                        % Feature dimensions.
wks_samples    = 0; 
mc_samples     = 2; 
gc_samples     = 0;
neigs          = 50;                                        % LB eigenvecs to be used.

feats1.compute_default_feautures(neigs, wks_samples, hks_samples, mc_samples, gc_samples);

%% Preparing mesh saliency
top_k = 30;
[ids, dists]     = knnsearch(mesh1.vertices, mesh1.vertices, 'K', top_k + 1);
neighborhoods.distances = dists;
neighborhoods.ids       = ids; 

vertex_function = feats1.F(:,2);

sigma = 2 * 0.003  * mesh1.max_euclidean_distance();
cut_off = 2 * sigma;

F1 = Mesh_Features.gaussian_weighted_average(vertex_function, neighborhoods, sigma, cut_off);

sigma = 2 * sigma ; 
cut_off = 2 * sigma;
F2 = Mesh_Features.gaussian_weighted_average(vertex_function, neighborhoods, sigma, cut_off);
%%

F = abs(F1 - F2);
mesh1.plot(vertex_function)
mesh1.plot(F);
mean(F-vertex_function)

















%% Load two Meshes and compute the F-map.
    num_eigs       = 15;
    wks_samples    = 15;
    hks_samples    = 10;
  
    meshfile       = '../data/kid_rodola/0001.isometry.1.off';
    mesh1          = Mesh(meshfile, 'rodola_1_1');        
    mesh1.set_vertex_areas('barycentric');    
    LB1            = Laplace_Beltrami(mesh1, mesh1.get_vertex_areas('barycentric'));
    [evals, evecs] = LB1.get_spectra(num_eigs);
%     save('../data/output/LB1', 'LB1');              
%     load('../data/output/LB1');    
    [evals, evecs] = LB1.get_spectra(num_eigs);
    
    [energies, sigma] = Mesh_Features.energy_sample_generator('log_linear', evals(2), evals(end), wks_samples);
    wks_sig           = Mesh_Features.wave_kernel_signature(evecs(:,2:end), evals(2:end), energies, sigma);    
    heat_time         = Mesh_Features.energy_sample_generator('log_sampled', evals(2), evals(end), hks_samples);
    hks_sig           = Mesh_Features.heat_kernel_signature(evecs(:,2:end), evals(2:end), heat_time);
            
%     source_probes     = LB1.project_functions(num_eigs, wks_sig, hks_sig);     
    source_probes     = [hks_sig wks_sig];
    
    meshfile          = '../data/kid_rodola/0002.isometry.1.off';
    mesh2             = Mesh(meshfile, 'rodola_2_1');    
    mesh2.set_default_vertex_areas('barycentric');    
    LB2               = Laplace_Beltrami(mesh2); 
    [evals, evecs]    = LB2.get_spectra(num_eigs);
%     save('../data/output/LB2', 'LB2');              
%     load('../data/output/LB2');    
    [evals, evecs] = LB2.get_spectra(num_eigs);
    
    [energies, sigma] = Mesh_Features.energy_sample_generator('log_linear', evals(2), evals(end), wks_samples);
    wks_sig           = Mesh_Features.wave_kernel_signature(evecs(:,2:end), evals(2:end), energies, sigma);    
    heat_time         = Mesh_Features.energy_sample_generator('log_sampled', evals(2), evals(end), hks_samples);
    hks_sig           = Mesh_Features.heat_kernel_signature(evecs(:,2:end), evals(2:end), heat_time);
    
%     target_probes     = LB2.project_functions(num_eigs, wks_sig, hks_sig);         
    target_probes     = [hks_sig wks_sig];

    F = Functional_Map(LB1, LB2);    
    F.compute_f_map('functions_only', 40, 40, source_probes, target_probes, 'normalize', 1);    
    gt_map = (1:mesh1.num_vertices)';   % Ground truth correspondences from Source_Mesh to Target_Mesh.
    
    C = textread('../data/kid_rodola/sym.txt', '%s', 'delimiter', ' ');  % Read symmetries:
    C = cell2mat(C); symmetries = str2num(C);            

    
%% Test 1: pseudo-inverse VS. Our theretical one.        
    Proj1 = LB1.evecs(10)' * LB1.A;
    Proj2 = pinv(LB1.evecs(10));    
    imagesc(Proj1);
    figure; imagesc(Proj2);
    all_close(Proj1, Proj2, 1e-6, +Inf)
    all_close(Proj1, Proj2, +Inf, 10000)
    all_close(Proj1, Proj2, +Inf, 100000)
    % Panos conclusion: our projector is scaled poorly (seems to be more than 10000 times smaller).
 
%% 2. Area-Conformal map
%  Similarly I think when computing e.g., the area difference map with your way that is independent of the basis we get
%  results different that when we exploit that the underlying one is the laplace beltrami. Can you see if there is
%  a bug or we understand how/where we might be different?

    D1 = F.area_difference();
    D2 = F.area_difference2();
    imagesc(D1); figure; imagesc(D2);
    all_close(D1, D2, 1e-2, +Inf)

% 2b.  Can you please add code for computing conformal maps (independent of basis).
%     DC    = F.conformal_difference();
    
%% 3. On groundtruth map: 
%  See Functional_Map.groundtruth_functional_map
%  we do not use the areas of the vertices of the from mesh which I think is weird. Is it?


% 4. Shall we normalize the projected probe functions to have unit euclidean norm, or unit norm wrt. Areas?        
% See the normalize argument on the compute fmap.
    X1 = F.compute_f_map('functions_only', 40, 40, source_probes, target_probes, 'normalize', 1);    
    [dists, random_points] = F.pairwise_distortion(gt_map, 'nsamples', 200, 'symmetries', symmetries);
    mean(dists)    
   
    X2 = F.compute_f_map('functions_only', 40, 40, source_probes, target_probes, 'area_normalize', 1);        
    [dists2, random_points2] = F.pairwise_distortion(gt_map, 'indices', random_points, 'symmetries', symmetries);           
    all(random_points == random_points2)
    mean(dists2)

%    Is mean(dists) or mean(dists2) better?
    
    

