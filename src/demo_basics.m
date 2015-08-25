%%  A Script demonstrating the basic functionalities of the FmapLib (Work in progress).
    clr;
    gitdir;
    cd FmapLib/src

%% Mini exposition draft.
%  Input mesh 1.
meshfile       = '../data/input/tosca/cat0.off';
mesh_name      = 'cat0';
mesh1          = Mesh(meshfile, mesh_name);
mesh1.set_default_vertex_areas('barycentric');              % Associate an area with each vertex via the 'barycentric' rule.

LB1            = Laplace_Beltrami(mesh1);                   % Create a cotangent scheme mesh laplacian.
feats1         = Mesh_Features(mesh1, LB1);                 % Mesh node features.

hks_samples    = 30;                                        % Feature dimensions.
wks_samples    = 30; 
mc_samples     = 0; 
gc_samples     = 0;
neigs          = 32;                                        % LB eigenvecs to be used.

feats1.compute_default_feautures(neigs, wks_samples, hks_samples, mc_samples, gc_samples);

%  Input mesh 2.
meshfile       = '../data/input/tosca/cat1.off';
mesh2          = Mesh(meshfile, 'cat1');  
mesh2.set_default_vertex_areas('barycentric');

LB2            = Laplace_Beltrami(mesh2);
feats2         = Mesh_Features(mesh2, LB2);
feats2.compute_default_feautures(neigs, wks_samples, hks_samples, mc_samples, gc_samples);

% Make a functional map.
fmap = Functional_Map(LB1, LB2);
fmap.compute_f_map('frobenius_square', neigs, neigs, feats1, feats2, 'lambda', 0);
fmap.plot();

fmap.compute_f_map('frobenius_square', neigs, neigs, feats1, feats2, 'lambda', 20);
figure; fmap.plot();

fmap.compute_f_map('frobenius', neigs, neigs, feats1, feats2, 'lambda', 20);
figure; fmap.plot();
    
    
    
    
%% Load two Meshes and initialize their Feauture/LB classes.
    meshfile       = '../data/input/kid_rodola/0001.isometry.1.off';
    mesh1          = Mesh(meshfile, 'rodola_1_1');        
    mesh1.set_default_vertex_areas('barycentric');    
    LB1            = Laplace_Beltrami(mesh1);
    feats1         = Mesh_Features(mesh1, LB1);
    
    meshfile          = '../data/input/kid_rodola/0002.isometry.1.off';
    mesh2             = Mesh(meshfile, 'rodola_2_1');
    mesh2.set_default_vertex_areas('barycentric');    
    LB2               = Laplace_Beltrami(mesh2); 
    feats2            = Mesh_Features(mesh2, LB2);
            
%%  Compute Features           
    neigs          = 200;
    wks_samples    = 100; 
    hks_samples    = 100;    
    mc_samples     = 100; 
    gc_samples     = 100;
%%    
    feats1.compute_default_feautures(neigs, wks_samples, hks_samples, mc_samples, gc_samples);
    feats2.compute_default_feautures(neigs, wks_samples, hks_samples, mc_samples, gc_samples);    
%%        
%   save('../data/output/demo_basics', 'LB1', 'LB2', 'feats1', 'feats2', 'mesh1', 'mesh2');
    load('../data/output/demo_basics');    
       
%%  Make Fmaps
    F1 = Functional_Map(LB1, LB2);
    X1 = F1.compute_f_map('functions_only', neigs, neigs, feats1, feats2);        
    
    F2 = Functional_Map(LB1, LB2);
    X2 = F2.compute_f_map('frobenius_square', neigs, neigs, feats1, feats2, 'lambda', 20);                
    
%% Evaluate maps.   
%  Load symmetries, prepare groundtruth, call pairwise_distortions.
    C = textread('../data/input/kid_rodola/sym.txt', '%s', 'delimiter', ' ');  % Read symmetries:
    C = cell2mat(C); symmetries = str2num(C);               
    gt_map = (1:mesh1.num_vertices)';                   % Ground truth correspondences from Source_Mesh to Target_Mesh.        
%%    
    [dists1, random_points] = F1.pairwise_distortion(gt_map, 'nsamples', 200, 'symmetries', symmetries);     
    mean(dists1) 
    hist(dists1);    
    [dists2, ~]             = F2.pairwise_distortion(gt_map, 'indices', random_points, 'symmetries', symmetries);    
    mean(dists2) 
    figure; hist(dists2);
    
    
%%  Add to unit-test    
    gt_map            = (1:mesh1.num_vertices)';   % Ground truth correspondences from Source_Mesh to Target_Mesh.    
    [~, source_basis] = LB1.get_spectra(neigs);
    [~, target_basis] = LB2.get_spectra(neigs);  
    X_opt             = Functional_Map.groundtruth_functional_map(source_basis, target_basis, gt_map, mesh2.get_vertex_areas('barycentric'));  %TODO-P, Ugly.
    % Evaluate the X_opt
    eval_points = 200;
    [dists,  random_points]  = Functional_Map.pairwise_distortion_of_map(X_opt, LB1, LB2, gt_map, 'nsamples', eval_points, 'symmetries', symmetries);
    mean(dists)
    hist(dists)
    

%     [dists2, random_points2] = Functional_Map.pairwise_distortion_of_map(X_opt, LB1, LB2, gt_map, 'indices', random_points);
%     assert(length(dists) == eval_points);
%     assert(all(dists==dists2) && all(random_points == random_points2));
%     length(unique(random_points)) == length(random_points)
    %%
    % Use symmetries.    
    C = textread('../data/input/kid_rodola/sym.txt', '%s', 'delimiter', ' ');  % Read symmetries:
    C = cell2mat(C); sym = str2num(C);            
       
    %%
    [dists3, random_points3] = Functional_Map.pairwise_distortion_of_map(X_opt, mesh1, mesh2, source_basis, target_basis, gt_map, 'indices', random_points, 'symmetries', sym);
    assert(all(dists3 <= dists2));

    %% Use a small number of eigenvectors to do the Fmap.
    [~, source_basis] = LB1.get_spectra(10, 'barycentric');
    [~, target_basis] = LB2.get_spectra(10, 'barycentric');    
    X_opt_small       = Functional_Map.groundtruth_functional_map(source_basis, target_basis, gt_map, mesh2.get_vertex_areas('barycentric'));
    [dists2, random_points2] = Functional_Map.pairwise_distortion_of_map(X_opt_small, mesh1, mesh2, source_basis, target_basis, gt_map, 'indices', random_points, 'fast');
    assert(mean(dists2) > mean(dists)) % TODO-P Change to something more reasonable.
    
    
    

    
    