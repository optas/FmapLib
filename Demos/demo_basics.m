%%  A Script demonstrating the basic functionalities of the FmapLib on Shapes (Work in progress).

%  Load mesh 1.
meshfile       = '../data/input/tosca_small/michael1.off';
mesh1          = Mesh(meshfile, 'mike1');
mesh1.set_default_vertex_areas('barycentric');              % Associate an area with each vertex via the 'barycentric' rule.

%  Load mesh 2.
meshfile       = '../data/input/tosca_small/michael2.off';
mesh2          = Mesh(meshfile, 'mike2');  
mesh2.set_default_vertex_areas('barycentric');

% Plot meshes.
mesh1.plot();
mesh2.plot();


% Compute a basis for functions defined on the mesh vertices (here, Laplace Beltrami).
LB1            = Laplace_Beltrami(mesh1);                   % Uses the cotangent scheme for the laplacian discretisation.
LB2            = Laplace_Beltrami(mesh2);

% Generate functions over the mesh vertices.
feats1         = Mesh_Features(mesh1, LB1);                 % Mesh node features.
feats2         = Mesh_Features(mesh2, LB2);

% Parameters for the function generation.
hks_samples    = 100;                                        % Feature dimensions.
wks_samples    = 100; 
mc_samples     = 100;  
gc_samples     = 100;
neigs          = 50;                                        % LB eigenvecs to be used.

feats1.compute_default_feautures(neigs, wks_samples, hks_samples, mc_samples, gc_samples);
feats2.compute_default_feautures(neigs, wks_samples, hks_samples, mc_samples, gc_samples);

%% Plot some feature functions.
close all;
feats1.index
mesh1.plot(feats1.F(:,10));
mesh1.plot(feats1.F(:,210));

%% Make a functional map and plot some derived measures.
fmap_method   = 'frobenius_square';
lambda        = 20;
map1          = Functional_Map(LB1, LB2);
map1.compute_f_map(fmap_method, neigs, neigs, feats1, feats2, 'lambda', lambda);

map1.plot_transferred_xyz();
map1.plot_area_distortion();
map1.plot_conformal_distortion();

%% Further Evaluate maps with geodesic criteria.
symmetries   = Mesh_IO.read_symmetries('../data/input/tosca_symmetries/michael.sym');
groundtruth = (1:mesh1.num_vertices)';          % Groundtruth node to node correspondences.
[dists_1, indices] = map1.pairwise_distortion(groundtruth, 'nsamples', 500, 'symmetries', symmetries);                    
mean(dists_1)    



