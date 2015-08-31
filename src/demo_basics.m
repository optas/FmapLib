%%  A Script demonstrating the basic functionalities of the FmapLib (Work in progress).
% clr;
% gitdir;
% cd FmapLib/src
close all; clear all; clc;
%% Mini exposition of the FmapLib.
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

% Compute a basis for functions defined on the mesh vertices (Laplace Beltrami).
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
neigs          = 100;                                        % LB eigenvecs to be used.

feats1.compute_default_feautures(neigs, wks_samples, hks_samples, mc_samples, gc_samples);
feats2.compute_default_feautures(neigs, wks_samples, hks_samples, mc_samples, gc_samples);
%%
close all;
mesh1.plot(feats1.F(:,301))                     % TODO-E, seems like first derived gc probe functions are really bad. why? fix?
figure; plot(feats1.F(:,301))                   % If not a bug/or easy fix, can we plot something more reallistic than the 'constant' function?
mesh1.plot(feats1.F(:,302))
figure; plot(feats1.F(:,302))                   % etc.

%% Make some functional maps.
% 1.
fmap_method   = 'frobenius';
lambda        = 20;
map1          = Functional_Map(LB1, LB2);
map1.compute_f_map(fmap_method, neigs, neigs, feats1, feats2, 'lambda', lambda);
map1.plot_transferred_xyz();
map1.plot_area_distortion();
map1.plot_conformal_distortion();

% 2.
gt_map = Functional_Map.groundtruth_functional_map(LB1, LB2, (1:mesh1.num_vertices)', neigs, neigs);
gt_map.plot_transferred_xyz();
gt_map.plot_area_distortion();
gt_map.plot_conformal_distortion();


%% Further Evaluate maps.   
fid = fopen('../data/input/tosca_symmetries/michael.sym'); % TODO-P add to IO.read_symmetries(); 
C   = textscan(fid, '%s', 'delimiter', '\n');   % Read symmetries
fclose(fid);
symmetries  = str2double(C{:});    
groundtruth = (1:mesh1.num_vertices)';          % Groundtruth node to node correspondences.
[dists_a, indices] = map1.pairwise_distortion(groundtruth, 'nsamples', 500, 'symmetries', symmetries);                    
mean(dists_a)    