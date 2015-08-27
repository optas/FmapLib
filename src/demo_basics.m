%%  A Script demonstrating the basic functionalities of the FmapLib (Work in progress).
% clr;
% gitdir;
% cd FmapLib/src
close all; clear all; clc;
%% Mini exposition of the FmapLib.
%  Load mesh 1.

% meshfile       = '../data/input/tosca_small/gorilla1.off';
% mesh_name      = 'gorilla1';

meshfile       = '../data/input/tosca_small/victoria0.off';
mesh_name      = 'vica0';

mesh1          = Mesh(meshfile, mesh_name);
mesh1.set_default_vertex_areas('barycentric');              % Associate an area with each vertex via the 'barycentric' rule.

sum(mesh1.barycentric_v_area == 0)                          % 1. TODO-E: debug. Everything works as expected, however the connectivity of gorilla1 is bad there are 
                                                            % 1328 vertices unlinked to a triangle so the area at the vertex is
                                                            % zero. The problem come from the triangulation.
                                                            

%%  Load mesh 2.
meshfile       = '../data/input/tosca_small/michael2.off';
mesh1          = Mesh(meshfile, 'mike2');
mesh1.set_default_vertex_areas('barycentric'); 

meshfile       = '../data/input/tosca_small/michael1.off';
mesh2          = Mesh(meshfile, 'mike1');  
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
hks_samples    = 50;                                        % Feature dimensions.
wks_samples    = 50; 
mc_samples     = 20; 
gc_samples     = 20;
neigs          = 32;                                        % LB eigenvecs to be used.

feats1.compute_default_feautures(neigs, wks_samples, hks_samples, mc_samples, gc_samples);
feats2.compute_default_feautures(neigs, wks_samples, hks_samples, mc_samples, gc_samples);

%% Make some functional maps.
fmap_method   = 'frobenius';
lambda        = 20;
map1          = Functional_Map(LB1, LB2);
map1.compute_f_map(fmap_method, neigs, neigs, feats1, feats2, 'lambda', lambda);


%% Test on ground truth fmap
fmap_truth = Functional_Map.groundtruth_functional_map(LB1.evecs(neigs), LB2.evecs(neigs), (1:mesh1.num_vertices)', mesh2.get_vertex_areas());
map_truth  = Functional_Map(LB1, LB2);
map_truth.set_fmap(fmap_truth);

map_truth.plot();
map_truth.plot_transferred_xyz();                                % 2. TODO-E Is this a good way to show the point-point correspondence?
                                                                 % I change the color so it looks better. A part from that everything seems okay.
                    

%% Maximum distortion
% 3. TODO-E Add Ways to show the locations with the most distorted
% angle/area elements.                 
area_diff = map_truth.area_difference2();
[~,~,V] = svd(area_diff);
max_dist = LB1.evecs(neigs) * V(:,1);
mesh1.plot(max_dist);
title('Maximum area distortion');

conf_diff = map_truth.conformal_difference2();
[~,~,V] = svd(conf_diff);
max_dist = LB1.evecs(neigs) * V(:,1);
mesh1.plot(max_dist);
title('Maximum conformal distortion');

%% Further Evaluate maps.   
fid = fopen('../data/input/tosca_symmetries/michael.sym'); % TODO-P add to IO.read_symmetries(); 
C   = textscan(fid, '%s', 'delimiter', '\n');   % Read symmetries
fclose(fid);
symmetries  = str2double(C{:});    
groundtruth = (1:mesh1.num_vertices)';          % Groundtruth node to node correspondences.
[dists_a, indices] = map1.pairwise_distortion(groundtruth, 'nsamples', 500, 'symmetries', symmetries);                    
mean(dists_a)    