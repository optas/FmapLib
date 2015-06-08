%%  A Script demonstrating the basic functionalities of the FmapLib (Work in progress).
gitdir;
cd FmapLib/src
%% Load a Mesh and calculate basic quantities.
meshfile  = '../data/kid_rodola/0001.isometry.1.off';
inmesh    = Mesh(meshfile, 'rodola_1_1');
inmesh.set_triangle_angles();
inmesh.set_vertex_areas('barycentric');
inmesh

%% Calculate the first 100 spectra, based on barycentric vertex areas.
LB             = Laplace_Beltrami(inmesh);
[evals, evecs] = LB.get_spectra(100, 'barycentric');
save('../data/output/mesh_and_LB', 'inmesh', 'LB');

%% Load Precomputed ones.
load('../data/output/mesh_and_LB', 'inmesh', 'LB');
[evals, evecs] = LB.get_spectra(100, 'barycentric');

%%
pairs = [1,2; 1,55; 1,100]';
% pairs must be passed as 2 x N
tic
D1 = comp_geodesics_pairs(inmesh.vertices(:,1), inmesh.vertices(:,2), inmesh.vertices(:,3), inmesh.triangles', pairs);
toc  
%%
%%
sources = [1];
tic
D2 = comp_geodesics_to_all(inmesh.vertices(:,1), inmesh.vertices(:,2), inmesh.vertices(:,3), inmesh.triangles', sources);
toc