%%  A Script demonstrating the basic functionalities of the FmapLib.
%%  (Work in progress)

%% Load a Mesh and calculate basic quantities.
meshfile  = '../data/kid_rodola/0001.isometry.1.off';
% meshfile  = '/Users/optas/Dropbox/Matlab_Projects/3D_Meshes/Data/kid_rodola/0001.isometry.1.off';
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



