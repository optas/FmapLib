%%  A Script demonstrating the basic functionalities of the FmapLib.
%%  (Work in progress)

%% Load a Mesh and calculate basic quantities.
meshfile  = '../data/kid_rodola/0001.isometry.1.off';
inmesh    = Mesh(meshfile, 'rodola_1_1');
inmesh.set_triangle_angles();
inmesh.set_vertex_areas('barycentric');
inmesh

%% Calculate the first 100 spectra based on barycentric vertex areas.
LB             = Laplace_Beltrami(inmesh);
[evals, evecs] = LB.get_spectra(100, 'barycentric');

%% Compute WKS for log sample
variance = 6;
[E, sigma] = Mesh_Feature.energy_sample_generator('log_linear', evals(1), evals(end), 10, variance);

%% Save objects
    save('mesh_and_LB', 'refMesh', 'LB');

