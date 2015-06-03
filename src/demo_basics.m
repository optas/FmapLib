%%  A Script demonstrating the basic functionalities of the FmapLib.
%%  (Work in progress)

%% Load a Mesh and calculate basic quantities.
meshfile  = '../data/kid_rodola/0001.isometry.1.off';
% shapeFile  = '/Users/optas/Dropbox/Matlab_Projects/3D_Meshes/Data/kid_rodola/0001.isometry.1.off';
inmesh    = Mesh(meshfile, 'rodola_1_1');
inmesh.set_triangle_angles();
inmesh.set_vertex_areas('barycentric');
inmesh

%% Calculate the first 100 spectra, based on barycentric vertex areas.
LB             = Laplace_Beltrami(inmesh);
[evals, evecs] = LB.get_spectra(100, 'barycentric');
% Save objects
% save('../data/output/mesh_and_LB', 'inmesh', 'LB');


%% WKS - HKS
nsamples = 100;
[energies, sigma] = Mesh_Features.energy_sample_generator('log_linear', evals(2), evals(end), nsamples);

wks_sig           = Mesh_Features.wave_kernel_signature(evecs(:,2:end), evals(2:end), energies, sigma);

hks_sig = Mesh_Features.heat_kernel_signature(evecs(:,2:end), evals(2:end), energies);



% wks_aubrey_sig    = Mesh_Features.wks_aubrey(evecs(:,2:end), evals(2:end), energies, sigma);  % Add to Unit Test.
% allclose(wks_sig, wks_aubrey_sig)
