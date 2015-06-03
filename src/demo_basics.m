%%  A Script demonstrating the basic functionalities of the FmapLib.
%%  (Work in progress)

%% Load a Mesh and calculate basic quantities.
% meshfile  = '../data/kid_rodola/0001.isometry.1.off';
meshfile  = '/Users/optas/Dropbox/Matlab_Projects/3D_Meshes/Data/kid_rodola/0001.isometry.1.off';
inmesh    = Mesh(meshfile, 'rodola_1_1');
inmesh.set_triangle_angles();
inmesh.set_vertex_areas('barycentric');
inmesh

%% Calculate the first 100 spectra, based on barycentric vertex areas.
LB             = Laplace_Beltrami(inmesh);
[evals, evecs] = LB.get_spectra(100, 'barycentric');
% Save objects
save('../data/output/mesh_and_LB', 'inmesh', 'LB');


%% WKS - HKS
nsamples = 20;
% [energies, sigma] = Mesh_Features.energy_sample_generator('log_linear', evals(2), evals(end), nsamples);
% wks_sig           = Mesh_Features.wave_kernel_signature(evecs(:,2:end), evals(2:end), energies, sigma);


energies          = Mesh_Features.energy_sample_generator('log_sampled', evals(2), evals(end), nsamples)
hks_sig           = Mesh_Features.heat_kernel_signature(evecs(:,2:end), evals(2:end), energies);
v_areas           = inmesh.get_vertex_areas('barycentric');
hks_sun           = Mesh_Features.hks_Sun(evecs, evals, v_areas, energies, 1);
allclose(hks_sig, hks_sun)

hks_sig           = Mesh_Features.heat_kernel_signature(evecs(:,2:end), evals(2:end), energies);
v_areas           = ones(inmesh.num_vertices, 1);
hks_sun           = Mesh_Features.hks_Sun(evecs, evals, v_areas, energies, 1);
allclose(hks_sig, hks_sun)

% wks_aubrey_sig    = Mesh_Features.wks_aubrey(evecs(:,2:end), evals(2:end), energies, sigma);  % Add to Unit Test.
% allclose(wks_sig, wks_aubrey_sig)
