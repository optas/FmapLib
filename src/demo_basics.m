%%  A Script demonstrating the basic functionalities of the FmapLib.
%%  (Work in progress)

%% Load a Mesh and calculate basic quantities.
meshfile  = '../data/kid_rodola/0001.isometry.1.off';
% meshfile  = '/Users/optas/Dropbox/Matlab_Projects/3D_Meshes/Data/kid_rodola/0001.isometry.1.off';
inmesh    = Mesh(meshfile, 'rodola_1_1');
inmesh.set_triangle_angles();
inmesh.set_vertex_areas('barycentric');
inmesh

%% Test vertex normal
inmesh.set_vertex_normals();
N = inmesh.vertex_normal;

figure;
trisurf(inmesh.triangles, inmesh.vertices(:,1), inmesh.vertices(:,2), inmesh.vertices(:,3), 1);
hold on;
quiver3(inmesh.vertices(:,1), inmesh.vertices(:,2), inmesh.vertices(:,3), N(:,1), N(:,2), N(:,3));
hold off;
axis equal;

%% Mean curvature
mean_curv = Mesh_Features.mean_curvature(inmesh, Laplace_Beltrami(inmesh), [0, 0.5, 2, 10]);

figure;
subplot(2, 2, 1); trisurf(inmesh.triangles, inmesh.vertices(:,1), inmesh.vertices(:,2), inmesh.vertices(:,3), mean_curv(:,1));
shading interp; axis equal; title('Mean Curvature');
subplot(2, 2, 2); trisurf(inmesh.triangles, inmesh.vertices(:,1), inmesh.vertices(:,2), inmesh.vertices(:,3), mean_curv(:,2));
shading interp; axis equal; title('Mean Curvature Smoothed');
subplot(2, 2, 3); trisurf(inmesh.triangles, inmesh.vertices(:,1), inmesh.vertices(:,2), inmesh.vertices(:,3), mean_curv(:,3));
shading interp; axis equal; title('Mean Curvature Smoothed');
subplot(2, 2, 4); trisurf(inmesh.triangles, inmesh.vertices(:,1), inmesh.vertices(:,2), inmesh.vertices(:,3), mean_curv(:,4));
shading interp; axis equal; title('Mean Curvature Smoothed');

%% Gauss curvature
gauss_curv = Mesh_Features.gaussian_curvature(inmesh, [0, 20, 200, 500]);

figure;
subplot(2, 2, 1); trisurf(inmesh.triangles, inmesh.vertices(:,1), inmesh.vertices(:,2), inmesh.vertices(:,3), gauss_curv(:,1));
shading interp; axis equal; title('Gauss Curvature');
subplot(2, 2, 2); trisurf(inmesh.triangles, inmesh.vertices(:,1), inmesh.vertices(:,2), inmesh.vertices(:,3), gauss_curv(:,2));
shading interp; axis equal; title('Gauss Curvature Smoothed');
subplot(2, 2, 3); trisurf(inmesh.triangles, inmesh.vertices(:,1), inmesh.vertices(:,2), inmesh.vertices(:,3), gauss_curv(:,3));
shading interp; axis equal; title('Gauss Curvature Smoothed');
subplot(2, 2, 4); trisurf(inmesh.triangles, inmesh.vertices(:,1), inmesh.vertices(:,2), inmesh.vertices(:,3), gauss_curv(:,4));
shading interp; axis equal; title('Gauss Curvature Smoothed');

%% Calculate the first 100 spectra, based on barycentric vertex areas.
LB             = Laplace_Beltrami(inmesh);
[evals, evecs] = LB.get_spectra(100, 'barycentric');
save('../data/output/mesh_and_LB', 'inmesh', 'LB');

%% Load Precomputed ones.
load('../data/output/mesh_and_LB', 'inmesh', 'LB');
[evals, evecs] = LB.get_spectra(100, 'barycentric');

%% Shape Descriptors
 % WKS - HKS - GPS

nsamples = 20;
[energies, sigma] = Mesh_Features.energy_sample_generator('log_linear', evals(2), evals(end), nsamples);
wks_sig           = Mesh_Features.wave_kernel_signature(evecs(:,2:end), evals(2:end), energies, sigma);
wks_aubrey_sig    = Mesh_Features.wks_Aubry(evecs(:,2:end), evals(2:end), energies, sigma);                 % TODO-P Add to Unit Test.
all_close(wks_sig, wks_aubrey_sig)

% 
% energies          = Mesh_Features.energy_sample_generator('log_sampled', evals(2), evals(end), nsamples);
% hks_sig           = Mesh_Features.heat_kernel_signature(evecs(:,2:end), evals(2:end), energies);
% v_areas           = inmesh.get_vertex_areas('barycentric');
% hks_sun           = Mesh_Features.hks_Sun(evecs, evals, v_areas, energies, 1);
% all_close(hks_sig, hks_sun)
% 
% hks_sig           = Mesh_Features.heat_kernel_signature(evecs(:,2:end), evals(2:end), energies);
% v_areas           = ones(inmesh.num_vertices, 1);
% hks_sun           = Mesh_Features.hks_Sun(evecs, evals, v_areas, energies, 1);
% all_close(hks_sig, hks_sun)
% 
% gps_sig           = Mesh_Features.global_point_signature(evecs(:,2:end), evals(2:end));

%% D2 Testing
