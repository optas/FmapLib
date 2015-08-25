% classdef Test_Functional_Map < matlab.unittest.TestCase
%     % Unit test verifying the expected behavior and functionality of the
%     % class 'Laplace_Beltrami'.
%     %
%     % Usage Example:
%     %                test1 = Test_Laplace_Beltrami(); run(test1);
%     %
%     % Note: As of 26/05, checks the storing-returning of previously
%     % computed spectra for only a small set of eigenvalue-values.
% 
%     properties
%         shape_file  = '../../Data/kid_rodola/0001.isometry.1.off';
%         ref_mesh    = [];
%         LB          = [];
%     end
% 
%     methods (TestClassSetup)
%         % This code runs before any class instance is created.
%         % See here: http://www.mathworks.com/help/matlab/matlab_prog/write-setup-and-teardown-code-using-classes.html
%         function initialize_mesh_and_LB(obj)
%             % Initializes the class instance to have a Mesh and its
%             % associated LB.
%             obj.ref_mesh  = Mesh(obj.shape_file);
%             obj.LB        = Laplace_Beltrami(obj.ref_mesh);            
%             eigs_num  = 200; area_type = 'barycentric';
%             [~, ~] = obj.LB.get_spectra(eigs_num, area_type);
%         end
%     end
% 
%     methods (Test)
%         function test_stored_spectra(obj)
%             area_type = 'barycentric';
%             eigs_num  = 200;
%             [evals, evecs] = obj.LB.get_spectra(eigs_num, area_type);
% 
%             for i=1:5
%                 k = randi(eigs_num); % Number of eigs to be retrieved.
%                 [lambda, phi] = obj.LB.get_spectra(k, area_type);
%                 obj.verifyTrue(isequal(lambda, evals(1:k)), 'Retrieving precomputed eigenvalues.') ;
%                 obj.verifyTrue(isequal(phi, evecs(:,1:k)),  'Retrieving precomputed eigenvectors.');
%             end
%         end



% %%  Add to unit-test    
%     gt_map            = (1:mesh1.num_vertices)';   % Ground truth correspondences from Source_Mesh to Target_Mesh.    
%     [~, source_basis] = LB1.get_spectra(neigs);
%     [~, target_basis] = LB2.get_spectra(neigs);  
%     X_opt             = Functional_Map.groundtruth_functional_map(source_basis, target_basis, gt_map, mesh2.get_vertex_areas('barycentric'));  %TODO-P, Ugly.
%     % Evaluate the X_opt
%     eval_points = 200;
%     [dists,  random_points]  = Functional_Map.pairwise_distortion_of_map(X_opt, LB1, LB2, gt_map, 'nsamples', eval_points, 'symmetries', symmetries);
%     mean(dists)
%     hist(dists)
%     
% 
% %     [dists2, random_points2] = Functional_Map.pairwise_distortion_of_map(X_opt, LB1, LB2, gt_map, 'indices', random_points);
% %     assert(length(dists) == eval_points);
% %     assert(all(dists==dists2) && all(random_points == random_points2));
% %     length(unique(random_points)) == length(random_points)
%     %%
%     % Use symmetries.    
%     C = textread('../data/input/kid_rodola/sym.txt', '%s', 'delimiter', ' ');  % Read symmetries:
%     C = cell2mat(C); sym = str2num(C);            
%        
%     %%
%     [dists3, random_points3] = Functional_Map.pairwise_distortion_of_map(X_opt, mesh1, mesh2, source_basis, target_basis, gt_map, 'indices', random_points, 'symmetries', sym);
%     assert(all(dists3 <= dists2));
% 
%     %% Use a small number of eigenvectors to do the Fmap.
%     [~, source_basis] = LB1.get_spectra(10, 'barycentric');
%     [~, target_basis] = LB2.get_spectra(10, 'barycentric');    
%     X_opt_small       = Functional_Map.groundtruth_functional_map(source_basis, target_basis, gt_map, mesh2.get_vertex_areas('barycentric'));
%     [dists2, random_points2] = Functional_Map.pairwise_distortion_of_map(X_opt_small, mesh1, mesh2, source_basis, target_basis, gt_map, 'indices', random_points, 'fast');
%     assert(mean(dists2) > mean(dists)) % TODO-P Change to something more reasonable.
%     
%     