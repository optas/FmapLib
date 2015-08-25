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