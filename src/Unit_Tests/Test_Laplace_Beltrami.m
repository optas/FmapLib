classdef Test_Laplace_Beltrami < matlab.unittest.TestCase
    % Unit test verifying the expected behavior and functionality of the
    % class 'Laplace_Beltrami'.
    %
    % Usage Example:
    %                test1 = Test_Laplace_Beltrami(); run(test1);
    %
    % Note: As of 26/05, checks the storing-returning of previously
    % computed spectra for only a small set of eigenvalue-values.

    properties
        shape_file  = '../../Data/kid_rodola/0001.isometry.1.off';
        ref_mesh    = [];
        LB          = [];
    end

    methods (TestClassSetup)
        % This code runs before any class instance is created.
        % See here: http://www.mathworks.com/help/matlab/matlab_prog/write-setup-and-teardown-code-using-classes.html
        function initialize_mesh_and_LB(obj)
            % Initializes the class instance to have a Mesh and its
            % associated LB.
            obj.ref_mesh  = Mesh(obj.shape_file);
            obj.LB        = Laplace_Beltrami(obj.ref_mesh);            
            eigs_num  = 200; area_type = 'barycentric';
            [~, ~] = obj.LB.get_spectra(eigs_num, area_type);
        end
    end

    methods (Test)
        function test_stored_spectra(obj)
            area_type = 'barycentric';
            eigs_num  = 200;
            [evals, evecs] = obj.LB.get_spectra(eigs_num, area_type);

            for i=1:5
                k = randi(eigs_num); % Number of eigs to be retrieved.
                [lambda, phi] = obj.LB.get_spectra(k, area_type);
                obj.verifyTrue(isequal(lambda, evals(1:k)), 'Retrieving precomputed eigenvalues.') ;
                obj.verifyTrue(isequal(phi, evecs(:,1:k)),  'Retrieving precomputed eigenvectors.');
            end
        end
        
        function test_project_functions(obj)            
            area_type = 'barycentric'; eigs_num  = 100;            
            [evals, evecs] = obj.LB.get_spectra(eigs_num, area_type);
            % Eigs does not give orthonormal vectors.
                     
            wks_samples    = 300; hks_samples    = 200;
            % Generate Mesh Feautures/Functions to project on LB basis.
            [energies, sigma] = Mesh_Features.energy_sample_generator('log_linear', evals(2), evals(end), wks_samples);
            wks_sig           = Mesh_Features.wave_kernel_signature(evecs(:,2:end), evals(2:end), energies, sigma);
            heat_time         = Mesh_Features.energy_sample_generator('log_sampled', evals(2), evals(end), hks_samples);
            hks_sig           = Mesh_Features.heat_kernel_signature(evecs(:,2:end), evals(2:end), heat_time);
            
            for i=1:5
                k = randi(eigs_num);     % Number of eigs to be retrieved.
                wf = randi(wks_samples); % Number of wks feautures to be used.
                hf = randi(hks_samples); % Number of hks feautures to be used.                

                res1              = obj.LB.project_functions(area_type, k, wks_sig(:, 1:wf), hks_sig(:, 1:hf));                                
                obj.verifyTrue(size(res1, 1) == k);
                obj.verifyTrue(size(res1, 2) == wf+hf);                   
                res2              = evecs(:, 1:k) \ [wks_sig(:, 1:wf) hks_sig(:, 1:hf)];
                obj.verifyTrue(all_close(res1,res2));

            end
        end
    end

end

