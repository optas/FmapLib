classdef Test_Graph < matlab.unittest.TestCase
    % Unit test verifying the expected behavior and functionality of the
    % class 'Graph'.
    %
    % Usage Example:
    %                test1 = Test_Graph(); run(test1);
    %
  
    properties
        gtypes = {'lattice', 'clique'};
        lattice
        clique            
    end

    methods (TestClassSetup)
        % This code runs before any class instance is created.
        % See here: http://www.mathworks.com/help/matlab/matlab_prog/write-setup-and-teardown-code-using-classes.html
        function initialize_some_graph_models(obj)            
            obj.lattice = Graph.generate('lattice', 10, 10);
            obj.clique  = Graph.generate('clique', 100);            
        end
    end

    methods (Test)
        function test_incidence_matrix(obj)            
            for g = obj.gtypes
                gname = g{:};
                I = obj.(gname).incidence_matrix();
                Li = I*I'; % Incidence derived Laplacian.
                Lap = Laplacian(obj.(gname), 'comb');
                obj.verifyTrue(all(all(Lap.L == Li)));                        
            end
        end
               
    end

end

