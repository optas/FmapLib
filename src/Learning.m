classdef Learning
    % Machine Learning Toolbox
    
    
    methods (Static)
        
        function [train, test, train_labels, test_labels] = sample_class_observations(test_percent, varargin)            
            train = {};
            test  = {};
            train_labels = [];
            test_labels  = [];
            for i=1:length(varargin)
                class_size = length(varargin{i});
                test_i     = randsample(class_size, ceil(test_percent * class_size));   % Vector with indexes to be included as training.
                train_i    = ~ismember(1:class_size, test_i);                           % Logical vector [0, 0, 1, 1].                
                test_num   = length(test_i);
                train_num  = sum(train_i);
                
                [train{end+1 : end+train_num}] = varargin{i}{train_i};               % MATLAB requires multiple outputs to be enclosed between [].                
                [test{end+1  : end+test_num}]  = varargin{i}{test_i};
                
                if length(varargin) > 1
                 train_labels(end+1 : end+train_num) = i;
                 test_labels(end+1 : end+test_num)   = i;
                end
            end            
        end

    end
    
end

