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
                
                [train{end+1 : end+train_num}] = varargin{i}{train_i};                  % MATLAB requires multiple outputs to be enclosed between [].                
                [test{end+1  : end+test_num}]  = varargin{i}{test_i};
                
                if length(varargin) > 1
                 train_labels(end+1 : end+train_num) = i;
                 test_labels(end+1 : end+test_num)   = i;
                end
            end            
        end
        
        function [scores] = fmap_classify_naively(test_data, train_data, train_labels, maps)
            classes_num  = length(unique(train_labels));
            tests        = length(test_data);
            scores       = zeros(tests, classes_num);
            
            for test = 1:tests
                test_name = test_data{test};
                targets = maps(test_name);
                for key = targets.keys
                    target_name = key{:};
                    pos = strmatch(target_name, train_data, 'exact');
                    if ~isempty(pos)   % Only the training examples contribute the classification score.
                        target_class = train_labels(pos);                        
                        src_tar = targets(target_name);     % source to target fmap.
                        tar_src = maps(target_name);        
                        tar_src = tar_src(test_name);       % target to source fmap                        
                        M = src_tar.fmap * tar_src.fmap;
                        scores(test, target_class) = scores(test, target_class) + sum(sum(abs(M))) - sum(abs(diag(M)));                        
                    end
                end
                
            end
            
        end
        
        function [pairs] = inter_class_pairs(train_labels, train_data, class_id)
            class_examples = find(train_labels == class_id);
            names          = train_data(class_examples);            % Cell with the names of the class members.
            class_size     = length(names);
            
            if class_size < 2
                error('The class must have at least 2 members to form a pair.')
            end
            pairs_num = class_size * (class_size-1);                % (i,j) and (j,i) pairs are formed.
            pairs     = cell(pairs_num , 2);
            m = 1;
            for i=1:length(names)
                for j=1:length(names)
                    if i ~= j
                        pairs{m, 1} = names{i};
                        pairs{m, 2} = names{j};
                        m = m + 1;
                        
                    end
                end
            end
            assert(m == pairs_num + 1);
        end
            
            
        
        function [maps_optim, maps_per_iter] = low_rank_training_of_fmaps(all_maps, class_id, train_data, train_labels)        
            class_examples = find(train_labels == class_id);
            names          = train_data(class_examples);            % Cell with the names of the class members.
            init_maps      = cell(length(names), length(names));
            for i=1:length(names)
                idict = all_maps(names{i});
                for j=1:length(names)
                    if i ~=j
                        init_maps{i,j} = idict(names{j}).fmap;
                    end
                end
            end            
            W = 1/length(names) .* ones(length(names),length(names)) % TODO-P setting a constant W 
            [maps_optim, maps_per_iter] = Functional_Map.low_rank_filtering(init_maps, W);
            
        end
        
        
        

    end
    
end

