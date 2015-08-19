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

        function [low_rank_maps] = low_rank_filtering_of_fmaps(init_maps)                     % TODO-P Static if nothing is added.
            % Convert the init_maps a cell array that will be passed in Functional_Map.low_rank_filtering.
            index_dict   = string_keys_to_ints(init_maps);            
            n            = length(init_maps.keys);
            init_maps_as_cell = cell(n, n);
            for src_name = init_maps.keys                
                src_dict = init_maps(src_name{:});
                src_ind  = index_dict(src_name{:});
                for trg_name = src_dict.keys                    
                    init_maps_as_cell{src_ind, index_dict(trg_name{:})} = src_dict(trg_name{:}).fmap;
                end
            end            
            init_maps_as_cell
            % Call the optimization routine.
            W = repmat(1/n, n, n); % TODO-P setting a constant W.
            [maps_optim, ~] = Functional_Map.low_rank_filtering(init_maps_as_cell, W);                        
            
            % Write back the resulting cell array as a contrainers.Map
            low_rank_maps = containers.Map;
            for src_name = init_maps.keys                
                low_rank_maps(src_name{:}) = containers.Map;                                
                src_dict_in  = init_maps(src_name{:});                
                src_ind      = index_dict(src_name{:});
                temp         = low_rank_maps(src_name{:});
                for trg_name = src_dict_in.keys                                   
                    temp(trg_name{:}) = src_dict_in(trg_name{:}).copy();                    
                    temp(trg_name{:}).set_fmap((maps_optim{src_ind, index_dict(trg_name{:})}));     % TODO-P: Simplify.
                end
            end            
        end
        
        
        function [weights] = feature_weights_of_class(fmaps)
            weights = containers.Map;
            for src_name = fmaps.keys
                src_dict = fmaps(src_name{:});                
                for trg_name = src_dict.keys
                    M  = src_dict(trg_name{:});
                    Fs = M.source_basis.project_functions(M.source_neigs, M.source_features.F);
                    Fs = divide_columns(Fs, sqrt(sum(Fs.^2)));                                   % Normalize features.
                    Ft = M.target_basis.project_functions(M.target_neigs, M.target_features.F);               
                    Ft = divide_columns(Ft, sqrt(sum(Ft.^2)));                
                    if weights.isKey(src_name{:})                    
                        weights(src_name{:}) = weights(src_name{:})  + sum(abs((M.fmap * Fs) - Ft), 1); % TODO-P parameterize with a scoring function.
                    else
                        weights(src_name{:}) = zeros(1, size(Fs, 2));
                    end                    
                end
                weights(src_name{:}) = weights(src_name{:}) ./ size(src_dict, 1);
                
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
        
        

    end
    
end

