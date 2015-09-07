clr;
gitdir;
[dp, cp] = get_project_paths('ImageJointUnderstanding');
               
%% Load the image collection.

image_folder = [dp 'Rubinstein_CVPR_13/Data/Airplane100/'];
all_image_files = rdir([image_folder, [filesep '**' filesep]], 'regexp(name, ''\.jpg$'')');                
images = cell(length(all_image_files), 1);
for i=1:length(images)
    full_path       = all_image_files(i).name;
    path_substrings = strsplit(full_path, filesep); 
    last_word       = path_substrings{end};
    image_name      = last_word(1:end-4);               % Relying on the fact that length('.jpg') == length('.obj') == 4.                                                               
    images{i}       = Image(full_path, image_name);    
    gt              = [image_folder 'GroundTruth/' image_name '.png'];
    images{i}.set_gt_segmentation(imread(gt));    
end

%% Create image graphs and their Laplacians
image_graphs = containers.Map;
for i = 1:length(images)
    im_dims = sprintf('%d_%d',images{i}.height, images{i}.weight);
    if ~ image_graphs.isKey(im_dims)                       
        image_graphs(im_dims) = Laplacian(Image_Graph(images{i}, 'lattice'), 'norm');
    end    
end

%%
for im = image_graphs.values()
    im{:}.get_spectra(2)
end


%%
param.imageSize = [256 256];
param.orientationsPerScale = [8 8 8 8];
param.numberBlocks = 4;
param.fc_prefilt = 4;

LMgist(image_folder, '', param, image_folder)


%% Find gist-nearest neighbors per image.
k_neighbors = 2;
[gist_descriptors, gist_nn] = get_gist_nn(images, k_neighbors);

%% Extract features on pixel-level.
features   = cell(length(images), 1);
image_trio = [10, 34, 79];
for i = 1:length(image_trio)    
    features{i} = faceFeatures(images{image_trio(i)}, {'sift', 'lbp', 'color', 'hog'});
%     features{i} = faceFeatures(images{image_trio(i)}, {'lbp', 'color'});
end
%% Extract Basis
L = cell(1);
evecs = cell(1);
evals = cell(1);
neigs = 64;
image_trio = [10, 34, 79];
for i = 1:length(image_trio)    
    i
    im_id = image_trio(i);    
    im = images{im_id};         
    [h, w, ~] = size(im);
    
    F  = features{i};      % PP        %%Forgot to normalize probe fs
    Fp = reshape(F, w*h, size(F,3));
    
    G = simple_graphs('lattice', h, w);
    [row, col] = find(G);
    all_dist   = zeros(length(row), 1);
    for k = 1:length(row)            
        all_dist(k) = norm(Fp(row(k),:) - Fp(col(k),:));
    end
    
    sigma = 2*median(all_dist)^2;
    
    all_dist = exp(-all_dist.^2 ./ sigma);
    
    for k = 1:length(row)            
        G(row(k), col(k)) = all_dist(k);
    end    
    L{i} = adjacency_to_laplacian(G, 'comb');       % You should try 
    [evecs{i}, evals{i}] = eigs(L{i}, neigs, 'SA');
end

%%
save('just_in_case', 'image_trio', 'L', 'evecs', 'evals', 'features')
%%
    
%% Make F-maps.
    Fproj = cell(1);
    for i = 1:length(image_trio)            
        F  = features{i};
        F = reshape(F, size(F,1)*size(F,2), size(F,3));
        F = divide_columns(F, sqrt(sum(F.^2)));
        Fproj{i}  = evecs{i}' * F;
    end    
    
    n_ims = length(image_trio);
    Fmaps = cell(n_ims, n_ims);
    regulizer_w = 15;
    for i = 1:length(image_trio)            
        for j = 1:length(image_trio)         
            if i ~= j
                Fmaps{i,j} = Functional_Map.sum_of_frobenius_norms(Fproj{i}, Fproj{j}, diag(evals{i}), diag(evals{j}), regulizer_w);
            end
        end
    end
    
%% Get object Proposals
params = load('external/rp-master/config/rp.mat');  % Default Params for RP method.
params = params.params;
params.approxFinalNBoxes = 100;
proposals = cell(n_ims, 1);

for i = 1:length(image_trio)          
    im_id = image_trio(i);    
    im = images{im_id};         
    [h, w, ~] = size(im);
    proposals{i} = RP(im, params);
end

%%
overlaps{i} = cell(n_ims, 1);
for i = 1:length(image_trio)          
    im_id = image_trio(i);
    overlaps{i} = overlap_with_mask(proposals{i}, gts{im_id});       
end

%% Extract feautures of proposals and project them in corresponding basis.
proposal_features = cell(n_ims,1);
for i = 1:n_ims
    im_id = image_trio(i);    
    im = images{im_id};         
    [h, w, ~] = size(im);   
    prop_i = proposals{i};
    proposal_features{i} = cell(length(prop_i ),1);    
    frame = zeros(h, w, 387);
    
    for j = 1:length(proposals{i})        
        xmin = prop_i(j, 1);
        ymin = prop_i(j, 2);
        xmax = prop_i(j, 3);
        ymax = prop_i(j, 4);
        if xmin==xmax || ymin==ymax   
            continue
            'moufa'
        end
        
        patch = im(ymin:ymax, xmin:xmax, :);
%         image(uint8(patch))               
        yo = faceFeatures(patch, {'sift', 'color'});
        frame(ymin:ymax, xmin:xmax, :) = yo;        
        F = reshape(frame, w*h, size(frame, 3));
        F  = divide_columns(F, sqrt(sum(F.^2)));
        proposal_features{i}{j}  = evecs{i}' * F;
        frame(ymin:ymax, xmin:xmax, :) = 0;
    end
    
end


%% Now rank every pair of proposals (no triplets involved)
scores = cell(n_ims, n_ims);
for i = 1:n_ims
        for j = i+1:n_ims                        
                scores{i,j} = cell(length(proposal_features{i}), length(proposal_features{j}));
                for pi = 1:length(proposal_features{i})          % Proposals for image_i
                    for pj = 1:length(proposal_features{j})      % Proposals for image_j                        
                        if ~isempty(proposal_features{i}{pi}) &&  ~isempty(proposal_features{j}{pj})                      
                            Fs = proposal_features{i}{pi};
                            Ft = proposal_features{i}{pj};
                            misalignment = sum(sum(abs((Fmaps{i,j} * Fs) - Ft), 1));     %  ./ size(Fs, 2) % alignment error (L1 dist) of probe functions.       
                            misalignment = misalignment + sum(sum(abs((Fmaps{j,i} * Ft) - Fs), 1));
                            scores{i,j,pi,pj} = misalignment
                        end
                    end
                end
        end
end


%% Do simple stats

for i = 1:n_ims
        for j = i+1:n_ims                                        
                X = [];
                Y = [];
                for pi = 1:length(proposal_features{i})          % Proposals for image_i
                    for pj = 1:length(proposal_features{j})      % Proposals for image_j                        
                        if ~isempty(proposal_features{i}{pi}) &&  ~isempty(proposal_features{j}{pj})                                                  
                            
                            X(end+1) = scores{i,j,pi,pj};
                            Y(end+1) = overlaps{i}(pi) + overlaps{j}(pj);
                            
                        end
                    end
                end
                rho    = corr(X',Y')
                
                perc_5 = prctile(X, 5);
                [ind]  = find(X < perc_5);
                correspond_y = Y(ind);
                sum(Y < mean(correspond_y)) / length(Y)                
                
                figure;
                hist(X);                
                figure;
                hist(Y);                    
                
                [xval, xpos] = sort(X);                
                corres_y = Y(xpos(1:5));                
                                                
                [yval, ypos] = sort(Y);
                
                
                min(find(yval == max(corres_y))) ./ length(yval)
               
%                 input('')
                
        end
end





%% Partition each image into tiles.
grid_dim   = [20, 20];                                                        % Image will be partioned in grid_dim tiles.
grid_graph = simple_graphs('grid', grid_dim(1), grid_dim(2));

%% Extract BOW-feature.
bow_dict = load('/Users/optas/Dropbox/matlab_projects/Image_Graphs/data/Centers_MSRC');
bow_dict = bow_dict.Centers;

bow_feats   = cell(1); 
for i=1:total_images    
    [F, ~]       = bag_of_words(images(:,:,:,i), bow_dict, grid_dim(1), grid_dim(2));  % Extract Bag Of Word Features
    bow_feats{i} = F;
end

%%
    save('../data/output/aeroplane_first_20_bow_feats_100_100_grid', 'bow_feats')

%%
    load('../data/output/aeroplane_first_20_bow_feats', 'bow_feats')
%%
    F1 = bow_feats{1};

%% Derive weights of Laplacian.
    pd1 = pdist(F1);
    sigma = 2 * median(pd1(:))^2;
    pd1 = exp(-(pd1 .^2) / sigma );
    pd1 = squareform(pd1);
    m   = size(pd1, 1);
    pd1(1: m+1: end) = 1;
    imagesc(pd1)    
    pd1 = pd1 .* grid_graph;        % Try checker-board graph.    
%% Laplacian basis.
    neigs  = 64;
    Lg1    = diag(sum(pd1)) - pd1;
    [U, V] = eigs(Lg1, neigs, 'SA');
    Fproj  = U' * F1;
        
%     Freconstruct = U * Fproj;    
%     norm(abs(Freconstruct - F1), 'fro') / norm(F1, 'fro')
    
%%
