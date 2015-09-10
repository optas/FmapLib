% seg.coords = RP(img, params_rp); %[xmin, ymin, xmax, ymax]
% seg.coords = [seg.coords; 1, 1, size(img,2), size(img,1)]; % add a whole box
% standarizeImage
% see vis_CorLoc and there in e.g., the loadView_seg

clr;
[dp, cp] = get_project_paths('ImageJointUnderstanding');


%% Load Rubinstein image collection.
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


%% Resize every image 
new_height = 64;
new_weight = 64;
num_eigs   = 32;
images_laplacians = cell(length(all_image_files), 1);

for i=1:length(images)
    images{i}.set_resized_image(new_height, new_weight);
           
    Fi = faceFeatures(images{i}.resized, {'color'});  % normalize probe fs
    [h, w, ~] = size(images{i}.resized);
    
    G  = Graph.generate('lattice', h, w);               % P fix non_writtable properties.
    
    [node_from, node_to] = find(G.A);                   % all edges between nodes (double counting)

    from_row = ceil(node_from / double(w));             % next four convert row-expanded nodes of pixel matrix, to 2d (i,j) indices.
    to_row   = ceil(node_to   / double(w));
    from_col = node_from - ((from_row - 1) * w );
    to_col   = node_to   - ((to_row - 1) * w );

    G.num_edges * 2 == length(node_from)   % TODO take care for bidir.
    edges = length(node_from);

    all_dist  = zeros(edges, 1);
    parfor k = 1:edges                       % traverse all edges
        all_dist(k) = norm(squeeze(Fi(from_row(k), from_col(k), :)) - squeeze(Fi(to_row(k), to_col(k), :)));      % distances based on feautures (Fi)   
    end

    sigma = 2*median(all_dist)^2           % convert distances to similarities
    all_dist = exp(-all_dist.^2 ./ sigma);
    
    A = G.A;

    for k = 1:edges
        A(node_from(k), node_to(k)) = all_dist(k);    
    end

    G = Graph(A, 0);
    Li_w = Laplacian(G, 'comb');
    Li_w.get_spectra(num_eigs);
    
    images_laplacians{i} = Li_w;    
end
    
%% Create image graphs, their Laplacians and compute their spectra.
image_graphs = containers.Map;
eigs_num     = 64;
for i = 1:length(images)
    im_dims = sprintf('%d_%d',images{i}.height, images{i}.weight);
    if ~ image_graphs.isKey(im_dims)                       
        image_graphs(im_dims) = Laplacian(Image_Graph(images{i}, 'lattice'), 'norm');
        image_graphs(im_dims).get_spectra(eigs_num);
    end    
end
save('voc_airplane_left_lattice_with_64_spectra.mat', 'image_graphs')
%%
load('voc_airplane_left_lattice_with_64_spectra.mat')

%% Laplacians that take into account the RGB diffrences as weights.
i = 1;
im_dims = sprintf('%d_%d', images{i}.height, images{i}.weight);
Li = image_graphs(im_dims);    
Fi = faceFeatures(images{i}.CData, {'color'});  % normalize probe fs
w = images{i}.weight;
h = images{i}.height;
[node_from, node_to] = find(Li.G.A);                   % all edges between nodes (double counting)

from_row = ceil(node_from / double(w));             % next four convert row-expanded nodes of pixel matrix, to 2d (i,j) indices.
to_row   = ceil(node_to   / double(w));
from_col = node_from - ((from_row - 1) * w );
to_col   = node_to   - ((to_row - 1) * w );

Li.G.num_edges *2 == length(node_from)   % TODO take care for bidir.
edges = length(node_from);

all_dist  = zeros(edges, 1);
for k = 1:edges                       % traverse all edges
    all_dist(k) = norm(squeeze(Fi(from_row(k), from_col(k), :)) - squeeze(Fi(to_row(k), to_col(k), :)));      % distances based on feautures (Fi)   
end

sigma = 2*median(all_dist)^2           % convert distances to similarities
all_dist = exp(-all_dist.^2 ./ sigma)

G = Graph.generate('lattice', h, w);                % P fix non_writtable properties.
A = G.A;

for k = 1:edges
    A(node_from(k), node_to(k)) = all_dist(k);    
end

G = Graph(A, 0);
Li_w = Laplacian(G, 'norm');

Li_w.get_spectra(3);

i = 1;
for eigs_num = [1,2,5,10,30,45,62];

    im_dims = sprintf('%d_%d', images{i}.height, images{i}.weight);
%     Li = image_graphs(im_dims);    
    Li = Li_w;

    Fi = faceFeatures(images{i}.CData, {'color'});

    [h, w, c] =  size(images{i}.CData);
    Fp = reshape(permute(Fi,[2,1,3]), h*w, c);   % Now Fp, contains the pixel-level features stack in a long vector row-wise.


    compressed    = Li.project_functions(eigs_num, Fp);
    reconstructed = Li.synthesize_functions(compressed);

    accu = mean(divide_columns(sqrt(sum(abs(Fp - reconstructed).^2)), sqrt(sum(Fp.^2))))   % percent in loss of accurasy. TODO-P
    
    reconstructed = reshape(reconstructed, w, h, c);
    reconstructed = permute(reconstructed, [2,1,3]);

    im_r_i        = Image(reconstructed);        
    f = im_r_i.plot;        
    title('Plotting reconstructed RGB. 4NN Laplacian basis with RGB weights.')
    xlabel(sprintf('Number of Eigenvalues: %d', eigs_num))
    ylabel(sprintf('Mean accuracy: %f', accu))
    saveas(f, sprintf('RGB_weighted_%d.png', eigs_num))    
end

%% Derive hogs for every image.
eigs_num  = 64;
hogs      = cell(length(all_image_files), 1);
hogs_proj = cell(length(all_image_files), 1);

for i=1:length(images)
    hogs{i} = faceFeatures(images{i}.CData, {'hog'});

    [h, w, c] =  size(images{i}.CData);
    
    Fp = reshape(permute(hogs{i},[2,1,3]), h*w, size(hogs{i}, 3));   % Now Fp, contains the pixel-level features stack in a long vector row-wise.
    Fp = divide_columns(Fp, sqrt(sum(Fp.^2)));                       % Rescale to unit norm.    
    
    im_dims = sprintf('%d_%d',images{i}.height, images{i}.weight);
    Li = image_graphs(im_dims);    
    hogs_proj{i} = Li.project_functions(eigs_num, Fp);

    reconstructed = Li.synthesize_functions(hogs_proj{i});
    accu = mean(sqrt(sum(abs(Fp - reconstructed).^2)))
end


%% Make F-maps.
    Fproj = cell(1);
    for i = 1:length(image_trio)            
        F  = features{i};
        F  = reshape(F, size(F,1)*size(F,2), size(F,3));
        F  = divide_columns(F, sqrt(sum(F.^2)));
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
    



image(im2uint8(im_r_i.CData))
% im_r_i.plot()
%%
            











%% Find gist-nearest neighbors per image.
param.imageSize = [256 256];
param.orientationsPerScale = [8 8 8 8];
param.numberBlocks = 4;
param.fc_prefilt = 4;
LMgist(image_folder, '', param, image_folder)

k_neighbors = 2;
[gist_descriptors, gist_nn] = get_gist_nn(images, k_neighbors);
    

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
