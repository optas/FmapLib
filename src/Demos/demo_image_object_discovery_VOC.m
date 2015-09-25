clear; clc;
[dp, cp]   =  get_project_paths('ImageJointUnderstanding');
            
%% Load VOC image collection with gist, gt_segmentation and object_proposals.
image_folder    = [dp 'VOC/2007/Train_Validate_Data/VOC2007_6x2/aeroplane_left/'];
all_image_files = rdir([image_folder, [filesep '**' filesep]], 'regexp(name, ''\.jpg$'')');
images          = cell(length(all_image_files), 1);
obj_proposals   = cell(length(all_image_files), 1);
gist_signatures = cell(length(all_image_files), 1);
for i=1:length(images)
    full_path         = all_image_files(i).name;
    path_substrings   = strsplit(full_path, filesep); 
    last_word         = path_substrings{end};
    image_name        = last_word(1:end-4);                             % Relying on the fact that length('.jpg') == length('.obj') == 4.                                                               
    images{i}         = Image(full_path, image_name);    
    
    gt_segmentation   = load([image_folder image_name '.mat']);
    gt_segmentation   = gt_segmentation.bbox_list;
    images{i}.set_gt_segmentation(gt_segmentation);
    
    obj_proposals{i}  = load([image_folder image_name '_seg.mat']);    
    obj_proposals{i}  = int16(obj_proposals{i}.feat.boxes);
  
    gist_signatures{i}  = load([image_folder image_name '_gist.mat']);    
    gist_signatures{i}  = gist_signatures{i}.gist_feat;
end
num_images = length(images);

%% Load pixel-leve ground-truth segmentations
pixel_gt_segmentation = cell(num_images, 1);
gt_seg_folder = [dp 'VOC/2007/Train_Validate_Data/SegmentationClass/'];
ending        = '.png';
class_id      =  1;    % Specific to aeroplanes - hack.
for i = 1:num_images
    try
        gt = imread([gt_seg_folder images{i}.name ending]);
    catch
        [gt_seg_folder images{i}.name ending]
        continue
    end
        
    mask = zeros(size(gt));
    mask(gt == class_id) = 1;
    pixel_gt_segmentation{i} = mask;    
end

%% Resize every image
new_height = 128;
new_width  = NaN;
for i= 1:num_images
    images{i}.set_resized_image(new_height, new_width);
end

%% Extract a laplacian basis.
radius     = 3;
eigs_num   = 64;
sigma_f    = 800 / (255^2);

image_laplacians = cell(length(all_image_files), 1);

for i = 11:num_images
    im_i       = images{i}.get_resized_image();
    new_heigtt = im_i.height;
    new_width  = im_i.width;
        
    G = Image_Graph(im_i, 'r_radius_connected',  radius);
    fprintf('Graph %d constructed.\n', i);
    
    Fi = im_i.color();
    sigma_s  =  2 * (0.1 * norm([new_width-1, new_height-1]))^2;    
    G.adjust_weights_via_feature_differences(Fi , 'normalized_cut', 'sigma_s', sigma_s, 'sigma_f', sigma_f);
    
    image_laplacians{i} = Laplacian(G, 'norm');
    image_laplacians{i}.get_spectra(eigs_num);
    fprintf('Laplacian %d constructed.\n', i);
end
%%
save([dp 'output/laplacians_aeroplane_left_norm'], 'image_laplacians');
%%
load([dp 'output/laplacians_aeroplane_left_norm'], 'image_laplacians');
%% Visualize Eigenvectors of single image.
im_id = 10; eig_to_vis = 20;
E = image_laplacians{im_id}.evecs(eigs_num);   
new_width  = images{im_id}.get_resized_image().width;
new_height = images{im_id}.get_resized_image().height;
figure();

for i = 1:eig_to_vis     
    eig_i = E(:, i);
    eig_i = reshape(eig_i, new_height, new_width);
%     imagesc(eig_i); axis('image'); axis off;    
    eig_i = (eig_i - min(eig_i(:))) ./ (max(eig_i(:)) - min(eig_i(:)));    % imshow expects values in [0,1] if they are floats.
    subplot(5,4,i); imshow(eig_i);    
end
suptitle('20 Smallest Eigenvectors (Min-Cut method).')    

%% Extract hog features: whole-image, (pixel-wise) and project them into Laplacian basis.
clc
hog_feats = cell(num_images ,1);
proj_hogs = cell(num_images, 1);

for i= 1:num_images
    im_i = images{i}.get_resized_image();   
    h = im_i.height; w = im_i.width;

%     hog_feats{i} = Image_Features.hog_signature(im_i);    
    hog_feats{i} = Image_Features.sift_signature(im_i);  
    
    Fp = reshape(hog_feats{i}, h*w, size(hog_feats{i}, 3));       % Now Fp, contains the pixel-level features stack in a long vector column-wise.
    Fp = divide_columns(Fp, sqrt(sum(Fp.^2)));                    % Rescale to unit norm.           
    proj_hogs{i}  = image_laplacians{i}.project_functions(eigs_num, Fp);    
    
%     reconstructed = image_laplacians{i}.synthesize_functions(proj_hogs{i}(:,1));
%     error_i = l2_norm(Fp(:,1) - reconstructed) ./ l2_norm(Fp(:,1));
%     imagesc(reshape(reconstructed, h, w));    
end
%%
% load([dp 'output/hogs_aeroplane_left'], 'hog_feats', 'proj_hogs');

%% Make all F-maps.
all_fmaps = cell(num_images, num_images);
regulizer_w = 15;
for i = 1:num_images
    evals_i = image_laplacians{i}.evals(eigs_num);    
%     evals_i = evals_i ./ max(evals_i);
    
    for j = 1:num_images
        if i ~= j
            evals_j = image_laplacians{j}.evals(eigs_num);            
%             evals_j = evals_j ./ max(evals_j);
            all_fmaps{i,j} = Functional_Map.sum_of_squared_frobenius_norms(proj_hogs{i}, proj_hogs{j}, evals_i, evals_j, regulizer_w);
        end        
    end
end

%% Make Iterative-Fmaps.
iter_fmaps = cell(num_images, num_images);
regulizer_w = 15;
for i = 1:2
    evals_i = image_laplacians{i}.evals(eigs_num);
%     evals_i = evals_i ./ max(evals_i);
    for j = 1:4
        if i ~= j
            evals_j = image_laplacians{j}.evals(eigs_num);                         
%             evals_j = evals_j ./ max(evals_j);
            [X, W, iter_ran]  = Functional_Map.iteratively_refined_fmap(proj_hogs{i}, proj_hogs{j}, evals_i, evals_j, regulizer_w);
            iter_fmaps{i,j} = X(i,j,end);
        end        
    end
end

%%
i = 6
h = images{i}.height;
w = images{i}.width;

images{i}.plot
plotBoxes(images{i}.gt_segmentation{1}, 'random', [], '-')

mask = pixel_gt_segmentation{i};
imshow(mask)

%%
for gt_j = images{i}.gt_segmentation                       
    gt = gt_j{:};
    bitmask = zeros(h, w);    
    bitmask(gt(2):gt(4), gt(1):gt(3)) = 1;     
    bitmask = mask .* bitmask;
    
%     overlaps(i,:,m) = max(overlaps(i,:,m), overlap_with_mask(obj_proposals{i}(top_patches(i, :, m), :), bitmask)');        
end

imshow(bitmask)
%%










% load([dp 'output/fmaps_aeroplane_left'], 'all_fmaps');

%% k-nn GIST initial image network
k_neighbors = 5;
gist_desc = cell2mat(gist_signatures);
[nns, gist_dists] = knnsearch(gist_desc, gist_desc, 'k', k_neighbors+1);
nns = nns(:, 2:k_neighbors+1);
gist_dists = gist_dists(:, 2:k_neighbors+1);


%% Extract and project features for every patch
method = 'zero_pad' ; % 'tiling'
patch_feat = cell(num_images, 1);
for i = 1:num_images    
    o_i  = obj_proposals{i};
    patch_feat{i} = cell(length(o_i), 1);
    
    im_i = images{i}.get_resized_image();        
    new_height = im_i.height;
    new_width  = im_i.width;
    
    for pi = 1:length(o_i)                   % Proposals for image_i
        if ~ Patch.are_valid_corners(o_i(pi, :),  images{i}) 
            continue;
        end
        corners_pi = Patch.find_new_corners(images{i}.height, images{i}.width, new_height, new_width, o_i(pi,:));        
        Fs = Patch.extract_patch_features(corners_pi , hog_feats{i}, method);                        
        Fs = reshape(Fs, new_height*new_width, size(Fs, 3));                              
        Fs = divide_columns(Fs, max(sqrt(sum(Fs.^2)), 1e-20));                    % Rescale to unit norm.    
        Fs = image_laplacians{i}.project_functions(eigs_num, Fs);
        patch_feat{i, pi} = Fs;
    end
end

load([dp 'output/patch_feats_' method], 'patch_feat')







%% Find top patches
top_p = 5;
number_iters = 5;
[new_nns, top_patches] = Patch.iterative_compute_patches(nns, patch_feat, all_fmaps, top_p, number_iters);

%%

save([dp 'output/top_patches_and_nns_tiling_minimum_square_rule'])

%% Evaluation 
overlaps = zeros(num_images, top_p, number_iters);
no_gt = 0;
for m = 1:number_iters
    for i = 1:num_images
        h = images{i}.height;
        w = images{i}.width;
        for gt_j = images{i}.gt_segmentation            
            if isempty(gt_j)
                no_gt = no_gt + 1;
                continue
            end        
            gt = gt_j{:};                
            bitmask = zeros(h, w);
            bitmask(gt(2):gt(4), gt(1):gt(3)) = 1;  
            overlaps(i,:,m) = max(overlaps(i,:,m), overlap_with_mask(obj_proposals{i}(top_patches(i, :, m), :), bitmask)');        
        end           
    end
end

save([dp 'output/top_overlapss_tiling_minimum_square_rule'], 'overlaps')
%%






%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                                          Basement
%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%% Get object Proposals
params = load([cp 'Randomized_Prim_Object_Proposal/config/rp.mat']);  % Default Params for RP method.
params = params.params;

params.approxFinalNBoxes = 100;    

proposals = cell(num_images, 1);

for i = 1:num_images
    im_id = image_trio(i);    
    im = images{im_id};         
    [h, w, ~] = size(im);
    proposals{i} = RP(im, params);
end

%%
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
    accu = mean(sqrt(sum(abs(Fp - reconstructed).^2)));
end


%% Find gist-nearest neighbors per image.
param.imageSize = [256 256];
param.orientationsPerScale = [8 8 8 8];
param.numberBlocks = 4;
param.fc_prefilt = 4;
LMgist(image_folder, '', param, image_folder)
k_neighbors = 2;
[gist_descriptors, gist_nn] = get_gist_nn(images, k_neighbors);
    

%% Extract BOW-feature.
bow_dict = load('/Users/optas/Dropbox/matlab_projects/Image_Graphs/data/Centers_MSRC');
bow_dict = bow_dict.Centers;

bow_feats   = cell(1); 
for i=1:total_images    
    [F, ~]       = bag_of_words(images(:,:,:,i), bow_dict, grid_dim(1), grid_dim(2));  % Extract Bag Of Word Features
    bow_feats{i} = F;
end


%%
