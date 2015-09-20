clear; clc;
[dp, cp] = get_project_paths('ImageJointUnderstanding');
            
%% Load VOC image collection with gist, gt_segmentation and object_proposals.
image_folder        = [dp 'VOC/2007/Train_Validate_Data/VOC2007_6x2/aeroplane_left/'];
all_image_files     = rdir([image_folder, [filesep '**' filesep]], 'regexp(name, ''\.jpg$'')');
images              = cell(length(all_image_files), 1);
obj_proposals       = cell(length(all_image_files), 1);
gist_signatures     = cell(length(all_image_files), 1);
orig_hog_patch_sigs = cell(length(all_image_files), 1);
for i=1:length(images)
    full_path         = all_image_files(i).name;
    path_substrings   = strsplit(full_path, filesep); 
    last_word         = path_substrings{end};
    image_name        = last_word(1:end-4);                             % Relying on the fact that length('.jpg') == length('.obj') == 4.                                                               
    images{i}         = Image(full_path, image_name);    
    
    gt_segmentation   = load([image_folder image_name '.mat']);
    gt_segmentation   = gt_segmentation.bbox_list;
    images{i}.set_gt_segmentation(gt_segmentation);
    
    obj_proposals{i}       = load([image_folder image_name '_seg.mat']);    
    orig_hog_patch_sigs{i} = obj_proposals{i}.feat.hist;
    obj_proposals{i}       = int16(obj_proposals{i}.feat.boxes);
    
    
    gist_signatures{i}  = load([image_folder image_name '_gist.mat']);    
    gist_signatures{i}  = gist_signatures{i}.gist_feat;
end
num_images = length(images);

%% Test 40, 43 images with iterative fmaps
%% Exctract image laplacians
pair = [40, 43];
new_height = 128;
new_width  = NaN;
radius     = 3;
eigs_num   = 64;
sigma_f   = 800 / (255^2);   
image_laplacians = cell(length(all_image_files), 1);
for i = pair
    images{i}.set_resized_image(new_height, new_width);        
    im_i = images{i}.get_resized_image();    
    new_width = im_i.width;
        
    G = Image_Graph(im_i, 'r_radius_connected',  radius);
    fprintf('Graph %d constructed.\n', i)

    Fi = im_i.color();
    sigma_s  =  2 * (0.1 * norm([new_width-1, new_height-1]))^2;    
    G.adjust_weights_via_feature_differences(Fi , 'normalized_cut', 'sigma_s', sigma_s, 'sigma_f', sigma_f);     % make it a seperate Graph    
    image_laplacians{i} = Laplacian(G.Gw, 'norm');
    image_laplacians{i}.get_spectra(eigs_num);
    fprintf('Laplacian %d constructed.\n', i);    
end

%% Visualize Small Eigenvectors
im_id = 43; eig_to_vis = 20;
E = image_laplacians{im_id}.evecs(eigs_num);   
new_width = images{im_id}.get_resized_image().width;
figure();
% set(f, 'DefaultFigureWindowStyle','docked')
for i = 1:eig_to_vis     
    eig_i = E(:,i);
    eig_i = reshape(eig_i, new_height, new_width);
    eig_i = (eig_i - min(eig_i(:))) ./ (max(eig_i(:)) - min(eig_i(:)));    % imshow expects values in [0,1] if they are floats.
    subplot(5,4,i); imshow(eig_i);    
end
suptitle('20 Smallest Eigenvectors (Min-Cut method).')    

%% Extract hog features over the Whole image for every pixel and project them into Laplacian basis.
hog_feats = cell(num_images ,1);
proj_hogs = cell(num_images, 1);
for i=pair
    im_i = images{i}.get_resized_image();        
    h = im_i.height;  w = im_i.width;    
    hog_feats{i} = Image_Features.hog_signature(im_i);        
    Fp = reshape(hog_feats{i}, h*w, size(hog_feats{i}, 3));       % Now Fp, contains the pixel-level features stack in a long vector column-wise.
    Fp = divide_columns(Fp, sqrt(sum(Fp.^2)));                    % Rescale to unit norm.           
    proj_hogs{i}    = image_laplacians{i}.project_functions(eigs_num, Fp);
end
%% Extract and project features for every patch
method = 'zero_pad' ;  
% method = 'tiling';
patch_feat = cell(num_images, 1);
for i = pair
    o_i  = obj_proposals{i};
    patch_feat{i} = cell(length(o_i), 1);    
    im_i = images{i}.get_resized_image();        
    new_height = im_i.height;
    new_width  = im_i.width;
    
    for pi = 1:length(o_i)                          % Proposals for image_i
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

%% Original Hog of patches, build matching.
raw_hog_matching = cell(num_images, num_images, 2);
for i = 1:length(pair)        
    for j = i+1:length(pair)
        I=pair(i);
        J=pair(j);
        [raw_hog_matching{I,J, 1}, raw_hog_matching{I,J,2}] = knnsearch(orig_hog_patch_sigs{J}, orig_hog_patch_sigs{I}, 'k', 2);
    end
end

%% Use original HOG matching for aligning probe functions (pathces)
probe_source = cell2mat(patch_feat(I,:));
probe_target = cell(1);
for i = 1:length(raw_hog_matching{I,J,1})
    match = raw_hog_matching{I,J,1}(i);
    probe_target{i} = patch_feat{J, match};
end
probe_target = cell2mat(probe_target);

%% Make all F-maps.
all_fmaps = cell(num_images, num_images);
regulizer_w = 20;
weight_mask = size(patch_feat{I,1},2);

for i = 1:length(pair)
    evals_i  = image_laplacians{pair(i)}.evals(eigs_num);
    for j = i+1:length(pair)
        if pair(i) ~= pair(j)
            evals_j = image_laplacians{pair(j)}.evals(eigs_num);            
            W = double(1) ./ raw_hog_matching{pair(i), pair(j), 2}(:,1);
            W = repmat(W', weight_mask, 1);
            W = W(:);            
            [all_fmaps{pair(i), pair(j)}.X, all_fmaps{pair(i), pair(j)}.W, ran_it] = Functional_Map.iteratively_refined_fmap(probe_source, probe_target, ...
                                                                             evals_i, evals_j, regulizer_w, 'max_iter', 10, 'weight_mask', weight_mask, 'weights', W);
        end        
    end
end
%% Get Fmap based - matching
X_final = all_fmaps{I,J}.X(:, :, end);
w_final = all_fmaps{I,J}.W(:, end);
w_final = w_final(1:weight_mask:end);
nns     = raw_hog_matching{I,J,1}(:,1);
dist    = w_final;
[sorted_dis, sorted_ind] = sort(dist, 'descend');
sorted_nns = nns(sorted_ind);                   % Sort the top-scoring matches.

%% Plot the Matchings
ps_name = 'Fmap_based_matching_zero_pad_hog_init_fro_10_iter';
for c = 1:length(sorted_dis)
    c
    if c == 1        
        f = figure('visible', 'off');        
        subplot(2,1,1); imshow(images{I}.CData);
        subplot(2,1,2); imshow(images{J}.CData);        
        suptitle('Original Images');
        print ( '-dpsc2', ps_name, f);        
        
        f = figure('visible', 'off');     % Plot weight evolution (over diff. iterations)
        plot(all_fmaps{I, J}.W(:,1)); hold on; 
        plot(all_fmaps{I, J}.W(:,2), 'g'); hold on;  
        plot(all_fmaps{I, J}.W(:,end), 'r')
        ylabel('Weight'); xlabel('Patches: each contributes 36 functions with same weight (plateau).')
        legend('init', '2-iter', sprintf('%d-iter(last)',ran_it))
        print ( '-dpsc2', ps_name, '-append', f);                
    end
    
    f = figure('visible', 'off');
    
    op_I = obj_proposals{I}(sorted_ind(c), :);
    subplot(2,1,1); imshow(images{I}.content_in_rectangle(op_I));    
    xlabel(sprintf('Cor-loc: %f', overlap_with_mask(op_I, bitmask_I)));
  
    
    op_J = obj_proposals{J}(sorted_nns(c), :);
    subplot(2,1,2); imshow(images{J}.content_in_rectangle(op_J));    
    xlabel(sprintf('Cor-loc: %f', overlap_with_mask(op_J, bitmask_J)));
        
        
    suptitle(sprintf('Fmap-based similarity %f .', sorted_dis(c)));    
    print ( '-dpsc2', ps_name, '-append', f )    
end
close all
%%


%% Plot pairs of matching (decreasing order of closeness) based on HOG
nns  = raw_hog_matching{I,J,1}(:,1);
dist = raw_hog_matching{I,J,2}(:,1);
[sorted_dis, sorted_ind] = sort(dist);
sorted_nns = nns(sorted_ind);                   % Sort the top-scoring correspondences.


%% Make bitmaps for GT
bitmask_I = zeros(images{I}.height, images{I}.width);
for gt_I = images{I}.gt_segmentation            
    if isempty(gt_I)
            continue
    end        
    gt = gt_I{:};                
    bitmask_I(gt(2):gt(4), gt(1):gt(3)) = 1;      
end

bitmask_J = zeros(images{J}.height, images{J}.width);
for gt_J = images{J}.gt_segmentation            
    if isempty(gt_J)
            continue
    end        
    gt = gt_J{:};                
    bitmask_J(gt(2):gt(4), gt(1):gt(3)) = 1;      
end            










