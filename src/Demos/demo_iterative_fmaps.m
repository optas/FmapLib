clear; clc;
[dp, cp]    = get_project_paths('ImageJointUnderstanding');
            
%% Load a pair of images that has similar foreground and drastically different background.
object_type         = 'winston_chair';
image_folder        = [dp 'Our_Constructions/Distinct_Pairs/' object_type '/'];
all_image_files     = rdir([image_folder, [filesep '**' filesep]], 'regexp(name, ''\.jpg$'')');
images              = cell(length(all_image_files), 1);
for i=1:length(images)
    full_path       = all_image_files(i).name;
    path_substrings = strsplit(full_path, filesep); 
    last_word       = path_substrings{end};
    image_name      = last_word(1:end-4);               % Relying on the fact that length('.jpg') == length('.obj') == 4.                                                               
    images{i}       = Image(full_path, image_name);    
    
%     gt              = [image_folder 'GroundTruth/' image_name '.png'];
%     images{i}.set_gt_segmentation(imread(gt));    
    
%     gt_segmentation   = load([image_folder image_name '.mat']);
%     gt_segmentation   = gt_segmentation.bbox_list;
%     images{i}.set_gt_segmentation(gt_segmentation);
       
end

num_images = length(images);

%% Load/compute image proposals
approx_prop_per_image = 100;
obj_proposals   = cell(length(all_image_files), 1);
params = load([cp 'Randomized_Prim_Object_Proposal/config/rp.mat']);  % Default Params for RP method.
params = params.params;
params.approxFinalNBoxes = approx_prop_per_image;    

for i=1:length(images)    
    obj_proposals{i} = RP(images{i}.CData, params);
end

%%
% 
% op_I = obj_proposals{I}
%     subplot(2,1,1);  images{I}.plot_patch(op_I);

images{1}.plot
plotBoxes(obj_proposals{1}, 'random', [], '-')
%%







%% Filterout bad proposals (too small or even not rectangular)
obj_filtered = cell(num_images, 1);
for i=1:length(images)    
    keep = [];
    for j=1:size(obj_proposals{i}, 1)
        if Patch.is_valid(obj_proposals{i}(j,:), images{i}) && Patch.is_within_limits(obj_proposals{i}(j,:), 0.9, 0.9, images{i}) 
            keep(end+1) = j;
        end
    end    
    obj_filtered{i} = obj_proposals{i}(keep,:);    
end
obj_proposals = obj_filtered;

%% Compute patch-wide HOG Feats
% patch_wide_hog = cell(num_images, 1);
% for i = 1:num_images
%     patch_wide_hog{i} = get_hog(obj_proposals{i}, images{i}.CData);
% end 
% %%
patch_wide_hog2 = cell(num_images, 1);
for i = 1:num_images
    seg.coords = obj_proposals{i};    
    patch_wide_hog2{i} = extract_segfeat_hog(images{i}.CData, seg);
    patch_wide_hog2{i} = patch_wide_hog2{i}.hist;
end 
patch_wide_hog = patch_wide_hog2;
%% Specifying a single pair of images;
I = 1;
J = 2;
pair = [I, J];
num_images = 2;
%% Exctract image laplacians
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
save([image_folder 'mincut_lapacians_128_NaN_64eigs'], 'image_laplacians')
%% or
new_height = 128;
new_width  = NaN;
eigs_num = 64;
for i = pair
    images{i}.set_resized_image(new_height, new_width);        
end
load([image_folder 'mincut_lapacians_128_NaN_64eigs'])

%% Visualize Small Eigenvectors of single image.
im_id = I; eig_to_vis = 20;
E = image_laplacians{I}.evecs(eigs_num);   
new_width = images{I}.get_resized_image().width;
figure();
for i = 1:eig_to_vis     
    eig_i = E(:,i);
    eig_i = reshape(eig_i, new_height, new_width);
    eig_i = (eig_i - min(eig_i(:))) ./ (max(eig_i(:)) - min(eig_i(:)));    % imshow expects values in [0,1] if they are floats.
    subplot(5,4,i); imshow(eig_i);    
end
suptitle('20 Smallest Eigenvectors (Min-Cut method).')    

%% Extract hog features over every pixel of the resized image and project them into Laplacian basis.
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


%% Black them mother fucker. 
hog_feats = cell(num_images ,1);
for i=pair
    im_i = images{i}.get_resized_image();        
    h = im_i.height;  w = im_i.width;    
    hog_feats{i} = ones(h,w);   
end

%% Extract and project features for every patch
method = 'zero_pad' ;  
% method = 'tiling';
patch_feat        = cell(num_images, 1);
for i = pair
    o_i           = obj_proposals{i};
    patch_feat{i} = cell(length(o_i), 1);
    im_i          = images{i}.get_resized_image();        
    new_height    = im_i.height;
    new_width     = im_i.width;
    
    for pi = 1:length(o_i)                          % Proposals for image_i        
        corners_pi = Patch.find_new_corners(images{i}.height, images{i}.width, new_height, new_width, o_i(pi,:));
        
        
        Fs = Patch.extract_patch_features(corners_pi , hog_feats{i}, method);                        
        
        
        
        Fs = reshape(Fs, new_height*new_width, size(Fs, 3));                              
        Fs = divide_columns(Fs, max(sqrt(sum(Fs.^2)), 1e-20));                    % Rescale to unit norm.   
        Fs = image_laplacians{i}.project_functions(eigs_num, Fs);
        patch_feat{i, pi} = Fs;
    end
end
feat_num = size(patch_feat{1,1}, 2);

%% Build hog based matchings and align the Probe Functions with it.
raw_hog_matching = cell(num_images, num_images, 2);
aligned_probes   = cell(num_images, num_images, 2);
for i = 1:length(pair)        
    for j = 1:length(pair)                      % NN is not symmetric thus we need both (i,j) and (j,i)
        if i ~= j
            pi = pair(i);
            pj = pair(j);
            [raw_hog_matching{pi, pj, 1}, raw_hog_matching{pi, pj, 2}] = knnsearch(patch_wide_hog{pj}, patch_wide_hog{pi}, 'k', 2);
                         
%             % Align probe functions
%             probe_source = cell2mat(patch_feat(pi,:));            
%             probe_target = cell(1);
%             for k = 1:length(raw_hog_matching{pi, pj, 1})
%                 match = raw_hog_matching{pi, pj, 1}(k);
%                 probe_target{k} = patch_feat{pj, match};
%             end
%             probe_target = cell2mat(probe_target);            
%             aligned_probes{pi, pj, 1} = probe_source;
%             aligned_probes{pi, pj, 2} = probe_target;            
        end
    end
end

%% Make all F-maps.
all_fmaps = cell(num_images, num_images);
regularizer_w = 20;
weight_mask = feat_num;

for i = 1:length(pair)
    evals_i  = image_laplacians{pair(i)}.evals(eigs_num);
    evals_i  = evals_i ./max(evals_i);
    
    for j = 1:length(pair)
        if pair(i) ~= pair(j)
            evals_j = image_laplacians{pair(j)}.evals(eigs_num);            
            evals_j = evals_j ./ max(evals_j);
            
            W = double(1) ./ raw_hog_matching{pair(i), pair(j), 2}(:,1);
            W = repmat(W', weight_mask, 1);
            W = W(:);        
            
            probe_source = aligned_probes{pair(i), pair(j), 1};
            probe_target = aligned_probes{pair(i), pair(j), 2};
            
            [all_fmaps{pair(i), pair(j)}.X, all_fmaps{pair(i), pair(j)}.W] = Functional_Map.iteratively_refined_fmap(probe_source, probe_target, ...
                                                                             evals_i, evals_j, regularizer_w, 'max_iter', 200, 'weight_mask', weight_mask) ; % , 'weights', W);
        end        
    end
end

ran_for = size(all_fmaps{I,J}.X, 3);
%%
ran_for  = size(all_fmaps{I,J}.X, 3);
diff = zeros(ran_for, 1);
for i=1:ran_for
    diff(i) = max(abs(all_fmaps{I,J}.W(:,i+1)  - all_fmaps{I,J}.W(:,i) ));
end
plot([1:ran_for], diff);

%% Get Fmap based - matching
I = 1;
J = 2;
X_final = all_fmaps{I,J}.X(:, :, end);
w_final = all_fmaps{I,J}.W(:, end);
w_final = w_final(1:weight_mask:end);
nns     = raw_hog_matching{I, J, 1}(:,1);
dist    = w_final;
[sorted_dis, sorted_ind] = sort(dist, 'descend');
sorted_nns = nns(sorted_ind);                   % Sort the top-scoring matches.
method = sprintf('Fmap_lambda_%d_uni_weights', regularizer_w) ;

%% Make bitmaps for GT

bitmask_I = images{I}.gt_segmentation;
bitmask_I(bitmask_I ~= 0) = 1;

bitmask_J = images{J}.gt_segmentation;
bitmask_J(bitmask_J ~= 0) = 1;
%  
% bitmask_I = zeros(images{I}.height, images{I}.width);
% for gt_I = images{I}.gt_segmentation            
%     if isempty(gt_I)
%             continue
%     end        
%     gt = gt_I{:};                
%     bitmask_I(gt(2):gt(4), gt(1):gt(3)) = 1;      
% end
% 
% bitmask_J = zeros(images{J}.height, images{J}.width);
% for gt_J = images{J}.gt_segmentation            
%     if isempty(gt_J)
%             continue
%     end        
%     gt = gt_J{:};                
%     bitmask_J(gt(2):gt(4), gt(1):gt(3)) = 1;      
% end            

%% If you want to plot pairs of matching (decreasing order of closeness) based on HOG
nns  = raw_hog_matching{I,J,1}(:,1);
dist = raw_hog_matching{I,J,2}(:,1);
[sorted_dis, sorted_ind] = sort(dist);
sorted_nns = nns(sorted_ind);                   % Sort the top-scoring correspondences.
method = 'HOG';

%% Plot the Matchings
clc
ps_name = [method '_based_match_' object_type '_' sprintf('%d_%d', I, J)]        

top_pairs = length(obj_proposals{1});
for c = 1:top_pairs
    c
    if c == 1        
        f = figure('visible', 'off');        
        subplot(2,1,1); imshow(images{I}.CData);
        subplot(2,1,2); imshow(images{J}.CData);        
        suptitle('Original Images');
        print ( '-dpsc2', ps_name, f);        
        
        if strcmp(method, 'Fmap')            
            f = figure('visible', 'off');     % Plot weight evolution (over diff. iterations)
            plot(all_fmaps{I, J}.W(:,1)); hold on; 
            plot(all_fmaps{I, J}.W(:,end-1), 'g'); hold on;  
            plot(all_fmaps{I, J}.W(:,end), 'r')
            ylabel('Weight'); xlabel('Patches: each contributes 36 functions with same weight (plateau).')
            legend('init', sprintf('%d-iter', ran_for-1), sprintf('%d-iter(last)', ran_for));
            print ( '-dpsc2', ps_name, '-append', f);                
        end
    end
    
    f = figure('visible', 'off');
    
    op_I = obj_proposals{I}(sorted_ind(c), :);    
    subplot(2,1,1);  images{I}.plot_patch(op_I);
%     xlabel(sprintf('Cor-loc: %f', overlap_with_mask(op_I, bitmask_I)));
      
    op_J = obj_proposals{J}(sorted_nns(c), :);
    subplot(2,1,2); images{J}.plot_patch(op_J);
%     xlabel(sprintf('Cor-loc: %f', overlap_with_mask(op_J, bitmask_J)));
        
    
    suptitle(sprintf([method ' based similarity %f .'], sorted_dis(c)));    
    print ( '-dpsc2', ps_name, '-append', f )    
end
close all
%%











