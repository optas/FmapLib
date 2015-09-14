clear; clc;
[dp, cp] = get_project_paths('ImageJointUnderstanding');
            
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
  
    gist_signatures{i} = load([image_folder image_name '_gist.mat']);    
    gist_signatures{i}  = gist_signatures{i}.gist_feat;
end

num_images = length(images);

%% Resize every image and extract a laplacian basis.
new_height = 10;
new_width  = 10;

radius     = 1;
eigs_num   = 2;

sigma_f   = 800 / (255^2);   

image_laplacians = cell(length(all_image_files), 1);

for i= 1:num_images
    images{i}.set_resized_image(new_height, new_width);    
    im_i = images{i}.get_resized_image();    
    new_width = im_i.width;
        
    G = Image_Graph(im_i, 'r_radius_connected',  radius);
    fprintf('Graph %d constructed.\n', i)

    
    Fi = im_i.color();
    sigma_s  =  2 * (0.1 * norm([new_width-1, new_height-1]))^2;    
    G.adjust_weights_via_feature_differences(Fi , 'normalized_cut', 'sigma_s', sigma_s, 'sigma_f', sigma_f);    
    
    image_laplacians{i} = Laplacian(G.Gw, 'norm');
    image_laplacians{i}.get_spectra(eigs_num);
    fprintf('Laplacian %d constructed.\n', i)
end
% save([dp 'output/laplacians_aeroplane_left_norm'], 'image_laplacians');


%% Extract hog feature (pixel-wise) and project them into Laplacian basis.
hog_feats = cell(num_images ,1);
proj_hogs = cell(num_images, 1);

for i= 1:num_images
    im_i = images{i}.get_resized_image();    
    
    h = im_i.height;
    w = im_i.width;
    
    hog_feats{i} = Image_Features.hog_signature(im_i);    
    
    Fp = reshape(hog_feats{i}, h*w, size(hog_feats{i}, 3));       % Now Fp, contains the pixel-level features stack in a long vector column-wise.
    Fp = divide_columns(Fp, sqrt(sum(Fp.^2)));                    % Rescale to unit norm.    
       
    proj_hogs{i}    = image_laplacians{i}.project_functions(eigs_num, Fp);
end
% save

%% Make all F-maps.
all_fmaps = cell(num_images, num_images);
regulizer_w = 15;
for i = 1:num_images
    evals_i = image_laplacians{i}.evals(eigs_num);
    for j = 1:num_images
        if i ~= j
            evals_j = image_laplacians{j}.evals(eigs_num);            
            all_fmaps{i,j} = Functional_Map.l1_and_frobenius_norms_cvx(proj_hogs{i}, proj_hogs{j}, evals_i, evals_j, regulizer_w);
        end        
    end
end

%% k-nn GIST initial image network
k_neighbors = 5;
gist_desc = cell2mat(gist_signatures);
[nns, gist_dists] = knnsearch(gist_desc, gist_desc, 'k', k_neighbors+1);
nns = nns(:, 2:k_neighbors+1);
gist_dists = gist_dists(:, 2:k_neighbors+1);


%% Extract and project features for every patch
tic
patch_feat = cell(num_images, 1);
for i = 1:num_images    
    o_i  = obj_proposals{i};
    patch_feat = cell(length(o_i), 1);
    
    im_i = images{i}.get_resized_image();        
    new_height = im_i.height;
    new_width  = im_i.width;
    
    for pi = 1:length(o_i)                   % Proposals for image_i
        if ~ Patch.are_valid_corners(o_i(pi, :),  images{i}) 
            continue;
        end
        corners_pi = Patch.find_new_corners(images{i}.height, images{i}.width, new_height, new_width, o_i(pi,:));
        
        Fs = Patch.extract_patch_features(corners_pi , hog_feats{i});        
        Fs = reshape(Fs, new_height*new_width, size(Fs, 3));        
        Fs = divide_columns(Fs, sqrt(sum(Fs.^2)));                    % Rescale to unit norm.    
        Fs = image_laplacians{i}.project_functions(eigs_num, Fs);
        patch_feat{i, pi} = Fs;
    end
end
toc

%% Rank every proposal (no triplets involved)

top_p   = 5;                                 % How many top-scoring patches to keep per image at each iteration
top_patches = zeros(num_images, top_p);

for i = 1:num_images
    im_i = images{i}.get_resized_image();    

    orig_height = images{i}.height;
    orig_width = images{i}.width;

    new_height = im_i.height;
    new_width  = im_i.width;

    o_i  = obj_proposals{i};

    score_i = zeros(length(o_i) , 1);
    
    for pi = 1:length(o_i)                   % Proposals for image_i
        if ~ Patch.are_valid_corners(o_i(pi, :),  images{i}) 
%                 o_i(pi, :)
%                 input('')
%                 break;
            continue;

        end
        corners_pi = Patch.find_new_corners(orig_height, orig_width, new_height, new_width, o_i(pi,:));
        
        


        for j = nns(i, :)
            misalignment = 0;             
            im_j = images{j}.get_resized_image();  
            orig_heightj = images{j}.height;
            orig_widthj = images{j}.width;
            new_heightj = im_j.height;
            new_widthj = im_j.width;

            o_j = obj_proposals{j};

            for pj = 1:length(o_j)      % Proposals for image_j                        
                if ~ Patch.are_valid_corners(o_j(pj, :),  images{j})                        
                    continue; 

                end
                corners_pj = Patch.find_new_corners(orig_heightj, orig_widthj, new_heightj, new_widthj, o_j(pj,:));
                Ft = Patch.extract_patch_features(corners_pj , hog_feats{j});
                Ft = reshape(Ft, new_heightj*new_widthj, size(Ft, 3));
                Ft = image_laplacians{i}.project_functions(eigs_num, Ft);

                misalignment = misalignment + sum(sum(abs((all_fmaps{i,j} * Fs) - Ft), 1)); % alignment error (L1 dist) of probe functions.       
                misalignment = misalignment + sum(sum(abs((all_fmaps{j,i} * Ft) - Fs), 1));
            end

            score_i(pi) = score_i(pi) + (misalignment / length(o_j));
        end        
    end
    [sort_scores, indices] = sort(score_i);      
    top_patches(i, :) = indices(1:top_p)            
end
%%
% Zimo corloc
% Zimo corloc with standout

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


%% Visualize an eigenvector coming from the mincut method
E = image_laplacians{3}.evecs(10);   
E = E(:,2);
E = reshape(E, new_height, new_width);
E = (E - min(E(:))) ./ (max(E(:)) - min(E(:)));    % imshow expects values in [0,1] if they are floats.
imshow(E)




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
    accu = mean(sqrt(sum(abs(Fp - reconstructed).^2)))
end


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



%% Extract BOW-feature.
bow_dict = load('/Users/optas/Dropbox/matlab_projects/Image_Graphs/data/Centers_MSRC');
bow_dict = bow_dict.Centers;

bow_feats   = cell(1); 
for i=1:total_images    
    [F, ~]       = bag_of_words(images(:,:,:,i), bow_dict, grid_dim(1), grid_dim(2));  % Extract Bag Of Word Features
    bow_feats{i} = F;
end


%%
