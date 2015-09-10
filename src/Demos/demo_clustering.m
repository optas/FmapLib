clr;
gitdir;
cd 'FmapLib/src';

%% Load the collection of Meshes and their semantic attributes (i.e., class of each represented mesh).
collection_name = 'Tosca';
collection_dir  = '../data/input/Tosca';
semantics_file  = '../data/input/TOSCA_class_attributes';
Tosca           = Mesh_Collection(collection_name, collection_dir, semantics_file);
Tosca.size()


%% Compute Laplacian Basis and Mesh Features.
neigs     = 64; 
area_type = 'barycentric';
for mesh_name = Tosca.meshes.keys                  % TODO-E: implement class iterator i.e., for m = Tosca {do }.
    mesh_name = mesh_name {:};
    m        = Tosca.meshes(mesh_name);    
    m.set_default_vertex_areas('barycentric');
    if sum(m.barycentric_v_area <= 0) ~= 0
        mesh_name
    end
end
                
Tosca.compute_laplace_beltrami_basis(neigs, area_type);

%% Compute Features.
hks_samples = 100; wks_samples = 100; mc_samples = 50; gc_samples = 50;
Tosca.compute_default_feautures(hks_samples, wks_samples, mc_samples, gc_samples);  % Computes hks, wks, mean/gauss curvature for each mesh.

%% Save - Load
% save('../data/output/mike_vica_collection', 'Tosca', '-v7.3')
% load('../data/output/tosca_collection');              
% loas('../data/output/mike_vica_collection');

%% Split dataset into training/testing samples according to their classes.
cats   = Tosca.meshes_with_semantic_condition('Class_2', 'Cat');                % See semantics file.
dogs   = Tosca.meshes_with_semantic_condition('Class_2', 'Dog');
humans = Tosca.meshes_with_semantic_condition('Class_1', 'Human');
mikes  = Tosca.meshes_with_semantic_condition('Class_3', 'Michael');
vicas  = Tosca.meshes_with_semantic_condition('Class_2', 'Female');

%%
% [train_data, test_data, train_labels, test_labels] = Learning.sample_class_observations(0.2, cats, dogs, humans);
[train_data, test_data, train_labels, test_labels] = Learning.sample_class_observations(0.2, mikes, vicas);

%% Compute Functional Maps
% All pairs between training and testing are used (expensive - improve by parfor + sampling).
pairs = cell(2* length(test_data) * length(train_data), 2);
p = 1;
for i=test_data
    for j=train_data
        pairs{p,1}   = i{:}; pairs{p,2}   = j{:};
        pairs{p+1,1} = j{:}; pairs{p+1,2} = i{:};
        p = p + 2;
    end
end
all_maps = Tosca.compute_fmaps(pairs, Tosca.raw_features, 'frobenius_square', 'lambda', 20);

%% 1. Naive Learning
scores = Learning.fmap_classify_naively(test_data, train_data, train_labels, all_maps);
%%

class_id      = 1;
inter_pairs   = Learning.inter_class_pairs(train_labels, train_data, class_id);
init_maps     = Tosca.compute_fmaps(inter_pairs, Tosca.raw_features, 'frobenius_square', 'lambda', 20);


%%
weights1      = Learning.feature_weights_of_class(init_maps);
%%
low_rank_maps = Learning.low_rank_filtering_of_fmaps(init_maps);
%%
weights2      = Learning.feature_weights_of_class(low_rank_maps);
%%
scores2       = Learning.fmap_classify_with_trained_weights(test_data, train_data, train_labels, weights);