%% Load the collection of Meshes and their semantic attributes (i.e., class of each represented mesh).
gitdir
cd 'FmapLib/src'

collection_name = 'Tosca';
collection_file = '../data/input/tosca';
semantics       = '../data/input/TOSCA_class_attributes';
Tosca           = Mesh_Collection(collection_name, collection_file, semantics);

%% Compute Laplacian Basis and Mesh Features
neigs     = 128; 
area_type = 'barycentric';
Tosca.compute_laplace_beltrami_basis(neigs, area_type);
Tosca.compute_default_feautures();                      % Computes hks, wks, mean/gauss curvature for each mesh.

%% Save - Load
% save('../data/output/tosca_collection', 'Tosca', '-v7.3')
% load('../data/output/tosca_collection');              


%% Split dataset into training/testing samples according to their classes.
cats   = Tosca.meshes_with_semantic_condition('Class_2', 'Cat');    % See semantics file.
dogs   = Tosca.meshes_with_semantic_condition('Class_2', 'Dog');
humans = Tosca.meshes_with_semantic_condition('Class_1', 'Human');
[train_data, test_data, train_labels, test_labels] = Learning.sample_class_observations(0.2, cats, dogs, humans);


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



%% Learning
Learning.fmap_classify_naively(test_data, train_data, train_labels, all_maps);






%%
offs = [];
mi   = [];
for i=1:2:length(all_maps)
    A = all_maps{i}.fmap * all_maps{i+1}.fmap;    
    offs(end+1) = sum(sum(abs(A))) - sum(abs(diag(A)));
    mi(end+1)   = norm(A-eye(size(A,1)));
end

