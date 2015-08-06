%%  A Script researching Ensemble Procedures for Functional Maps.
    clr;
    gitdir;
    cd FmapLib/src
        
%% Load two Meshes and compute their LBs.
    meshfile       = '../data/input/tosca/dog1.off';
    mesh1          = Mesh(meshfile, 'dog1');
    mesh1.set_default_vertex_areas('barycentric');    
    LB1            = Laplace_Beltrami(mesh1);       
    feats1         = Mesh_Features(mesh1, LB1);
    
    meshfile       = '../data/input/tosca/dog2.off';
    mesh2          = Mesh(meshfile, 'dog2');
    mesh2.set_default_vertex_areas('barycentric');
    LB2            = Laplace_Beltrami(mesh2);
    feats2         = Mesh_Features(mesh2, LB2);

%%  Compute some Mesh features that will be used to produce Fmaps.
    neigs          = 300;       % Eigenvectors that will be used in computing wks/hks.
    wks_samples    = 100;
    hks_samples    = 100;    
    mc_samples     = 100; 
    gc_samples     = 100;
    feats1.compute_default_feautures(neigs, wks_samples, hks_samples, mc_samples, gc_samples);
	feats2.compute_default_feautures(neigs, wks_samples, hks_samples, mc_samples, gc_samples);
%%
%     save('../data/output/ensembles/demo_ensemble', 'mesh1', 'mesh2', 'LB1', 'LB2', 'feats1', 'feats2');              
%     load('../data/output/ensembles/demo_ensemble');    

%% Ensembling.    
    fmap        = Functional_Map(LB1, LB2);
    num_of_maps = 30;
    lambda      = 100;                           % Regulizer meta-parameter.
    ensembles   = cell(num_of_maps, 1);
    all_feats   = size(feats1.F, 2);
    min_feats   = ceil(0.1 * all_feats);
    
    for i=1:num_of_maps        
        feats_i     = floor(min_feats + rand() * (all_feats - min_feats));
        sample_feat = randsample(all_feats, feats_i)';
        source_f.F  = feats1.F(:, sample_feat);                             % TODO-P remove .F and make it a class.
        target_f.F  = feats2.F(:, sample_feat);       
        assert(size(source_f.F, 2) == size(target_f.F, 2) && size(target_f.F, 2) == feats_i);
        fmap.compute_f_map('frobenius_square', neigs, neigs, source_f, target_f, 'lambda', lambda);
  
        Fs = fmap.source_basis.project_functions(fmap.source_neigs, source_f.F);
        Fs = divide_columns(Fs, sqrt(sum(Fs.^2)));                                   % Normalize features.
        Ft = fmap.target_basis.project_functions(fmap.target_neigs, target_f.F);               
        Ft = divide_columns(Ft, sqrt(sum(Ft.^2)));                
              
        weight_i = sum(sum(abs((fmap.fmap * Fs) - Ft), 1)) ./ feats_i;

        ensembles{i}.map   = fmap.fmap;              
        ensembles{i}.error = weight_i;
        
    end
    
%% Aggregating Ensembles
%   1. Just the average of them.
    Xmean = zeros(size(ensembles{1}.map));
    for i=1:num_of_maps            
        Xmean = Xmean + ensembles{i}.map;
    end
    Xmean = (1/(num_of_maps)) .* Xmean;
    Fmean = Functional_Map(LB1, LB2);
    Fmean.set_fmap(Xmean);

%   2. Linear combination of fmaps based on weights.    
    Xlin = zeros(size(ensembles{1}.map));
    for i=1:num_of_maps            
        Xlin = Xlin + ((1/ensembles{i}.error) * ensembles{i}.map);
    end
    Flin = Functional_Map(LB1, LB2);
    Flin.set_fmap(Xlin ./ num_of_maps);    

%   3. Combine the learned errors to create a map that emphasized them.
    total_error =  0; 
    for i=1:num_of_maps            
        total_error = total_error + ensembles{i}.error;
        
    end
    Xlin2 = zeros(size(ensembles{1}.map));
    for i=1:num_of_maps            
        Xlin2 = Xlin2 + ((1 - (ensembles{i}.error / total_error)) * ensembles{i}.map);
    end
    Flin2 = Functional_Map(LB1, LB2);
    Flin2.set_fmap(Xlin2);


    
% Evaluating distortions.
    % Load symmetries.    
    fid = fopen('../data/input/tosca_symmetries/dog.sym'); % TODO-P add to IO.read_symmetries();
    C   = textscan(fid, '%s', 'delimiter', '\n');          % Read symmetries:
    fclose(fid);
    symmetries   = str2double(C{:});    
    
    % Make groundtuth correspondence (for tosca, within same class i node matched i).
    groundtruth = (1:mesh1.num_vertices)';      
    
    [dists_m, indices]       = Fmean.pairwise_distortion(groundtruth, 'nsamples', 400, 'symmetries', symmetries);            
    [dists_l, ~]             = Flin.pairwise_distortion(groundtruth, 'indices', indices, 'symmetries', symmetries);        
    [dists_l2, ~]            = Flin2.pairwise_distortion(groundtruth, 'indices', indices, 'symmetries', symmetries);        
    
    % Use All probe functions.
    fmap.compute_f_map('frobenius_square', neigs, neigs, feats1, feats2, 'lambda', lambda);                      
    [dists_a, ~] = fmap.pairwise_distortion(groundtruth, 'indices', indices, 'symmetries', symmetries);
            
    % Use the optimal/groundtruth map.
    X_opt               = Functional_Map.groundtruth_functional_map(LB1.evecs(neigs), LB2.evecs(neigs), groundtruth, diag(LB2.A));       
    [dists_o,  ~]       = Functional_Map.pairwise_distortion_of_map(X_opt, LB1, LB2, groundtruth, 'indices', indices, 'symmetries', symmetries);
    
    % Use only the wks/hks features.
    feats1_intr.F = feats1.F(:, wks_samples + hks_samples);                                                   % TODO-P utilize feats.index.
    feats2_intr.F = feats2.F(:, wks_samples + hks_samples);
    fmap.compute_f_map('frobenius_square', neigs, neigs, feats1_intr, feats2_intr, 'lambda', lambda);            % Use hks/wks only.
    [dists_i, ~] = fmap.pairwise_distortion(groundtruth, 'indices', indices, 'symmetries', symmetries);
    
    mean(dists_o)
    mean(dists_l)    
    mean(dists_i)
    mean(dists_a)
    mean(dists_m)
    
    

          
%     for i=1:num_of_maps        
%         [dists, ~]       = ensembles{i}.map.pairwise_distortion(groundtruth, 'indices', indices, 'symmetries', symmetries);        
%         [dists2, ~]      = Fmean.pairwise_distortion(groundtruth, 'indices', indices, 'symmetries', symmetries);        
%         all_dists{i}.ens = dists;
%         all_dists{i}.agg = dists2;
%     
%     end        
%     for i=1:num_of_maps       %TODO-P save sampling numbers
%         fprintf('Ensmble with %d feautures. Mean: %f STD:%f \n', ensembles{i}.nfeat, mean(all_dists{i}.ens), std(all_dists{i}.ens))
%         fprintf('A_Fmap Mean:%f STD:%f \n', mean(all_dists{i}.agg), std(all_dists{i}.agg))
%         ensembles{i}.nfeat
%     end

   
    
