% Load IC EEG-fMRI networks and vizualize them with brain net viewer
%
% 15-04-2016 Jonathan Wirsich CEMEREM
% 01-02-2017 Adapt to Zika Data Beckman
%
% Version adapted for
% Jonathan Wirsich, Enrico Amico 2020
%
% Wirsich, J., Amico, E., Giraud A.L. Goñi, J, Sadaghiani S.,2020 
% Multi-timescale hybrid components of the functional brain connectome: A bimodal EEG-fMRI decomposition
% Network Neuroscience (2020) 4 (3): 658–677. https://doi.org/10.1162/netn_a_00135

%% includes
%add brain connectivity toolbox (download from https://sites.google.com/site/bctnet/)
% addpath('TODO/lib/2019_03_03_BCT/')
%add brainnetviewer (download from https://www.nitrc.org/projects/bnv/)
%tested version: 2015-08-07
brainnet_dir = 'TODO/lib/BrainNetViewer_20150807';
addpath(brainnet_dir)

%% Configuration
%load corrdinates
reflect_folder = fileparts(mfilename('fullpath'));

coords = importdata(fullfile(reflect_folder, '/../data/destr_coords_simple_nosubc.txt'));

%% local objects
destrLabels=importdata(fullfile(reflect_folder, '/../data/destr_label_nosubc.txt'));

%load yeoorder
load(fullfile(reflect_folder, '../data/aparc_a2009_yeoRS7_148reg_eeg_nosubc.mat'))
coords = coords(yeoOrder, :);
destrLabels = destrLabels(yeoOrder);

path_local = fullfile(reflect_folder,'/../data/resu/PC75IC10_matrix_perc/');

%frequency Paris
labels{1} = 'Frequency_Main_PC75IC10';
%ICN Paris
labels{2} = 'ICN_Main_PC75IC10';
%frequency Marseille
labels{3} = 'Frequency_Generalization_PC75IC10';
%ICN Marseille
labels{4} = 'ICN_Generalization_PC75IC10';

%% main
regions = 148;

for label_it = 1:length(labels)
    %load nbs adja-matrix
    filename = ['99_percentile_' labels{label_it} '_mask_T_eeg'];

    tmp = importdata([path_local filename '.mat']);

    adja = tmp;
    % adja = tmp.mask_T_eeg;
    % adja = tmp.mask_T_eeg;
    %take only largest nbs network
    adja(adja>1) = 0;
    %get nodal degrees
    degr = degrees_und(adja);

    %load labels
    %build node file: x,y,z,#subnet,node-weight,label 
    %see brainnetvioewer documentation for format description
    nodes = cell(length(degr(degr>0)),6);
    count = 1;
    for r1 = 1:regions
        %TODO squeeze
        if degr(r1) > 0

    %         display([l1 '-' int2str(degr(r1))])
    %         listdegr{end+1} = [l1 '-' int2str(degr(r1))];
            %coords
            nodes{count, 1} = coords(r1,1);
            nodes{count, 2} = coords(r1,2);
            nodes{count, 3} = coords(r1,3);
            %#subnet
            nodes{count, 4} = 1;
            %weighting
            nodes{count, 5} = degr(r1);

            l1 = destrLabels{r1};
            nodes{count, 6} = l1;

            count = count+1;
        end
    end

    %build adja of particpating edges only
    adja_cut = adja(degr>0, degr>0);

    %save nodes
    fileID = fopen([path_local filename '.node'],'w');
    formatSpec = '%f %f %f %f %f %s \n';
    [nrows,ncols] = size(nodes);
    for row = 1:length(nodes)
        fprintf(fileID,formatSpec,nodes{row,:});
    end
    fclose(fileID);

    %save edges file
    save([path_local filename '.edge'], 'adja_cut', '-ascii');

    %call brainnetviewer
    %BrainNet_MapCfg(filename1, filename2...)
    BrainNet_MapCfg(fullfile(brainnet_dir, '/Data/SurfTemplate/BrainMesh_ICBM152.nv'), [path_local filename '.node'], [path_local filename '.edge'])
    input('Press Enter For next Network')
end