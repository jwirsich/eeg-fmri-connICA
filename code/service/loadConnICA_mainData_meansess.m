function connICA_matrix = loadConnICA_mainData_meansess(inpath, atlas, modal)
%load data from main dataset
% 
% Jonathan Wirsich, Enrico Amico 2020
%
% Wirsich, J., Amico, E., Giraud A.L. Goñi, J, Sadaghiani S.,2020 
% Multi-timescale hybrid components of the functional brain connectome: A bimodal EEG-fMRI decomposition
% Network Neuroscience (2020) 4 (3): 658–677. https://doi.org/10.1162/netn_a_00135

data_path = [inpath 'data'];

if strcmp(atlas, 'destrieux')
    regions = 148;
    load([data_path filesep 'eeg_fmri_connectomes_destrieux_scrubbed_mainData.mat']);
    load([data_path filesep 'aparc_a2009_yeoRS7_148reg_eeg.mat']);
end

count = 0;

connICA_matrix = zeros(1, sum(1:regions-1));

%try to load precalculated
for i = 1:length(subj)
    
    X_mean = zeros(length(subj(i).sess),size(connICA_matrix,2));
    
    for s = 1:length(subj(i).sess)
        conn1 = Connectome(getConnect(subj(i).sess(s), modal), regions);
        
        conn_mat = conn1.getMatrix;
        
        %convert EEg connectomes to matching index
        if ~strcmp(modal, 'fMRI')
            conn_mat = f_compute_SC_corr(conn_mat);
        end
        
        %convert to triu mask
        mask_ut = triu(true(regions,regions),1);
        X = conn_mat(mask_ut);
        
        X_mean(s,:) = X;
    end
    X_mean = mean(X_mean);
    count = count + 1;
    connICA_matrix(count, :) = X_mean;
    
end

end