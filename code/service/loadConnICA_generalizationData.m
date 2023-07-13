function connICA_matrix = loadConnICA_generalizationData(inpath, modal, doTransfo)
%load data from Generalization dataset
% 
% Jonathan Wirsich, Enrico Amico 2020
%
% Wirsich, J., Amico, E., Giraud A.L. Goñi, J, Sadaghiani S.,2020 
% Multi-timescale hybrid components of the functional brain connectome: A bimodal EEG-fMRI decomposition
% Network Neuroscience (2020) 4 (3): 658–677. https://doi.org/10.1162/netn_a_00135

    data_path = [inpath 'data'];
    regions = 148; %destrieux
    dat = load([data_path filesep 'eeg_fmri_connectomes_destrieux_scrubbed_generalizationData_origorder_triu.mat']);

    connICA_matrix = zeros(1, sum(1:regions-1));
    mask_triu = triu(true(regions,regions),1);
    subj_no = size(dat.allRegress_orig, 1);
    modal_ind = {'fMRI' 'delta' 'theta' 'alpha' 'beta' 'gamma'};

    for i = 1:subj_no 
        conn_mat = zeros(regions,regions);

        X = squeeze(dat.allRegress_orig(i, contains(modal_ind,modal)==1, :));
        conn_mat(mask_triu)= X;
        conn_mat = conn_mat + conn_mat';

        %convert EEG data to macthing index
        if doTransfo == 1
            if ~strcmp(modal, 'fMRI')
                conn_mat = f_compute_SC_corr(conn_mat);
            end
        end

        connICA_matrix(i, :) = conn_mat(mask_triu);
    end

end