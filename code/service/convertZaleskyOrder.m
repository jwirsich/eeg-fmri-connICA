% convert connectome data from [Wirsich 2017, NIMG] to standard Freesurfer order
%
% Jonathan Wirsich, Enrico Amico 2020
%
% Wirsich, J., Amico, E., Giraud A.L. Goñi, J, Sadaghiani S.,2020 
% Multi-timescale hybrid components of the functional brain connectome: A bimodal EEG-fMRI decomposition
% Network Neuroscience (2020) 4 (3): 658–677. https://doi.org/10.1162/netn_a_00135
%
reflect_folder = fileparts(mfilename('fullpath'));
data_path = fullfile(reflect_folder, '/../../data');
regions = 148;
subj_no = 14;
load([data_path filesep 'eeg_fmri_connectomes_destrieux_scrubbed_generalizationData']);
%load Modules derived in Wirsich 2017 NeruoImage using aproach from Zalesky 2014 PNAS
load([data_path filesep 'rsn_modules_zaleski.mat']);
%construct permutation vector
perm_modules = double.empty;
for sub = 1:length(subnets)
    perm_modules = [perm_modules subnets{sub}];
end

[aa orig_order] = sort(perm_modules);
modal_ind = {'fMRI' 'delta' 'theta' 'alpha' 'beta' 'gamma'};
mask_triu = triu(true(regions,regions),1);
allRegress_orig = zeros(subj_no,length(modal_ind),nnz(mask_triu));

for j = 1:length(modal_ind)
    for i = 1:subj_no 
        conn1 = Connectome(getConnect_generalization(squeeze(allRegress(i,:,:)), modal_ind{j}), regions);
        conn_mat = conn1.getMatrix;
        conn_mat = conn_mat(orig_order, orig_order);
        allRegress_orig(i,j,:) = conn_mat(mask_triu);
    end
end

%you may save result to
%eeg_fmri_connectomes_destrieux_scrubbed_generalizationData_origorder_triu.mat