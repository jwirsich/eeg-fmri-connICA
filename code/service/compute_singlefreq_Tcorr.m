function [T_corr_fmri,T_corr_eeg,A_corr,Trait_names] = compute_singlefreq_Tcorr(comp_multifreq,singleFreq_dir,freq_index)
% function compute_singlefreq_Tcorr
% compare a specific frequency trait + weights found on multiband
% with the single frequency ones to check for matches
% T_corr = numtraits x 2; correlation between traits and weights
% freq_index = 1 delta, 2 theta, 3 alpha,4 beta, 5 gamma
%
%
% Jonathan Wirsich, Enrico Amico 2020
%
% Wirsich, J., Amico, E., Giraud A.L. Goñi, J, Sadaghiani S.,2020 
% Multi-timescale hybrid components of the functional brain connectome: A bimodal EEG-fMRI decomposition
% Network Neuroscience (2020) 4 (3): 658–677. https://doi.org/10.1162/netn_a_00135
%
List_traits= dir(fullfile(singleFreq_dir,'*T*mat'));
num_traits = length(List_traits);
Trait_names = cell(1,num_traits);
T_corr_fmri = zeros(num_traits,2);
T_corr_eeg = zeros(num_traits,2);
A_corr = zeros(num_traits,1);
for t=1:num_traits
    Trait_names{t} = List_traits(t).name;
    load(fullfile(singleFreq_dir,Trait_names{t}));
    A_sF = comp.avg.Amean;
    num_subj = length(A_sF);
    freq_range = ((freq_index-1)*num_subj)+1:(num_subj*freq_index);
    A_mF = comp_multifreq.avg.Amean(freq_range);
    T_vec_sF_fmri = comp.avg.vectorfMRI';
    T_vec_sF_eeg = comp.avg.vectorEEG';    
    T_vec_mF_fmri = comp_multifreq.avg.vectorfMRI';
    T_vec_mF_eeg = comp_multifreq.avg.vectorEEG';
    Th_eeg_sF = prctile(abs(T_vec_sF_eeg),95);
    Th_eeg_mF = prctile(abs(T_vec_mF_eeg),95);
    
    Th_fmri_sF = prctile(abs(T_vec_sF_fmri),95);
    Th_fmri_mF = prctile(abs(T_vec_mF_fmri),95);
    
    T_corr_fmri(t,1) = corr(T_vec_sF_fmri,T_vec_mF_fmri);
    T_corr_eeg(t,1) = corr(T_vec_sF_eeg,T_vec_mF_eeg);
    T_corr_fmri(t,2) = corr(abs(T_vec_sF_fmri(abs(T_vec_sF_fmri)>Th_fmri_sF)),abs(T_vec_mF_fmri(abs(T_vec_mF_fmri)>Th_fmri_mF)));
    T_corr_eeg(t,2) = corr(abs(T_vec_sF_eeg(abs(T_vec_sF_eeg)>Th_eeg_sF)),abs(T_vec_mF_eeg(abs(T_vec_mF_eeg)>Th_eeg_mF)));
    
    A_corr(t) = corr(A_mF',A_sF');
end
disp('fMRI (2 cols, no mask - mask) - EEG (2 cols no mask - mask)) - A or weights (1 col)')
disp([T_corr_fmri T_corr_eeg A_corr])

