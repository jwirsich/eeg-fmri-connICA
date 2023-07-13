function [ICC_a_fbands,ICC_a_subj, p_a_fbands, p_a_subj] = f_ICC_edgewise_v1(A,configs)
% Compute ICC
% Evaluate ICC on the a_mean (mean over ICA runs) weights of hybrid connICA 
% need to Add ICC to the path
% 
% Jonathan Wirsich, Enrico Amico 2020
%
% Wirsich, J., Amico, E., Giraud A.L. Goñi, J, Sadaghiani S.,2020 
% Multi-timescale hybrid components of the functional brain connectome: A bimodal EEG-fMRI decomposition
% Network Neuroscience (2020) 4 (3): 658–677. https://doi.org/10.1162/netn_a_00135

numFreq = configs.nFreq;
numSubj = configs.numSubj;
data4icc = reshape(A,[numSubj,numFreq]); % judges are frequency bands, transposed judges == subjects
rows2delete_A = isnan(sum(data4icc,2));
data4icc(rows2delete_A,:) = [];
[ICC_a_fbands,~,~,~,~,~,p_a_fbands] = ICC(data4icc','1-1') ;
[ICC_a_subj,~,~,~,~,~,p_a_subj] = ICC(data4icc,'1-1');    

return;