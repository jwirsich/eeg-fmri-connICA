%% This batch runs PCA cleaning on orig_matrix = num functional connectome x num functional edges (upper triangular, as in connICA),
%% the FCs are organized in this order: Subj1_Vis1, Subj1_Vis2....SubjN_Vis1, SubjN_Vis2.
%
% Jonathan Wirsich, Enrico Amico 2020
%
% Wirsich, J., Amico, E., Giraud A.L. Goñi, J, Sadaghiani S.,2020 
% Multi-timescale hybrid components of the functional brain connectome: A bimodal EEG-fMRI decomposition
% Network Neuroscience (2020) 4 (3): 658–677. https://doi.org/10.1162/netn_a_00135
%
function [PCA_clean_matrix,numPCAComps] = run_PCA_cleaning_singlerun(connICA_matrix,VarExp) 
%% PCA cleaning
numFCs = size(connICA_matrix,1);
numPCAComps = size(connICA_matrix,1);
[COEFFS, SCORE, latent] = pca(connICA_matrix','NumComponents',numPCAComps); 
% COEFFS are the EIGENVECTORS, SCORE are the PROJECTIONS OF THE DATA IN THE EIGENVECS SUBSPACE 
variance = cumsum(latent)./sum(latent); 
variance = variance(1:numPCAComps); %explained variance with the selected num of PCA comps
figure, plot(1:numFCs,variance,'ok');
numPCAComps=find(variance>=VarExp,1); % perform PCA on the data retaining 90 % of the variance
[COEFFS, SCORE, latent] = pca(connICA_matrix','NumComponents',numPCAComps); 
PCA_clean_matrix = SCORE * COEFFS'; % PCA reconstructed demeaned data
PCA_clean_matrix = bsxfun(@plus, PCA_clean_matrix,mean(connICA_matrix')); % plug the mean back
PCA_clean_matrix = PCA_clean_matrix'; %back to subjects x edges form
