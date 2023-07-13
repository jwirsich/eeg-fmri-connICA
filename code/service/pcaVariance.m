% PCA variance check (half fMRI half EEG);
% PCA is imposed by the combined hybrid matrix
% 
% Jonathan Wirsich, Enrico Amico 2020
%
% Wirsich, J., Amico, E., Giraud A.L. Goñi, J, Sadaghiani S.,2020 
% Multi-timescale hybrid components of the functional brain connectome: A bimodal EEG-fMRI decomposition
% Network Neuroscience (2020) 4 (3): 658–677. https://doi.org/10.1162/netn_a_00135
%
%connICA_matrix is your hybrid full matrix
[COEFFS,SCORE,~,~,PCA.VARIANCE] = pca(connICA_matrix','NumComponents',size(connICA_matrix,1));
R2_FC = zeros(size(connICA_matrix,1),size(connICA_matrix,1)); 
R2_SC = zeros(size(connICA_matrix,1),size(connICA_matrix,1)); %this is fMRI (no of subjects x no PCA comps)

for s=1:size(connICA_matrix,1) %for each subject
    disp(s);
    Y_FC = connICA_matrix(s,1:nnz(configs.masktriu))';
    Y_SC = connICA_matrix(s,nnz(configs.masktriu)+1:end)';
    X_FC = [];
    X_SC = [];
    for j=1:size(connICA_matrix,1) %for each PCA comp
        if j==1
            X_FC = SCORE(1:nnz(configs.masktriu),j);
            X_SC = SCORE(nnz(configs.masktriu)+1:end,j);
        else
            X_FC = [X_FC SCORE(1:nnz(configs.masktriu),j)];
            X_SC = [X_SC SCORE(nnz(configs.masktriu)+1:end,j)];
        end
        [~,~,~,~,stats_FC] = regress(Y_FC,[ones(size(X_FC)) X_FC]);
        [~,~,~,~,stats_SC] = regress(Y_SC,[ones(size(X_SC)) X_SC]);
        R2_FC(s,j) = stats_FC(1);
        R2_SC(s,j) = stats_SC(1);
   end
end
ExpVar(1,:) = cumsum(PCA.VARIANCE)';
ExpVar_FC(1,:) = mean(R2_FC);
ExpVar_SC(1,:) = mean(R2_SC);
%should give figure S1 of the hybrid