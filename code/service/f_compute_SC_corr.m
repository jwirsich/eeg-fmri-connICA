function SC_corr = f_compute_SC_corr(SC)
% Compute matching index 
% see Rubinov and Sporns 2010 and Amico and Goni NetNeuroSc 2018
% 
% Jonathan Wirsich, Enrico Amico 2020
%
% Wirsich, J., Amico, E., Giraud A.L. Goñi, J, Sadaghiani S.,2020 
% Multi-timescale hybrid components of the functional brain connectome: A bimodal EEG-fMRI decomposition
% Network Neuroscience (2020) 4 (3): 658–677. https://doi.org/10.1162/netn_a_00135
%
N = size(SC,1);
SC_corr = zeros(N,N);
for i = 1:N-1
    for j = i+1:N
        vec_i = SC(i,:);
        vec_i([i,j])=[];
        vec_j = SC(j,:);
        vec_j([i,j])=[];
        SC_corr(i,j) = corr(vec_i',vec_j');
    end
end
SC_corr = SC_corr + SC_corr';
        