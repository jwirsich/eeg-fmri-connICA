%  Extract RSN blocks from Desikan labelling
%
% Jonathan Wirsich, Enrico Amico 2020
%
% Wirsich, J., Amico, E., Giraud A.L. Goñi, J, Sadaghiani S.,2020 
% Multi-timescale hybrid components of the functional brain connectome: A bimodal EEG-fMRI decomposition
% Network Neuroscience (2020) 4 (3): 658–677. https://doi.org/10.1162/netn_a_00135
%
regions = 68;
sameRSN = zeros(regions,regions);
for i=1:length(yeoROIs_eeg)
    for j=1:length(yeoROIs_eeg)
        if(yeoROIs_eeg(i)==yeoROIs_eeg(j))
            sameRSN(i,j) = yeoROIs_eeg(i);
            sameRSN(j,i) = yeoROIs_eeg(i);
        end
    end
end

n_steps = 0;
for i=1:max(yeoROIs_eeg)
    n_steps = n_steps + nnz(yeoROIs_eeg==i);
    ticks_space(i) = n_steps - nnz(yeoROIs_eeg==i)/2;
end
RSN_labels = {'VIS' 'SM' 'DA' 'VA' 'L' 'FP' 'DMN' 'SUBC'};
yticks = ticks_space;
figure, imagesc(sameRSN(yeoOrder_eeg,yeoOrder_eeg)); axis square; colorbar;
set(gca,'XTick',yticks,'XTickLabel', RSN_labels); xtickangle(45);
set(gca,'YTick',yticks,'YTickLabel', RSN_labels);
set(gca,'fontsize',16)