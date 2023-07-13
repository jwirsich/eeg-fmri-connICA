%
% Jonathan Wirsich, Enrico Amico 2020
%
% Wirsich, J., Amico, E., Giraud A.L. Goñi, J, Sadaghiani S.,2020 
% Multi-timescale hybrid components of the functional brain connectome: A bimodal EEG-fMRI decomposition
% Network Neuroscience (2020) 4 (3): 658–677. https://doi.org/10.1162/netn_a_00135

load('/media/jwirsich/beckman/jonathan/EEG_fMRI_fromParis/connICA/hybrid_expl/IC10_PCA80/HybridT1.mat')
tmp = comp;
load('/media/jwirsich/beckman/jonathan/EEG_fMRI_fromParis/connICA/hybrid_expl/IC10_PCA80/HybridT10.mat')
tmp2 = comp;

comp.avg.Amean = -(tmp.avg.Amean+tmp2.avg.Amean)/2;
comp.avg.vectorfMRI = -(tmp.avg.vectorfMRI+tmp2.avg.vectorfMRI)/2;
comp.avg.vectorEEG = -(tmp.avg.vectorfMRI+tmp2.avg.vectorEEG)/2;
comp.avg.matrixEEG = -(tmp.avg.matrixEEG+tmp2.avg.matrixEEG)/2;
comp.avg.matrixfMRI = -(tmp.avg.matrixfMRI+tmp2.avg.matrixfMRI)/2;

compIndex = 1;

figures = [];
figures(end+1) = figure;
a_mean = comp.avg.Amean;
ub_fMRI = prctile(comp.avg.vectorfMRI,99);
lb_fMRI = prctile(comp.avg.vectorfMRI,1);
ub_EEG = prctile(comp.avg.vectorEEG,99);
lb_EEG = prctile(comp.avg.vectorEEG,1);
load('aparc_a2009_yeoRS7_148reg_eeg.mat')
yeoROIs_eeg = yeoROIs;

regions = 148;
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
    onset_ticks(i) = n_steps+0.5;
end
RSN_labels = {'VIS' 'SM' 'DA' 'VA' 'L' 'FP' 'DMN' 'SUBC'};
yticks = ticks_space;

freq_cols = {'b' 'c' 'g' 'y' 'r'};
hold on;
for f=1:configs.nFreq
    Subj_index = (configs.numSubj*(f-1)+1):f*configs.numSubj; 
    bar(Subj_index,a_mean(Subj_index),freq_cols{f});   
    title(sprintf('weights connICA comp %d',compIndex));
    xlabel('# Subjects'); ylabel('Weights');
end

xlabel('subjects');
title(sprintf('connICA comp %d',compIndex))
figures(end+1) = figure;%(configs.h2);
subplot(1,2,1); 
imagesc(comp.avg.matrixfMRI,[lb_fMRI,ub_fMRI]); colormap jet; colorbar; axis square;
set(gca,'xtick',[]); set(gca,'ytick',[]); xlabel('regions'); ylabel('regions');
title(sprintf('connICA comp %d fMRI',compIndex))
title(sprintf('connICA comp %d',compIndex))
%label RSNs
set(gca,'XTick',yticks,'XTickLabel', RSN_labels); xtickangle(45);
set(gca,'YTick',yticks,'YTickLabel', RSN_labels);
set(gca,'fontsize',12)
hline(onset_ticks(1:end-1), 'k')
vline(onset_ticks(1:end-1), 'k')
hold on
subplot(1,2,2); 
imagesc(comp.avg.matrixEEG,[lb_EEG,ub_EEG]); colormap jet; colorbar; axis square;
set(gca,'xtick',[]); set(gca,'ytick',[]); xlabel('regions'); ylabel('regions');
title(sprintf('connICA comp %d EEG',compIndex))
%label RSNs
set(gca,'XTick',yticks,'XTickLabel', RSN_labels); xtickangle(45);
set(gca,'YTick',yticks,'YTickLabel', RSN_labels);
set(gca,'fontsize',12)
hline(onset_ticks(1:end-1), 'k')
vline(onset_ticks(1:end-1), 'k')