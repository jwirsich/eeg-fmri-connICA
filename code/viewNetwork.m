% View connectivity weightings of final IC components of Wirsich et al. 2020
%
% Jonathan Wirsich, Enrico Amico 2020
%
% Wirsich, J., Amico, E., Giraud A.L. Goñi, J, Sadaghiani S.,2020 
% Multi-timescale hybrid components of the functional brain connectome: A bimodal EEG-fMRI decomposition
% Network Neuroscience (2020) 4 (3): 658–677. https://doi.org/10.1162/netn_a_00135

reflect_folder = fileparts(mfilename('fullpath'));
load(fullfile(reflect_folder, '../data/aparc_a2009_yeoRS7_148reg_eeg_nosubc.mat'))

%load selected ICA
%frequency Main
selected{1} = fullfile(reflect_folder, '../data/ICs/main/IC10_PCA75/HybridT1.mat');
labels{1} = 'Frequency_Main_PC75IC10';
isTurned(1) = 0;
%ICN Main
selected{2} = fullfile(reflect_folder, '../data/ICs/main/IC10_PCA75/HybridT7.mat');
labels{2} = 'ICN_Main_PC75IC10';
isTurned(2) = 1;
%frequency Generalization
selected{3} = fullfile(reflect_folder, '../data/ICs/generalization/IC10_PCA75/HybridT3.mat');
labels{3} = 'Frequency_Generalization_PC75IC10';
isTurned(3) = 1;
%ICN Generalization
selected{4} = fullfile(reflect_folder, '../data/ICs/generalization/IC10_PCA75/HybridT8.mat');
labels{4} = 'ICN_Generalization_PC75IC10';
isTurned(4) = 1;

figures = [];
fig_labels = [];

for sel = 1:4
load(selected{sel})
    
N = 148;
mask = triu(true(N,N),1);

if isTurned(sel)
    T_eeg = -comp.avg.matrixEEG; % convention: keep positive weights
    T_fmri = -comp.avg.matrixfMRI; % convention: keep positive weights
else
    T_eeg = comp.avg.matrixEEG; % convention: keep positive weights
    T_fmri = comp.avg.matrixfMRI; % convention: keep positive weights
end

% define 5,95 procentile boundaries
%lb_eeg = prctile(T_eeg(mask),5);
ub_eeg = prctile(abs(T_eeg(mask)),99);
%lb_fmri = prctile(T_fmri(mask),5);
ub_fmri = prctile(abs(T_fmri(mask)),99);
% create eeg fmri masks
mask_T_eeg = zeros(size(T_eeg));
mask_T_eeg(T_eeg>ub_eeg)= 1;
%mask_T_eeg(T_eeg<lb_eeg)= -1;
mask_T_fmri = zeros(size(T_fmri));
mask_T_fmri(T_fmri>ub_fmri)= 1;
%mask_T_fmri(T_fmri<lb_fmri)= -1;

yeoROIs_eeg = yeoROIs;
yeoOrder_eeg = yeoOrder;
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
%     figure, imagesc(sameRSN(yeoOrder_eeg,yeoOrder_eeg)); axis square; colorbar;
%     set(gca,'XTick',yticks,'XTickLabel', RSN_labels); xtickangle(45);
%     set(gca,'YTick',yticks,'YTickLabel', RSN_labels);
%     set(gca,'fontsize',12)
    
%fMRI
fig_labels{end+1} = ['fMRI_' (labels{sel})];
figures{end+1} = figure('name', fig_labels{end})
imagesc(T_fmri); axis square; colormap jet; colorbar
set(gca,'XTick',yticks,'XTickLabel', RSN_labels); xtickangle(45);
set(gca,'YTick',yticks,'YTickLabel', RSN_labels);
set(gca,'fontsize',12)
hline(onset_ticks(1:end-1), 'k')
vline(onset_ticks(1:end-1), 'k')

%EEG
fig_labels{end+1} = ['EEG_' (labels{sel})];
figures{end+1} = figure ('name', fig_labels{end}), 
imagesc(T_eeg); axis square; colormap jet; colorbar;
set(gca,'XTick',yticks,'XTickLabel', RSN_labels); xtickangle(45);
set(gca,'YTick',yticks,'YTickLabel', RSN_labels);
set(gca,'fontsize',12)
hline(onset_ticks(1:end-1), 'k')
vline(onset_ticks(1:end-1), 'k')

fig_labels{end+1} = ['99_percentile_' (labels{sel})];
figures{end+1} = figure ('name', fig_labels{end}), 
subplot(1,2,1);spy(mask_T_fmri); axis square
set(gca,'XTick',yticks,'XTickLabel', RSN_labels); xtickangle(45);
set(gca,'YTick',yticks,'YTickLabel', RSN_labels);
set(gca,'fontsize',12)
hline(onset_ticks(1:end-1), 'k')
vline(onset_ticks(1:end-1), 'k')

subplot(1,2,2);spy(mask_T_eeg); axis square
set(gca,'XTick',yticks,'XTickLabel', RSN_labels); xtickangle(45);
set(gca,'YTick',yticks,'YTickLabel', RSN_labels);
set(gca,'fontsize',12)
hline(onset_ticks(1:end-1), 'k')
vline(onset_ticks(1:end-1), 'k')

%save maskT_fMRI / EEG

% save(['.../connICA/resu/PC75IC10_matrix_perc/99_percentile_' (labels{sel}) '_mask_T_fmri'], 'mask_T_fmri')
% save(['.../connICA/resu/PC75IC10_matrix_perc/99_percentile_' (labels{sel}) '_mask_T_eeg'] , 'mask_T_eeg')

end

%print figures
% for i = 1:length(figures)
%    filename = ['/media/jwirsich/DATA/projects/eeg_fmri_connICA/connICA/resu/PC75IC10_matrix_perc/' ...
%            fig_labels{i}];
%    print(figures(i), '-dpng', filename)
% end