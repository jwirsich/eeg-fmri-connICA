function viewICAcomp(comp, configs, label, isTurnSign, savePath)
% view ICA component (and save figure)
%
% Jonathan Wirsich, Enrico Amico 2020
%
% Wirsich, J., Amico, E., Giraud A.L. Goñi, J, Sadaghiani S.,2020 
% Multi-timescale hybrid components of the functional brain connectome: A bimodal EEG-fMRI decomposition
% Network Neuroscience (2020) 4 (3): 658–677. https://doi.org/10.1162/netn_a_00135

    compIndex = 1;

    figures = [];
    labels=[];

    labels{end+1} = [label '_aMean'];
    figures(end+1) = figure('name', labels{end});
    if ~isTurnSign
        a_mean = comp.avg.Amean;
    else
        a_mean = -comp.avg.Amean;
        comp.avg.matrixfMRI = -comp.avg.matrixfMRI;
        comp.avg.matrixEEG = -comp.avg.matrixEEG;
        comp.avg.vectorfMRI = -comp.avg.vectorfMRI;
        comp.avg.vectorEEG = -comp.avg.vectorEEG;
    end
    ub_fMRI = prctile(comp.avg.vectorfMRI,95);
    lb_fMRI = prctile(comp.avg.vectorfMRI,5);
    ub_EEG = prctile(comp.avg.vectorEEG,95);
    lb_EEG = prctile(comp.avg.vectorEEG,5);

    % ub_fMRI = prctile(comp.avg.vectorfMRI,100);
    % lb_fMRI = prctile(comp.avg.vectorfMRI,0);
    % ub_EEG = prctile(comp.avg.vectorEEG,100);
    % lb_EEG = prctile(comp.avg.vectorEEG,0);

    reflect_folder = fileparts(mfilename('fullpath'));
    tmp = load(fullfile(reflect_folder, '/../../data/aparc_a2009_yeoRS7_148reg_eeg_nosubc.mat'));
    yeoROIs_eeg = tmp.yeoROIs;

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
    labels{end+1} = [label 'connICA_comp'];
    figures(end+1) = figure('name', labels{end}, 'Renderer', 'painters', 'Position', [10 10 900 600]);%(configs.h2);
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

    if ~strcmp(savePath, '')
        %print figures
        for i = 1:length(figures)
           filename = [savePath ...
                   labels{i}];
           print(figures(i), '-dpng', filename)
        end
    end

end