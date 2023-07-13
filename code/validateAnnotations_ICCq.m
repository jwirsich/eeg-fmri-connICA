% calculate ICC and q for all ICs
%
% Jonathan Wirsich, Enrico Amico 2020
%
% Wirsich, J., Amico, E., Giraud A.L. Goñi, J, Sadaghiani S.,2020 
% Multi-timescale hybrid components of the functional brain connectome: A bimodal EEG-fMRI decomposition
% Network Neuroscience (2020) 4 (3): 658–677. https://doi.org/10.1162/netn_a_00135

reflect_folder = fileparts(mfilename('fullpath'));
load(fullfile(reflect_folder, '/../data/aparc_a2009_yeoRS7_148reg_eeg_nosubc.mat'))
yeoOrder_eeg = yeoOrder;
yeoROIs_eeg = yeoROIs(yeoOrder_eeg);
% addpath('/media/jwirsich/DATA/projects/git/multimodal-private/eeg-fmri-dyna')
% close all;
%loop through folders
paths{1} = fullfile(reflect_folder, '/../data/ICs/main/');
paths{2} = fullfile(reflect_folder, '/../data/ICs/generalization/');

path_data_name{1} = 'main';
path_data_name{2} = 'generalization';

icn_stack_group = cell(2,4,4);
freq_stack_group = cell(2,4,4);

figures = [];
labels = [];

for path_it = 1:2

disp(paths{path_it})    
    
path = paths{path_it};
dirs = dir(path);
dirIndex = find([dirs.isdir]);
%cut directories /. and /..
dirIndex = dirIndex(3:end);

IC_PCA_matrix_icn = zeros(4);
IC_PCA_matrix_freq = zeros(4);

icn_stack = cell(4,4);
freq_stack = cell(4,4);

count = 0;

for i = 1:length(dirIndex)
    IC_PCA_no = sscanf(dirs(dirIndex(i)).name, 'IC%d_PCA%d');
    %IC 5 to 20
    %pca 75 to 90 
    %convert to index
    IC_PCA_no(1) = IC_PCA_no(1)/5;
    IC_PCA_no(2) = (IC_PCA_no(2)-70)/5;
        
    %validate if file exists
    if exist([path dirs(dirIndex(i)).name filesep 'select.txt'], 'file') == 2
        %get annotations
        select = importdata([path dirs(dirIndex(i)).name filesep 'select.txt']);
        
        %find cannonicals
        dim = size(select.data);
        freq = cell(0,0);
        icn = cell(0,0);
        for i2 = 1:dim(1)
            if select.data(i2,4) == 1
                if select.data(i2,1) == 1
                    %load cannoncicals
                    freq{end+1} = load([path dirs(dirIndex(i)).name filesep 'HybridT' int2str(select.data(i2,2))]);
                    %flip signs 
                    if select.data(i2,3) == 1
                        tmp = freq{end};
                        tmp.comp.avg.matrixEEG = -tmp.comp.avg.matrixEEG;
                        tmp.comp.avg.matrixfMRI = -tmp.comp.avg.matrixfMRI;
                        tmp.comp.avg.vectorEEG = -tmp.comp.avg.vectorEEG;
                        tmp.comp.avg.vectorfMRI = -tmp.comp.avg.vectorfMRI;
                        tmp.comp.avg.Amean = -tmp.comp.avg.Amean;
                        freq{end} = tmp;
                    end
                else
                    icn{end+1} = load([path dirs(dirIndex(i)).name filesep 'HybridT' int2str(select.data(i2,2))]);
                    disp([int2str(length(icn)) ': ' path dirs(dirIndex(i)).name filesep 'HybridT' int2str(select.data(i2,2))])
                    %flip signs 
                    if select.data(i2,3) == 1
                        tmp = icn{end};
                        tmp.comp.avg.matrixEEG = -tmp.comp.avg.matrixEEG;
                        tmp.comp.avg.matrixfMRI = -tmp.comp.avg.matrixfMRI;
                        tmp.comp.avg.vectorEEG = -tmp.comp.avg.vectorEEG;
                        tmp.comp.avg.vectorfMRI = -tmp.comp.avg.vectorfMRI;
                        tmp.comp.avg.Amean = -tmp.comp.avg.Amean;
                        icn{end} = tmp;
                    end
                end
            end
        end
        %in case more than one canoncial merge
        aggreagation_level_icn = length(icn);
        if length(icn) > 1
            avg.comp.avg.matrixEEG = 0;
            avg.comp.avg.matrixfMRI = 0;
            avg.comp.avg.vectorEEG = 0;
            avg.comp.avg.vectorfMRI = 0;
            avg.comp.avg.Amean = 0;
            for it_icn = 1:length(icn)
                tmp = icn{it_icn};
                avg.comp.avg.matrixEEG = avg.comp.avg.matrixEEG + tmp.comp.avg.matrixEEG/length(icn);
                avg.comp.avg.matrixfMRI = avg.comp.avg.matrixfMRI + tmp.comp.avg.matrixfMRI/length(icn);
                avg.comp.avg.vectorEEG = avg.comp.avg.vectorEEG + tmp.comp.avg.vectorEEG/length(icn);
                avg.comp.avg.vectorfMRI = avg.comp.avg.vectorfMRI + tmp.comp.avg.vectorfMRI/length(icn);
                avg.comp.avg.Amean = avg.comp.avg.Amean + tmp.comp.avg.Amean/length(icn);
            end
            avg.configs = tmp.configs;
            icn_final{1} = avg;
            IC_PCA_matrix_icn(IC_PCA_no(1), IC_PCA_no(2)) = 2;
        elseif ~isempty(icn)
            icn_final = icn;
            IC_PCA_matrix_icn(IC_PCA_no(1), IC_PCA_no(2)) = 1;
        else
            icn_final = icn;
            IC_PCA_matrix_icn(IC_PCA_no(1), IC_PCA_no(2)) = 0;
        end
        
        aggreagation_level_freq = length(freq);
        if length(freq) > 1
            avg.comp.avg.matrixEEG = 0;
            avg.comp.avg.matrixfMRI = 0;
            avg.comp.avg.vectorEEG = 0;
            avg.comp.avg.vectorfMRI = 0;
            for it_freq = 1:length(freq)
                tmp = freq{it_freq};
                avg.comp.avg.matrixEEG = avg.comp.avg.matrixEEG + tmp.comp.avg.matrixEEG/length(freq);
                avg.comp.avg.matrixfMRI = avg.comp.avg.matrixfMRI + tmp.comp.avg.matrixfMRI/length(freq);
                avg.comp.avg.vectorEEG = avg.comp.avg.vectorEEG + tmp.comp.avg.vectorEEG/length(freq);
                avg.comp.avg.vectorfMRI = avg.comp.avg.vectorfMRI + tmp.comp.avg.vectorfMRI/length(freq);
                avg.comp.avg.Amean = avg.comp.avg.Amean + tmp.comp.avg.Amean/length(freq);
            end
            avg.configs = tmp.configs;
            freq_final{1} = avg;
            IC_PCA_matrix_freq(IC_PCA_no(1), IC_PCA_no(2)) = 2;
        elseif ~isempty(freq)
            freq_final = freq;
            IC_PCA_matrix_freq(IC_PCA_no(1), IC_PCA_no(2)) = 1;
        else
            freq_final = freq;
            IC_PCA_matrix_freq(IC_PCA_no(1), IC_PCA_no(2)) = 0;
        end
        
        if ~isempty(freq_final)
            freq_stack(IC_PCA_no(1), IC_PCA_no(2)) = freq_final;
            freq_stack_group(path_it, IC_PCA_no(1), IC_PCA_no(2)) = freq_final;
        end
        if ~isempty(icn_final)
            icn_stack(IC_PCA_no(1), IC_PCA_no(2)) = icn_final;
            icn_stack_group(path_it, IC_PCA_no(1), IC_PCA_no(2)) = icn_final;
        end
        
        %visualize and validate
        if ~isempty(freq_final)
            label = [dirs(dirIndex(i)).name ' Frequency ' int2str(aggreagation_level_freq)];
            viewICAcomp(freq_final{1}.comp, freq_final{1}.configs, label, 0, '')
            [ICC_a_fbands,ICC_a_subj, p_a_fbands, p_a_subj] = f_ICC_edgewise_v1(freq_final{1}.comp.avg.Amean,freq_final{1}.configs);
            
            [ci, q] = hardwiredModularityLouvain(freq_final{1}.comp.avg.matrixfMRI, 1, yeoROIs_eeg, 'negative_asym');
            [ci, q_eeg] = hardwiredModularityLouvain(freq_final{1}.comp.avg.matrixEEG, 1, yeoROIs_eeg, 'negative_asym');
            
            disp([label ': ICC freq ' num2str(ICC_a_fbands) ' (p=' num2str(p_a_fbands)  ') - ICC subj ' ...
                num2str(ICC_a_subj) ' (p=' num2str(p_a_subj)  ')'])
            disp([label ': q fMRI ' num2str(q) ' - q EEG ' num2str(q_eeg)])
        end
        if ~isempty(icn_final)
            label = [dirs(dirIndex(i)).name ' ICN ' int2str(aggreagation_level_icn)];
%             viewICAcomp(icn_final{1}.comp, icn_final{1}.configs, label, 0);
            [ICC_a_fbands,ICC_a_subj, p_a_fbands, p_a_subj] = f_ICC_edgewise_v1(icn_final{1}.comp.avg.Amean,icn_final{1}.configs);

            [ci, q] = hardwiredModularityLouvain(icn_final{1}.comp.avg.matrixfMRI, 1, yeoROIs_eeg, 'negative_asym');
            [ci, q_eeg] = hardwiredModularityLouvain(icn_final{1}.comp.avg.matrixEEG, 1, yeoROIs_eeg, 'negative_asym');
            
            disp([label ': ICC freq ' num2str(ICC_a_fbands) ' (p=' num2str(p_a_fbands)  ') - ICC subj ' ...
                num2str(ICC_a_subj) ' (p=' num2str(p_a_subj)  ')'])
            
            %TODO pvalue for q!!!
            disp([label ': q fMRI ' num2str(q) ' - q EEG ' num2str(q_eeg)])
        end
        
    else
        select = [];
    end

end

ylabels = {'IC5', 'IC10', 'IC15', 'IC20'};
xlabels = {'PCA75', 'PCA80', 'PCA85', 'PCA90'};
ticks_space = [1 2 3 4];
yticks = ticks_space;


%get correlation matrix

%TODO put back in later
% figure('name', 'freq')
% imagesc(IC_PCA_matrix_freq)
% set(gca,'XTick',yticks,'XTickLabel', xlabels); xtickangle(45);
% set(gca,'YTick',yticks,'YTickLabel', ylabels);
% figure('name', 'icn')
% imagesc(IC_PCA_matrix_icn)
% set(gca,'XTick',yticks,'XTickLabel', xlabels); xtickangle(45);
% set(gca,'YTick',yticks,'YTickLabel', ylabels);


%loop over bot datasets and get correlation

correlfMRI = zeros(16);
correlEEG = zeros(16);

correlfMRI_freq = zeros(16);
correlEEG_freq = zeros(16);

count1 = 0;

for i1 = 1:4
    for i2 = 1:4
        count1 = count1+1;
        count2 = 0;
        
        stack_labels{count1} = ['IC' int2str(i1*5) 'PCA' int2str(i2*5+70)];
    
        for it1 = 1:4
            for it2 = 1:4
                count2 = count2+1;
%                 if ~isempty(freq_stack{i1, i2}) && ~isempty(freq_stack{it1, it2})
                if ~isempty(icn_stack{i1, i2}) && ~isempty(icn_stack{it1, it2})
                    cor = corr(icn_stack{i1, i2}.comp.avg.vectorfMRI', icn_stack{it1, it2}.comp.avg.vectorfMRI');
                    correlfMRI(count1, count2) = cor;
                    
                    cor_eeg = corr(icn_stack{i1, i2}.comp.avg.vectorEEG', icn_stack{it1, it2}.comp.avg.vectorEEG');
                    correlEEG(count1, count2) = cor_eeg;
                end
                if ~isempty(freq_stack{i1, i2}) && ~isempty(freq_stack{it1, it2})
                    cor_freq = corr(freq_stack{i1, i2}.comp.avg.vectorfMRI', freq_stack{it1, it2}.comp.avg.vectorfMRI');                    
                    correlfMRI_freq(count1, count2) = cor_freq;
                    
                    cor_freq_eeg = corr(freq_stack{i1, i2}.comp.avg.vectorEEG', freq_stack{it1, it2}.comp.avg.vectorEEG');
                    correlEEG_freq(count1, count2) = cor_freq_eeg;
                    
%                     cor = corr(icn_stack{i1, i2}.comp.avg.Amean', icn_stack{it1, it2}.comp.avg.Amean');
                end
            end
        end
    end
end

%should be a loop fanncy fancy ;)
labels{end+1} = [path_data_name{path_it} '_inter_correl_ICN_fMRI'];
figures(end+1) = figure('name', labels{end});
imagesc(correlfMRI); colormap parula; colorbar; axis square; 
caxis([-0.2 1]);
yticks = 1:16;
set(gca,'XTick',yticks,'XTickLabel', stack_labels); xtickangle(45);
set(gca,'YTick',yticks,'YTickLabel', stack_labels);

labels{end+1} = [path_data_name{path_it} '_inter_correl_ICN_EEG'];
figures(end+1) = figure('name', labels{end});
imagesc(correlEEG); colormap parula; colorbar; axis square; 
caxis([-0.2 1]);
yticks = 1:16;
set(gca,'XTick',yticks,'XTickLabel', stack_labels); xtickangle(45);
set(gca,'YTick',yticks,'YTickLabel', stack_labels);

labels{end+1} = [path_data_name{path_it} '_inter_correl_freq_fMRI'];
figures(end+1) = figure('name', labels{end});
imagesc(correlfMRI_freq); colormap parula; colorbar; axis square; 
caxis([-0.2 1]);
yticks = 1:16;
set(gca,'XTick',yticks,'XTickLabel', stack_labels); xtickangle(45);
set(gca,'YTick',yticks,'YTickLabel', stack_labels);

labels{end+1} = [path_data_name{path_it} '_inter_correl_freq_EEG'];
figures(end+1) = figure('name', labels{end});
imagesc(correlEEG_freq); colormap parula; colorbar; axis square; 
caxis([-0.2 1]);

% caxis([0 0.6]);

yticks = 1:16;
set(gca,'XTick',yticks,'XTickLabel', stack_labels); xtickangle(45);
set(gca,'YTick',yticks,'YTickLabel', stack_labels);

end

count1 = 0;
correlfMRI_ICN = zeros(16);
correlEEG_ICN = zeros(16);
correlfMRI_freq = zeros(16);
correlEEG_freq = zeros(16);

p_correlfMRI_ICN = zeros(16);
p_correlEEG_ICN = zeros(16);
p_correlfMRI_freq = zeros(16);
p_correlEEG_freq = zeros(16);

for i1 = 1:4
    for i2 = 1:4
        count1 = count1+1;
        count2 = 0;
        
        stack_labels{count1} = ['IC' int2str(i1*5) 'PCA' int2str(i2*5+70)];
    
        for it1 = 1:4
            for it2 = 1:4
                count2 = count2+1;
                if ~isempty(icn_stack_group{1, i1, i2}) && ~isempty(icn_stack_group{2, it1, it2})
                    [cor, pval] = corr(icn_stack_group{1, i1, i2}.comp.avg.vectorfMRI', icn_stack_group{2, it1, it2}.comp.avg.vectorfMRI');
                    
%                     cor = corr(icn_stack_group{1, i1, i2}.comp.avg.Amean', icn_stack_group{2, it1, it2}.comp.avg.Amean');
                    correlfMRI_ICN(count1, count2) = cor;
                    p_correlfMRI_ICN(count1, count2) = pval;
                    
                    [corEEG, pval_eeg] = corr(icn_stack_group{1, i1, i2}.comp.avg.vectorEEG', icn_stack_group{2, it1, it2}.comp.avg.vectorEEG');
                    correlEEG_ICN(count1, count2) = corEEG;
                    p_correlEEG_ICN(count1, count2) = pval_eeg;
                end
                
                if ~isempty(freq_stack_group{1, i1, i2}) && ~isempty(freq_stack_group{2, it1, it2})
                    [cor, pval] = corr(freq_stack_group{1, i1, i2}.comp.avg.vectorfMRI', freq_stack_group{2, it1, it2}.comp.avg.vectorfMRI');
%                     cor = corr(freq_stack_group{1, i1, i2}.comp.avg.vectorfMRI', freq_stack_group{2, it1, it2}.comp.avg.vectorfMRI');
                    correlfMRI_freq(count1, count2) = cor;
                    p_correlfMRI_freq(count1, count2) = pval;
                    
                    [corEEG, pval_eeg] = corr(freq_stack_group{1, i1, i2}.comp.avg.vectorEEG', freq_stack_group{2, it1, it2}.comp.avg.vectorEEG');
%                     cor = corr(freq_stack_group{1, i1, i2}.comp.avg.vectorfMRI', freq_stack_group{2, it1, it2}.comp.avg.vectorfMRI');
                    correlEEG_freq(count1, count2) = corEEG;
                    p_correlEEG_freq(count1, count2) = pval_eeg;
                end
                                
            end
        end
    end
end

labels{end+1} = 'group_fMRI_ICN';
figures(end+1) = figure('name', labels{end});
imagesc(correlfMRI_ICN); colormap parula; colorbar; axis square; 
caxis([-0.2 1]);
yticks = 1:16;
set(gca,'XTick',yticks,'XTickLabel', stack_labels); xtickangle(45);
set(gca,'YTick',yticks,'YTickLabel', stack_labels);

labels{end+1} = 'group_EEG_ICN';
figures(end+1) = figure('name', labels{end});
imagesc(correlEEG_ICN); colormap parula; colorbar; axis square; 
caxis([-0.2 1]);
yticks = 1:16;
set(gca,'XTick',yticks,'XTickLabel', stack_labels); xtickangle(45);
set(gca,'YTick',yticks,'YTickLabel', stack_labels);

labels{end+1} = 'group_EEG_freq';
figures(end+1) = figure('name', labels{end});
imagesc(correlEEG_freq); colormap parula; colorbar; axis square; 
caxis([-0.2 1]);
yticks = 1:16;
set(gca,'XTick',yticks,'XTickLabel', stack_labels); xtickangle(45);
set(gca,'YTick',yticks,'YTickLabel', stack_labels);

labels{end+1} = 'group_fMRI_freq';
figures(end+1) = figure('name', labels{end});
imagesc(correlfMRI_freq); colormap parula; colorbar; axis square; 
caxis([-0.2 1]);
yticks = 1:16;
set(gca,'XTick',yticks,'XTickLabel', stack_labels); xtickangle(45);
set(gca,'YTick',yticks,'YTickLabel', stack_labels);


disp('ICN-FG-fMRI')
max_fMRI_fg = max(correlfMRI_ICN(correlfMRI_ICN>0));
[row,col] = find(correlfMRI_ICN==max_fMRI_fg);
disp(['Max-fMRI-FG Main: ' stack_labels{row} ' General:' stack_labels{col} ...
    ' - ' num2str(max_fMRI_fg) ', p=' num2str(p_correlfMRI_ICN(row, col))])
min_fMRI_fg = min(correlfMRI_ICN(correlfMRI_ICN>0));
[row,col] = find(correlfMRI_ICN==min_fMRI_fg);
disp(['Min-fMRI-FG Main: ' stack_labels{row} ' General:' stack_labels{col} ...
    ' - ' num2str(min_fMRI_fg) ', p=' num2str(p_correlfMRI_ICN(row, col))])
disp(num2str(mean(correlfMRI_ICN(correlfMRI_ICN>0))));

disp('ICN-FG-EEG')
max_EEG_fg = max(correlEEG_ICN(correlEEG_ICN>0));
[row,col] = find(correlEEG_ICN==max_EEG_fg);
disp(['Max-EEG-FG Main: ' stack_labels{row} ' General:' stack_labels{col} ...
    ' - ' num2str(max_EEG_fg) ', p=' num2str(p_correlEEG_ICN(row, col))])
min_EEG_fg = min(correlEEG_ICN(correlEEG_ICN>0));
[row,col] = find(correlEEG_ICN==min_EEG_fg);
disp(['Min-EEG-FG Main: ' stack_labels{row} ' General:' stack_labels{col} ...
    ' - ' num2str(min_EEG_fg) ', p=' num2str(p_correlEEG_ICN(row, col))])
disp(num2str(mean(correlEEG_ICN(correlEEG_ICN>0))));

disp('VIS-FS-fMRI')
max_fMRI_fs = max(correlfMRI_freq(correlfMRI_freq>0));
[row,col] = find(correlfMRI_freq==max_fMRI_fs);
disp(['Max-fMRI-FS Main: ' stack_labels{row} ' General:' stack_labels{col} ...
    ' - ' num2str(max_fMRI_fs) ', p=' num2str(p_correlfMRI_freq(row, col))])
min_fMRI_fs = min(correlfMRI_freq(correlfMRI_freq>0));
[row,col] = find(correlfMRI_freq==min_fMRI_fs);
disp(['Min-fMRI-FS Main: ' stack_labels{row} ' General:' stack_labels{col} ...
    ' - ' num2str(min_fMRI_fs) ', p=' num2str(p_correlfMRI_freq(row, col))])
disp(num2str(mean(correlfMRI_freq(correlfMRI_freq>0))));

disp('VIS-FS-EEG')
max_EEG_fs = max(correlEEG_freq(correlEEG_freq>0));
[row,col] = find(correlEEG_freq==max_EEG_fs);
disp(['Max-EEG-FS Main: ' stack_labels{row} ' General:' stack_labels{col} ...
    ' - ' num2str(max_EEG_fs) ', p=' num2str(p_correlEEG_freq(row, col))])
min_EEG_fs = min(correlEEG_freq(correlEEG_freq>0));
[row,col] = find(correlEEG_freq==min_EEG_fs);
disp(['Min-EEG-FS Main: ' stack_labels{row} ' General:' stack_labels{col} ...
    ' - ' num2str(min_EEG_fs) ', p=' num2str(p_correlEEG_freq(row, col))])
disp(num2str(mean(correlEEG_freq(correlEEG_freq>0))));



%chose one component test accross datasets

% print figures
% for i = 1:length(figures)
%    filename = ['/media/jwirsich/DATA/projects/eeg_fmri_connICA/connICA/resu/datasetfit/' ...
%            labels{i}];
%    print(figures(i), '-dpng', filename)
% end