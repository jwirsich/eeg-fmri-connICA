% Loop through preselected ICs (frequency and ICN network)
%
% Jonathan Wirsich, Enrico Amico 2020
%
% Wirsich, J., Amico, E., Giraud A.L. Goñi, J, Sadaghiani S.,2020 
% Multi-timescale hybrid components of the functional brain connectome: A bimodal EEG-fMRI decomposition
% Network Neuroscience (2020) 4 (3): 658–677. https://doi.org/10.1162/netn_a_00135

%loop through folders
reflect_folder = fileparts(mfilename('fullpath'));
paths{1} = fullfile(reflect_folder, '/../data/ICs/main/');
paths{2} = fullfile(reflect_folder, '/../data/ICs/generalization/');

icn_stack_group = cell(2,4,4);
freq_stack_group = cell(2,4,4);

for path_it = 1:2

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
%             viewICAcomp(freq_final{1}.comp, freq_final{1}.configs, label, 0)
        end
        if ~isempty(icn_final)
            label = [dirs(dirIndex(i)).name ' ICN ' int2str(aggreagation_level_icn)];
%             viewICAcomp(icn_final{1}.comp, icn_final{1}.configs, label, 0);
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
figure('name', 'freq')
imagesc(IC_PCA_matrix_freq)
set(gca,'XTick',yticks,'XTickLabel', xlabels); xtickangle(45);
set(gca,'YTick',yticks,'YTickLabel', ylabels);
figure('name', 'icn')
imagesc(IC_PCA_matrix_icn)
set(gca,'XTick',yticks,'XTickLabel', xlabels); xtickangle(45);
set(gca,'YTick',yticks,'YTickLabel', ylabels);
%loop over bot datasets and get correlation

correlfMRI = zeros(16);

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
%                     cor = corr(freq_stack{i1, i2}.comp.avg.vectorfMRI', freq_stack{it1, it2}.comp.avg.vectorfMRI');
                    cor = corr(icn_stack{i1, i2}.comp.avg.vectorfMRI', icn_stack{it1, it2}.comp.avg.vectorfMRI');
%                     cor = corr(icn_stack{i1, i2}.comp.avg.vectorEEG', icn_stack{it1, it2}.comp.avg.vectorEEG');
%                     cor = corr(icn_stack{i1, i2}.comp.avg.Amean', icn_stack{it1, it2}.comp.avg.Amean');
                    correlfMRI(count1, count2) = cor;
                end
            end
        end
    end
end

figure('name', int2str(path_it))
imagesc(correlfMRI); 
% caxis([0 0.6]);

yticks = 1:16;
set(gca,'XTick',yticks,'XTickLabel', stack_labels); xtickangle(45);
set(gca,'YTick',yticks,'YTickLabel', stack_labels);

end

count1 = 0;
correlfMRI = zeros(16);
for i1 = 1:4
    for i2 = 1:4
        count1 = count1+1;
        count2 = 0;
        
        stack_labels{count1} = ['IC' int2str(i1*5) 'PCA' int2str(i2*5+70)];
    
        for it1 = 1:4
            for it2 = 1:4
                count2 = count2+1;
%                 if ~isempty(freq_stack_group{1, i1, i2}) && ~isempty(freq_stack_group{2, it1, it2})
                if ~isempty(icn_stack_group{1, i1, i2}) && ~isempty(icn_stack_group{2, it1, it2})
%                     cor = corr(freq_stack_group{1, i1, i2}.comp.avg.vectorfMRI', freq_stack_group{2, it1, it2}.comp.avg.vectorfMRI');
                     cor = corr(icn_stack_group{1, i1, i2}.comp.avg.vectorfMRI', icn_stack_group{2, it1, it2}.comp.avg.vectorfMRI');
%                     cor = corr(icn_stack_group{1, i1, i2}.comp.avg.vectorEEG', icn_stack_group{2, it1, it2}.comp.avg.vectorEEG');
%                     cor = corr(icn_stack_group{1, i1, i2}.comp.avg.Amean', icn_stack_group{2, it1, it2}.comp.avg.Amean');
                    correlfMRI(count1, count2) = cor;
                end
            end
        end
    end
end

figure('name', 'group')
imagesc(correlfMRI)
yticks = 1:16;
set(gca,'XTick',yticks,'XTickLabel', stack_labels); xtickangle(45);
set(gca,'YTick',yticks,'YTickLabel', stack_labels);


%chose one component test accross datasets