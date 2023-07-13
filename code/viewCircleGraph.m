% Show cirular graphs of IC components
%
% Jonathan Wirsich, Enrico Amico 2020
%
% Wirsich, J., Amico, E., Giraud A.L. Goñi, J, Sadaghiani S.,2020 
% Multi-timescale hybrid components of the functional brain connectome: A bimodal EEG-fMRI decomposition
% Network Neuroscience (2020) 4 (3): 658–677. https://doi.org/10.1162/netn_a_00135
%
% includes
% Please download from https://github.com/paul-kassebaum-mathworks/circularGraph
% TODO addpath('.../lib/circularGraph-master')

% load T_mask before exceuting (from view network.m)
reflect_folder = fileparts(mfilename('fullpath'));
path_local = fullfile(reflect_folder,'/../data/resu/PC75IC10_matrix_perc/');
%load destrieux labels
destrLabels=importdata(fullfile(reflect_folder, '/../data/destr_label_nosubc.txt'));
figures = [];
labels = [];

%frequency Main dataset
ica_labels{1} = 'Frequency_Main_PC75IC10';
%ICN Main dataset
ica_labels{2} = 'ICN_Main_PC75IC10';
%frequency Generalization dataset
ica_labels{3} = 'Frequency_Generalization_PC75IC10';
%ICN Generalization dataset
ica_labels{4} = 'ICN_Generalization_PC75IC10';

%iterate labels
for lab_it = 1:length(ica_labels)
    label_name = ['99_percentile_' (ica_labels{lab_it}) '_mask_T_eeg'];

    load([path_local label_name])

    path_git = '/media/jwirsich/DATA/projects/git/';

    %load yeoorder
    load(fullfile(reflect_folder, '../data/aparc_a2009_yeoRS7_148reg_eeg_nosubc.mat'))

    %put adja in yeoorder
    %masks are already in yeoorder!!!
    % adja_yeo = mask_T_fmri;
    adja_yeo = mask_T_eeg;
    yeoROIs_eeg = yeoROIs;
    yeoOrder_eeg = yeoOrder;

    regions =  148;
    sameRSN = zeros(regions,regions);
    for i=1:length(yeoROIs_eeg)
        for j=1:length(yeoROIs_eeg)
            if(yeoROIs_eeg(i)==yeoROIs_eeg(j))
                sameRSN(i,j) = yeoROIs_eeg(i);
                sameRSN(j,i) = yeoROIs_eeg(i);
            end
        end
    end

    temp = sameRSN(yeoOrder_eeg,yeoOrder_eeg);
    colorMap_circ = zeros(regions, 3);
    for i = 1:regions
        switch temp(i,i)
            case 1
                colorMap_circ(i, :) = [0 1 0];
    %             colorMap_circ(i, :) = [1 0 0];
    %             colorMap_circ(i, :) = [0.5 0 0.5];
    %             colorMap_circ(i, :) = [0 1 1];
            case 2
                colorMap_circ(i, :) = [0 0 1];
    %             colorMap_circ(i, :) = [0 0 1];
    %             colorMap_circ(i, :) = [0 1 0];
            case 3
                colorMap_circ(i, :) = [0.5 0.5 0.5];
    %             colorMap_circ(i, :) = [0 0 1];
    %             colorMap_circ(i, :) = [0 1 0];
            case 4
    %             colorMap_circ(i, :) = [1 1 0];
                colorMap_circ(i, :) = [0 1 1 ];
            case 5
    %             colorMap_circ(i, :) = [0 1 1 ];
                colorMap_circ(i, :) = [0.5 0 0.5];
            case 6
                colorMap_circ(i, :) = [1 1 0];
    %             colorMap_circ(i, :) = [1 0.5 0];
    %             colorMap_circ(i, :) = [1 0 1];
    %             colorMap_circ(i, :) = [0.5 0.5 0.5];
            case 7
                colorMap_circ(i, :) = [1 0 0];
        end
    end

    %     RSN_labels = {'VIS' 'SM' 'DA' 'VA' 'L' 'FP' 'DMN' 'SUBC'};
    RSN_labels = {'VIS' 'SM' 'DA' 'VA' 'L' 'FP' 'DMN'};

    labels_circ = destrLabels(yeoOrder_eeg);
    last = '';
    for i = 1:regions
        if ~strcmp(last, RSN_labels{temp(i,i)})
            labels_circ{i} = RSN_labels{temp(i,i)};
        else
            labels_circ{i} = '';
        end
        last = RSN_labels{temp(i,i)};
    %     labels_circ{i} = [labels_circ{i}];
    %     labels_circ{i} = [RSN_labels{temp(i,i)} '-' labels_circ{i}];
    %     labels_circ{i} = [''];
    end

    %put adja in same rsn (non-yeo order)
    temp(adja_yeo>0) = -1;

    n_steps = 0;
    for i=1:max(yeoROIs_eeg)
        n_steps = n_steps + nnz(yeoROIs_eeg==i);
        ticks_space(i) = n_steps - nnz(yeoROIs_eeg==i)/2;
    end

    yticks = ticks_space;
    labels{end+1} = [label_name 'yeo_adja'];
    figures(end+1) = figure('name', labels{end}); imagesc(temp); axis square; colorbar;
    set(gca,'XTick',yticks,'XTickLabel', RSN_labels); xtickangle(45);
    set(gca,'YTick',yticks,'YTickLabel', RSN_labels);
    set(gca,'fontsize',16)

    %build a circular graph
    labels{end+1} = [label_name '_circ_yeo_adja'];
    figures(end+1) = figure('name', labels{end});
    %TODO load library before

    noRSNs = max(yeoROIs_eeg);
    sizeRSN = zeros(noRSNs, 1);
    for r = 1:noRSNs
        sizeRSN(r) = sum(yeoROIs_eeg==r);
    end

    %shift defaultmode
    shift_dmn_first = mod((1:regions) + sizeRSN(end),regions);
    shift_dmn_first(shift_dmn_first==0) = regions;

    circularGraph(adja_yeo, 'ColorMap', ...
        colorMap_circ, 'Label', labels_circ);

    adja = adja_yeo;
    countRSNs = zeros(noRSNs);
    countRSN=0;
    %count RSN labels
    temp_same = sameRSN(yeoOrder_eeg,yeoOrder_eeg);
    for r1 = 1:regions-1
        %read connections of adja
        for r2 = r1+1:regions
            if temp_same(r1,r2)>0 && adja(r1,r2)==1
                countRSN = countRSN+1;
                countRSNs(temp_same(r1,r2), temp_same(r1,r2)) =  countRSNs(temp_same(r1,r2), temp_same(r1,r2))+1;
            elseif adja(r1,r2)==1
                countRSNs(temp_same(r1,r1), temp_same(r2,r2)) =  countRSNs(temp_same(r1,r1), temp_same(r2,r2))+1;               
                countRSNs(temp_same(r2,r2), temp_same(r1,r1)) =  countRSNs(temp_same(r2,r2), temp_same(r1,r1))+1;               
            end
        end
    end

    disp(sum(adja_yeo(:)>0))

    %noramalize by number of possible connections: no_edges(RSN1)*no_edges(RSN2)
    % countRSNs_norm = zeros(noRSNs);
    % for r1 = 1:noRSNs
    %     for r2 = r1:noRSNs
    %         countRSNs_norm(r1, r2) = countRSNs(r1,r2)/(sizeRSN(r1));
    %         countRSNs_norm(r2, r1) = countRSNs_norm(r1, r2);
    % %             disp([RSN_labels{r1} '-' RSN_labels{r2} ': ' num2str(countRSNs_norm(r1, r2))])
    %         disp([RSN_labels{r1} '-' RSN_labels{r2} ': ' num2str(countRSNs(r1, r2))])
    %     end
    % end

    yticks = [1:8];
    labels{end+1} = [label_name '-RSNs'];
    figures(end+1) = figure('name', labels{end});
    imagesc(countRSNs);  axis square; colorbar;
    set(gca,'XTick',yticks,'XTickLabel', RSN_labels); xtickangle(45);
    set(gca,'YTick',yticks,'YTickLabel', RSN_labels);

    % labels{end+1} = [bands{b} '-RSNs -normed'];
    % figures(end+1) = figure('name', labels{end});
    % imagesc(countRSNs_norm);  axis square; colorbar;
    % set(gca,'XTick',yticks,'XTickLabel', RSN_labels); xtickangle(45);
    % set(gca,'YTick',yticks,'YTickLabel', RSN_labels);

    %count intra
    %count inter

    % print figures
    % for i = 1:length(figures)
    %    filename = ['/media/jwirsich/DATA/projects/eeg_fmri_connICA/connICA/resu/PC75IC10_circ/' ...
    %            labels{i}];
    %    print(figures(i), '-dpng', filename)
    % end

end