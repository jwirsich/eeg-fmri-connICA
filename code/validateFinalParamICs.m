% get ICC and q for final ICs
%
% Jonathan Wirsich, Enrico Amico 2020
%
% Wirsich, J., Amico, E., Giraud A.L. Goñi, J, Sadaghiani S.,2020 
% Multi-timescale hybrid components of the functional brain connectome: A bimodal EEG-fMRI decomposition
% Network Neuroscience (2020) 4 (3): 658–677. https://doi.org/10.1162/netn_a_00135

reflect_folder = fileparts(mfilename('fullpath'));
load(fullfile(reflect_folder, '/../data/aparc_a2009_yeoRS7_148reg_eeg_nosubc.mat'))
paths{1} = fullfile(reflect_folder, '/../data/ICs/main/');
paths{2} = fullfile(reflect_folder, '/../data/ICs/generalization/');

% finalparam = 'IC10_PCA80';
finalparam = 'IC10_PCA75';

yeoOrder_eeg = yeoOrder;
%data is already transformed so transform the yeo idx too!
yeoROIs_eeg = yeoROIs(yeoOrder_eeg);

%Paris
ICs{1} = [1 2 6 7]; % IC10PC75
% ICs{1} = [10 8]; % IC10PC80 selected (1 of 2 rsn)
% ICs{1} = [1 2 3 4 5 7 8 9 10];
%Marseille
ICs{2} = [3 8]; % IC10PC75
% ICs{2} = [4 3]; % IC10PC80
% ICs{2} = [1 3 4 6 7 8 9 10];

for path_it = 1:2
    disp(paths{path_it})
    
    for i = 1:length(ICs{path_it})
        IC_idx = ICs{path_it};
        label = int2str(IC_idx(i));
        tmp_all{path_it, i} = load([paths{path_it} finalparam filesep 'HybridT' int2str(IC_idx(i))]);
        
%         disp([paths{path_it} finalparam filesep 'HybridT' int2str(IC_idx(i))]);
        
        %define if to turn %TODO this is the better test in
        %validateAnnotations_ICCq.m it is read out from the .txt
        if mean(tmp_all{path_it, i}.comp.avg.Amean) < 0 
           tmp_all{path_it, i}.comp.avg.Amean = -tmp_all{path_it, i}.comp.avg.Amean;
           tmp_all{path_it, i}.comp.avg.matrixfMRI = -tmp_all{path_it, i}.comp.avg.matrixfMRI;
           tmp_all{path_it, i}.comp.avg.matrixEEG = -tmp_all{path_it, i}.comp.avg.matrixEEG;
           tmp_all{path_it, i}.comp.avg.vectorEEG = -tmp_all{path_it, i}.comp.avg.vectorEEG;
           tmp_all{path_it, i}.comp.avg.vectorfMRI = -tmp_all{path_it, i}.comp.avg.vectorfMRI;
        end
        tmp = tmp_all{path_it, i};
        
        [ICC_a_fbands,ICC_a_subj, p_a_fbands, p_a_subj] = f_ICC_edgewise_v1(tmp.comp.avg.Amean,tmp.configs);

        [ci, q] = hardwiredModularityLouvain(tmp.comp.avg.matrixfMRI, 1, yeoROIs_eeg, 'negative_asym');
        [ci, q_eeg] = hardwiredModularityLouvain(tmp.comp.avg.matrixEEG, 1, yeoROIs_eeg, 'negative_asym');
        [cor_EEGfMRI p] = corr(tmp.comp.avg.vectorfMRI', tmp.comp.avg.vectorEEG');
        disp([label ': q fMRI ' num2str(q) ' - q EEG ' num2str(q_eeg) ...
            '/ ICC freq ' num2str(ICC_a_fbands) ' (p= ' num2str(p_a_fbands) ') - ICC subj ' num2str(ICC_a_subj) ...
            ' (p= ' num2str(p_a_subj) ') Corr: ' num2str(cor_EEGfMRI) ' p:' num2str(p)])
        
        %scatter EEG vs. fMRI and correlate
%         imagesc
%         scatter(tmp.comp.avg.vectorfMRI, tmp.comp.avg.vectorEEG)
    end
end

for i1 = 1:length(ICs{1})
    IC_idx1 = ICs{1};
    for i2 = 1:length(ICs{2})
        IC_idx2 = ICs{2};
        label = [int2str(IC_idx1(i1)) '-' int2str(IC_idx2(i2))];
        [cor_EEGEEG, p_EEG] = corr(tmp_all{1, i1}.comp.avg.vectorEEG', tmp_all{2, i2}.comp.avg.vectorEEG');
        [cor_fMRIfMRI, p_fMRI] = corr(tmp_all{1, i1}.comp.avg.vectorfMRI', tmp_all{2, i2}.comp.avg.vectorfMRI');
        disp([label ' cor EEG-EEG: ' num2str(cor_EEGEEG) ' p=' num2str(p_EEG) ' cor fMRI-fMRI: ' num2str(cor_fMRIfMRI) ' p=' num2str(p_EEG)])
    end
end