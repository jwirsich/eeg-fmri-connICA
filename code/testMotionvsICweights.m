% test if ICA weights are influenced by motion (framewise displacement)
%
% Jonathan Wirsich, Enrico Amico 2020
%
% Wirsich, J., Amico, E., Giraud A.L. Goñi, J, Sadaghiani S.,2020 
% Multi-timescale hybrid components of the functional brain connectome: A bimodal EEG-fMRI decomposition
% Network Neuroscience (2020) 4 (3): 658–677. https://doi.org/10.1162/netn_a_00135
%

reflect_folder = fileparts(mfilename('fullpath'));
load(fullfile(reflect_folder, '../data/aparc_a2009_yeoRS7_148reg_eeg_nosubc.mat'))

paths{1} = fullfile(reflect_folder, '../data/ICs/main/');
paths{2} = fullfile(reflect_folder, '../data/ICs/generalization/');


fd_files{1} = fullfile(reflect_folder, '../data/eeg_fmri_connectomes_destrieux_scrubbed_mainData.mat');
ICs{1} = [1 2];

fd_files{2} = fullfile(reflect_folder, '../data/fd_generalizationData.mat');
ICs{2} = [3 8];

finalparam = 'IC10_PCA75';

for path_it = 1:2
  
    disp(paths{path_it})
    IC_idx = ICs{path_it};
    for it_ic = 1:length(ICs)
        disp(['HybridT' int2str(IC_idx(it_ic))]);
        load([paths{path_it} finalparam filesep 'HybridT' int2str(IC_idx(it_ic))]);
        %load FD
        tmp = load(fd_files{path_it});
        subj = tmp.subj;
    
        mean_fd = zeros(length(subj),1);
        scrubbed = zeros(length(subj),1);

        for i = 1:length(subj)
            stacked_fd = [];
            %load IC
            for s = 1:length(subj(i).sess)
                stacked_fd = [stacked_fd; subj(i).sess(s).fd_power];
            end
            %get scrubbed 
            scrubbed(i) = sum(stacked_fd > 0.5);
            mean_fd(i) = mean(stacked_fd);
        end

        %27-09-2019 TODO no error when enumeration was not found
        bands = enumeration('Bands');
        for b = 1:length(bands)
            disp(string(bands(b)));
            [r,p] = corr(comp.avg.Amean(length(subj)*(b-1)+1:length(subj)*b)', mean_fd, 'type', 'Spearman');
            disp(['r(fd-amean)= ' num2str(r) ' p=' num2str(p)]);

            [r,p] = corr(comp.avg.Amean(length(subj)*(b-1)+1:length(subj)*b)', scrubbed, 'type', 'Spearman');

            disp(['r(scrubbed-amean)= ' num2str(r) ' p=' num2str(p)]);

            figure('name', int2str(b))
            scatter(comp.avg.Amean(length(subj)*(b-1)+1:length(subj)*b)', mean_fd)
        end

        disp('All bands together')

        [r,p] = corr(comp.avg.Amean', [mean_fd; mean_fd; mean_fd; mean_fd; mean_fd], 'type', 'Spearman');
        disp(['r(fd-amean)= ' num2str(r) ' p=' num2str(p)]);

        [r,p] = corr(comp.avg.Amean', [scrubbed; scrubbed; scrubbed; scrubbed; scrubbed], 'type', 'Spearman');
        disp(['r(scrubbed-amean)= ' num2str(r) ' p=' num2str(p)]);

        figure('name', int2str(b))
        scatter(comp.avg.Amean, [mean_fd; mean_fd; mean_fd; mean_fd; mean_fd])
        scatter(comp.avg.Amean, [scrubbed; scrubbed; scrubbed; scrubbed; scrubbed])
        
   end
end