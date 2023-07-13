% View input mean connectome of connectomes used in the connICA matrix
%
% Jonathan Wirsich, Enrico Amico 2020
%
% Wirsich, J., Amico, E., Giraud A.L. Goñi, J, Sadaghiani S.,2020 
% Multi-timescale hybrid components of the functional brain connectome: A bimodal EEG-fMRI decomposition
% Network Neuroscience (2020) 4 (3): 658–677. https://doi.org/10.1162/netn_a_00135

reflect_folder = fileparts(mfilename('fullpath'));
data_path = fullfile(reflect_folder, '../data/');
regions = 148;

% load([data_path 'connICA_generalizationData_allfreq'])
load([data_path 'connICA_mainData_allfreq'])

figure
imagesc([connICA.connICA_matrixfMRI connICA.connICA_matrixDelta])
vline(10878.5, 'black')

colormap('parula')
% hline([14.5 28.5 42.5 56.5], {'black', 'black', 'black', 'black', 'black'})
% vline(10878.5, 'black')

figure
imagesc([connICA.connICA_matrixfMRI connICA.connICA_matrixTheta])
vline(10878.5, 'black')
colormap('parula')

figures = [];
labels = [];

[figures(end+1), labels{end+1}] = viewConnectome4Triu(mean(connICA.connICA_matrixfMRI(:,:)), regions, 'fMRI', 1);
[figures(end+1), labels{end+1}] = viewConnectome4Triu(mean(connICA.connICA_matrixDelta(:,:)), regions, 'delta', 1);
[figures(end+1), labels{end+1}] = viewConnectome4Triu(mean(connICA.connICA_matrixTheta(:,:)), regions, 'theta', 1);
[figures(end+1), labels{end+1}] = viewConnectome4Triu(mean(connICA.connICA_matrixAlpha(:,:)), regions, 'alpha', 1);
[figures(end+1), labels{end+1}] = viewConnectome4Triu(mean(connICA.connICA_matrixBeta(:,:)), regions, 'beta', 1);
[figures(end+1), labels{end+1}] = viewConnectome4Triu(mean(connICA.connICA_matrixGamma(:,:)), regions, 'gamma', 1);

% saveFigures('.../connICA/resu/meanconnect_Main/', figures, labels)
