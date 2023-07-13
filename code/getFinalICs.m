%% Get and visualize the final EEG-fMRI ICs
%
% Jonathan Wirsich, Enrico Amico 2020
%
% Wirsich, J., Amico, E., Giraud A.L. Goñi, J, Sadaghiani S.,2020 
% Multi-timescale hybrid components of the functional brain connectome: A bimodal EEG-fMRI decomposition
% Network Neuroscience (2020) 4 (3): 658–677. https://doi.org/10.1162/netn_a_00135

%load final ICs
reflect_folder = fileparts(mfilename('fullpath'));
paths{1} = [reflect_folder filesep '..' filesep 'data' filesep  'ICs' filesep 'main' filesep];
paths{2} = [reflect_folder filesep '..' filesep 'data' filesep  'ICs' filesep 'generalization' filesep];

%final selected IC according to Wirsich et al. 2020, NetNeuro
finalparam = 'IC10_PCA75';

load([reflect_folder filesep '..' filesep 'data' filesep 'aparc_a2009_yeoRS7_148reg_eeg_nosubc.mat'])

%configure path if you want to save figures
out_path = '';

tmp = load([paths{1} finalparam filesep 'HybridT1']);
viewICAcomp(tmp.comp, tmp.configs, 'Main Frequency', 0, out_path)
tmp = load([paths{1} finalparam filesep 'HybridT7']);
viewICAcomp(tmp.comp, tmp.configs, 'Main RSN', 1, out_path)

tmp = load([paths{2} finalparam filesep 'HybridT3']);
viewICAcomp(tmp.comp, tmp.configs, 'Generalization Frequency', 1, out_path)
tmp = load([paths{2} finalparam filesep 'HybridT8']);
viewICAcomp(tmp.comp, tmp.configs, 'Generalization RSN', 1, out_path)