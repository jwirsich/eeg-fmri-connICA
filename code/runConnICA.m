% Calculate EEG-fMRI ICA from multimodal connectomes
%
% Jonathan Wirsich, Enrico Amico 2020
%
% Wirsich, J., Amico, E., Giraud A.L. Goñi, J, Sadaghiani S.,2020 
% Multi-timescale hybrid components of the functional brain connectome: A bimodal EEG-fMRI decomposition
% Network Neuroscience (2020) 4 (3): 658–677. https://doi.org/10.1162/netn_a_00135

%imports
%download and link FastICA from https://research.ics.aalto.fi/ica/fastica/code/dlcode.shtml
addpath('TODO/FastICA_25')

reflect_folder = fileparts(mfilename('fullpath'));
git_path = [reflect_folder filesep '..' filesep];
configs.path = git_path;

%define outpath where to write components
outpath = 'TODO';
%define if looking at main or generalization data
isMainData = 0;
%preload matching index
preloadMatchingIndex = 1;
%explore several IC/PCA combinations
exploreICspace = 0;
%Hybrid: use all Frequencies
isHybrid = 1;
%if not hybrid define frequency band 1-5: delta-gamma
singleFrequency = 0;

if preloadMatchingIndex == 1
    if isMainData == 1
        load([git_path 'data' filesep 'connICA_mainData_allfreq.mat'])
    else
        load([git_path 'data' filesep 'connICA_generalizationData_allfreq.mat'])
    end
    connICA_matrixfMRI = connICA.connICA_matrixfMRI;
    connICA_matrixDelta = connICA.connICA_matrixDelta;
    connICA_matrixTheta = connICA.connICA_matrixTheta;
    connICA_matrixAlpha = connICA.connICA_matrixAlpha;
    connICA_matrixBeta = connICA.connICA_matrixBeta;
    connICA_matrixGamma = connICA.connICA_matrixGamma;
else
    if isMainData == 1
        connICA_matrixfMRI = loadConnICA_meansess(git_path, 'destrieux', 'fMRI');
        connICA_matrixDelta = loadConnICA_meansess(git_path, 'destrieux', 'delta');
        connICA_matrixTheta = loadConnICA_meansess(git_path, 'destrieux', 'theta');
        connICA_matrixAlpha = loadConnICA_meansess(git_path, 'destrieux', 'alpha');
        connICA_matrixBeta = loadConnICA_meansess(git_path, 'destrieux', 'beta');
        connICA_matrixGamma = loadConnICA_meansess(git_path, 'destrieux', 'gamma');
    else
        connICA_matrixfMRI = loadConnICA_generalizationData(git_path, 'fMRI');
        connICA_matrixDelta = loadConnICA_generalizationData(git_path, 'delta');
        connICA_matrixTheta = loadConnICA_generalizationData(git_path, 'theta');
        connICA_matrixAlpha = loadConnICA_generalizationData(git_path, 'alpha');
        connICA_matrixBeta = loadConnICA_generalizationData(git_path, 'beta');
        connICA_matrixGamma = loadConnICA_generalizationData(git_path, 'gamma');
    end
end

if isHybrid
    connICA_matrix = [connICA_matrixfMRI connICA_matrixDelta; connICA_matrixfMRI connICA_matrixTheta; ...
                      connICA_matrixfMRI connICA_matrixAlpha; connICA_matrixfMRI connICA_matrixBeta; ... 
                      connICA_matrixfMRI connICA_matrixGamma];
else
    switch singleFrequency
        case 1
            connICA_matrix = [connICA_matrixfMRI connICA_matrixDelta];
        case 2 
            connICA_matrix = [connICA_matrixfMRI connICA_matrixTheta];
        case 3 
            connICA_matrix = [connICA_matrixfMRI connICA_matrixAlpha];
        case 4
            connICA_matrix = [connICA_matrixfMRI connICA_matrixBeta];
        case 5
            connICA_matrix = [connICA_matrixfMRI connICA_matrixGamma];
        otherwise
            throw(MException('connICA:freqNotDefined', ['Frequency ' int2str(singleFrequency) 'not defined']))
    end
end

%number of subjects
if isMainData == 1
    configs.numSubj = 26;
else
    configs.numSubj = 14;
end
configs.nFreq = 5;

connICA.connICA_matrixfMRI = connICA_matrixfMRI;
connICA.connICA_matrixDelta = connICA_matrixDelta;
connICA.connICA_matrixTheta = connICA_matrixTheta;
connICA.connICA_matrixAlpha = connICA_matrixAlpha;
connICA.connICA_matrixBeta = connICA_matrixBeta;
connICA.connICA_matrixGamma = connICA_matrixGamma;

connICA.config = configs;

if isMainData
    save([git_path 'data' filesep 'connICA_mainData_allfreq.mat'], 'connICA');
else
    save([git_path 'data' filesep 'connICA_generalizationData_allfreq.mat'], 'connICA');
end

%alter connICA parameters for parameter exploration
if exploreICspace == 1
    noICs = [5 10 15 20];
    pcathrshlds = [0.75 0.8 0.85 0.9];
else
    noICs = [10];
    pcathrshlds = [0.75];
end

for it = 1:length(noICs)
    configs.numOfIC = noICs(it);
    for it2 = 1:length(pcathrshlds)
        configs.VarExp = pcathrshlds(it2);
        
        configs.output = [outpath ...
            'IC' int2str(noICs(it)) '_PCA' int2str(pcathrshlds(it2)*100)];
       
        stack = cell(0,0);
        try
            calc_connICA_hybrid(connICA_matrix, configs);
        catch exception
            disp(['IC' int2str(noICs(it)) '_PCA' num2str(pcathrshlds) ' skipped due to ' exception.message])
            stack{end+1} = exception.stack;
            rethrow(exception)
        end
    end
end