%% batch example connICA (Amico et al., NeuroImage 2017)
%% Enrico Amico & Joaquin Goni, Purdue University
%% PLEASE CITE US!
% If you are using this code for your research, please kindly cite us:
% Mapping the functional connectome traits of levels of consciousness.
% Amico, E. ... & Goñi, J. (2017). NeuroImage, 148, 201-211.
% http://www.sciencedirect.com/science/article/pii/S1053811917300204
%
%% IMPORTANT: FastICA package IS NEEDED!
% Please download FastICA package
% https://research.ics.aalto.fi/ica/fastica/code/dlcode.shtml
% Version used in Wirsich, J., Amico, E., Giraud A.L. Goñi, J, Sadaghiani S.,2020 
% Multi-timescale hybrid components of the functional brain connectome: A bimodal EEG-fMRI decomposition, Netw Neurosci (just accepted), https://doi.org/10.1162/netn_a_00135
% for an updated version try: https://engineering.purdue.edu/ConnplexityLab/publications/connICA_toolbox_v1.1.tar.gz

%requires connICA matrix in workspace

% ConnICA params
configs.numRegions = 68; % Desikan
configs.mask = triu(true(configs.numRegions,configs.numRegions),1); % upper triangular mask
configs.epsilon = 0.0001; % epsilon default is 0.0001
configs.maxNumIterations = 1000; % maxNumIterations default is 1000
configs.maxFinetune = 100; % maxFinetune default is 100
configs.numConn = size(connICA_matrix,1);%number of connectome (subjects or several sessions)
configs.numRuns = 100; %100;%run several time for robust  
configs.numOfIC = 15; %independent component number (max ICC or take 15-20)    
% load(fullfile(pwd,'..','data','Glasser_ch2_yeo_RS7.mat'),'yeoOrder');
load('aparc_aseg_yeoRS7_68reg_eeg.mat')
configs.yeoOrder = yeoOrder_eeg; %TODO plug in later
icasig = nan(size(connICA_matrix,2),configs.numOfIC,configs.numRuns,'single');%2 is for colume %FC_traits single for optimizing
A = nan(size(connICA_matrix,1),configs.numOfIC,configs.numRuns,'single');%(A number of football players * component) number %mixing matrix= subjects weights (each subjects)
flags.PCAclean = 1;
if flags.PCAclean==1
    VarExp = 0.75; % var different for EEG or fMRI
    [Pre_ICA_matrix,numPCAComps] = run_PCA_cleaning_singlerun2(connICA_matrix,VarExp);
else
    Pre_ICA_matrix = connICA_matrix;
end

for i=1:configs.numRuns
    disp(i)
    [icasig_onerun,A_onerun,~] = fastica(Pre_ICA_matrix,'approach','symm','numOfIC',configs.numOfIC,'verbose','off',...
        'epsilon',configs.epsilon,'maxNumIterations',configs.maxNumIterations,...
        'maxFinetune',configs.maxFinetune);%running fastica
    A(:,:,i) = single(A_onerun);%weight
    icasig(:,:,i) = single(icasig_onerun');%component composition 1(component number)*EDGES
end
configs.minFreq = 0.75;%minum frequency required
configs.corrMin = 0.75;%correlation between FC_traits (matrix(different fastica time)
configs.corrMin_A = 0.75;%correlation between subject weights %don't go below 50%
[comp_match_run1, freq] = run_robust_connICA_fast(A,icasig,configs);%comp_mach_run1 is list of components that is robust with respect to the first(best) run to the connica
Prova=(freq>configs.minFreq);%test & index(num_conponent) of the number of robust components (ex comp1, comp2, comp16... is robust)
Index=find(Prova(:,1));
disp(Index) % number of robust components

%% Put back the robust traits in matrix form
for k=1:length(Index)
    compIndex=Index(k);%choose the component(order of component) we want to look at 
    figure;
    %subplot(2,4,subfig); hold on;
    icasig_comp = [];
    a0 = A(:,comp_match_run1(compIndex,1),1); %this is used as reference so that the same component across runs is always positively correlated (not inverted)
    for i=1:configs.numRuns
        if comp_match_run1(compIndex,i)>0
            a = A(:,comp_match_run1(compIndex,i),i);
            icasig_comp_one = squeeze(icasig(:,comp_match_run1(compIndex,i),i));
            if corr(a,a0)>0
                plot(a);
                icasig_comp(i,:) = icasig_comp_one;
            else
                plot(-a);
                icasig_comp(i,:) = -icasig_comp_one;
            end
        end
    end

    clear comp;
    comp.avg.vector = [];
    comp.avg.vector = mean(icasig_comp); %avg per column (across runs)
    comp.avg.matrix = zeros(configs.numRegions,configs.numRegions);
    comp.avg.matrix(configs.mask)=comp.avg.vector; %fill upper triangular
    comp.avg.matrix = comp.avg.matrix + (comp.avg.matrix'); %symmetrize matrix
    comp.avg.matrix = comp.avg.matrix(configs.yeoOrder,configs.yeoOrder); % apply yeoOrder


    xlabel('subjects');
    title(sprintf('connICA comp %d',compIndex))
    figure;%(configs.h2);
    %subplot(2,4,subfig); 
    imagesc(comp.avg.matrix,[-3,3]); colormap jet; colorbar; axis square;
    set(gca,'xtick',[]); set(gca,'ytick',[]); xlabel('regions'); ylabel('regions');
    title(sprintf('connICA comp %d',compIndex))
end

disp(['Number of selected PCs: ' num2str(numPCAComps)])