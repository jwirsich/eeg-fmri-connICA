function calc_connICA_hybrid (connICA_matrix, configs)
%% call Hybrid ICA worker
%
% connICA_matrix: Orig_matrix is 2D (nsubj,FCedges)
%
% Jonathan Wirsich, Enrico Amico 2020
%
% Wirsich, J., Amico, E., Giraud A.L. Goñi, J, Sadaghiani S.,2020 
% Multi-timescale hybrid components of the functional brain connectome: A bimodal EEG-fMRI decomposition
% Network Neuroscience (2020) 4 (3): 658–677. https://doi.org/10.1162/netn_a_00135

    % ConnICA params
    output = configs.output;
    configs.numRegions = 148; % Destrieux
    configs.mask = triu(true(configs.numRegions,configs.numRegions),1); % upper triangular mask
    configs.epsilon = 0.0001; % epsilon default is 0.0001
    configs.maxNumIterations = 1000; % maxNumIterations default is 1000
    configs.maxFinetune = 100; % maxFinetune default is 100
    configs.numConn = size(connICA_matrix,1);%number of connectome (subjects or several sessions)
    configs.numRuns = 500; %100;%run several time for robust  
    tmp = load([configs.path 'data' filesep 'aparc_a2009_yeoRS7_148reg_eeg.mat']);
    configs.yeoOrder = tmp.yeoOrder; %get yeo7 order of Destrieux regions
    
    icasig = nan(size(connICA_matrix,2),configs.numOfIC,configs.numRuns,'single');%2 is for colume %FC_traits single for optimizing
    A = nan(size(connICA_matrix,1),configs.numOfIC,configs.numRuns,'single');%(A number of football players * component) number %mixing matrix= subjects weights (each subjects)
    flags.PCAclean = 1;
    if flags.PCAclean==1
        VarExp = configs.VarExp; % var different for EEG or fMRI
        [Pre_ICA_matrix,numPCAComps] = run_PCA_cleaning_singlerun(connICA_matrix,VarExp);
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
    figures = [];
    ICC_a_fbands=nan(1,length(Index));
    ICC_a_subj=nan(1,length(Index));
    for k=1:length(Index)
        compIndex=Index(k);%choose the component(order of component) we want to look at 
        %figure;
        %subplot(2,4,subfig); hold on;
        icasig_comp = [];
        a0 = A(:,comp_match_run1(compIndex,1),1); %this is used as reference so that the same component across runs is always positively correlated (not inverted)
        a_runs = nan(configs.numRuns,size(A,1));
        for i=1:configs.numRuns
            if comp_match_run1(compIndex,i)>0
                a = A(:,comp_match_run1(compIndex,i),i);
                icasig_comp_one = squeeze(icasig(:,comp_match_run1(compIndex,i),i));
                if corr(a,a0)>0
                    %plot(a);
                    a_runs(i,:) = a';
                    icasig_comp(i,:) = icasig_comp_one;
                else
                    %plot(-a);
                    a_runs(i,:) = -a';
                    icasig_comp(i,:) = -icasig_comp_one;
                end
            end
        end

        clear comp;
        vector_length = size(icasig_comp,2);
        comp.avg.vectorfMRI = [];
        comp.avg.vectorfMRI = mean(icasig_comp(:,1:vector_length/2)); %avg per column (across runs)
        ub_fMRI = prctile(comp.avg.vectorfMRI,99);
        lb_fMRI = prctile(comp.avg.vectorfMRI,1);
        comp.avg.matrixfMRI = zeros(configs.numRegions,configs.numRegions);
        comp.avg.matrixfMRI(configs.mask)=comp.avg.vectorfMRI; %fill upper triangular
        comp.avg.matrixfMRI = comp.avg.matrixfMRI + (comp.avg.matrixfMRI'); %symmetrize matrix
        comp.avg.matrixfMRI = comp.avg.matrixfMRI(configs.yeoOrder,configs.yeoOrder); % apply yeoOrder
        comp.avg.vectorEEG = [];
        comp.avg.vectorEEG = mean(icasig_comp(:,vector_length/2+1:end)); %avg per column (across runs)
        comp.avg.matrixEEG = zeros(configs.numRegions,configs.numRegions);
        comp.avg.matrixEEG(configs.mask)=comp.avg.vectorEEG; %fill upper triangular
        comp.avg.matrixEEG = comp.avg.matrixEEG + (comp.avg.matrixEEG'); %symmetrize matrix
        comp.avg.matrixEEG = comp.avg.matrixEEG(configs.yeoOrder,configs.yeoOrder); % apply yeoOrder
        ub_EEG = prctile(comp.avg.vectorEEG,99);
        lb_EEG = prctile(comp.avg.vectorEEG,1);
        %% ICC evaluation of the roubst trait
        comp.avg.Aruns = a_runs;
        a_mean = nanmean(a_runs);
        comp.avg.Amean = a_mean;
        [ICC_a_fbands(k),ICC_a_subj(k)] = f_ICC_edgewise_v1(a_mean,configs);

        dispFigures = 0;
        if dispFigures == 1

            figures(end+1) = figure;
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
            figures(end+1) = figure;%(configs.h2);
            subplot(1,2,1); 
            imagesc(comp.avg.matrixfMRI,[lb_fMRI,ub_fMRI]); colormap jet; colorbar; axis square;
            set(gca,'xtick',[]); set(gca,'ytick',[]); xlabel('regions'); ylabel('regions');
            title(sprintf('connICA comp %d fMRI',compIndex))
            title(sprintf('connICA comp %d',compIndex))
            hold on
            subplot(1,2,2); 
            imagesc(comp.avg.matrixEEG,[lb_EEG,ub_EEG]); colormap jet; colorbar; axis square;
            set(gca,'xtick',[]); set(gca,'ytick',[]); xlabel('regions'); ylabel('regions');
            title(sprintf('connICA comp %d EEG',compIndex))
        end

        save(fullfile(output,['HybridT' int2str(Index(k)) '.mat']),'comp','configs','Index');
    end

    disp(['Number of selected PCs: ' num2str(numPCAComps)])

    % print figures
    % for i = 1:length(figures)
    %     filename = [output filesep 'figure' num2str(i)];
    %     print(figures(i), '-dpng', filename)
    % end
end