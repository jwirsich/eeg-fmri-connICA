%% [comp_match_run1, freq] = run_robust_connICA_fast(A,icasig,configs)
% Extract the most robust ICA traits across runs, as explained in Amico et
% al., Neuroimage 2017.
% In summary, for the trait to be considered as robust it has to correlate higher than a certain threshold across runs both in terms of trait and weights (defined by user in configs.corrMin and configs.corrMin_A); 
% also it has to appear at least in a predefined of the total runs (set by user in configs.minFreq).
% Inputs:
% -A : 3D vector of ICA weights
% -icasig : 3D vector of ICA traits
% -configs: structure of params.
    % - configs.numOfIC: number of ICA components
    % - configs.numRuns: numnber of ICA runs
    % - configs.corrMin: number of minimum correlation between single-run vectors of traits for robustness selection  
    % - configs.corrMin_A = number of minimum correlation between single-run  vectors of weights for robustness selection
    % - configs.minFreq = number of minimum frequency across runs for robustness selection
% 
% -Outputs
% The robust traits are ordered and (if necessary) rearranged according to the first (best) run
% -comp_match_run1 = vector of the indices of the robust traits across runs 
% -freq = frequency of occurence of traits per runs
%
% example of use:
% [comp_match_run1, freq] = run_robust_connICA_fast(A,icasig,configs)
%% PLEASE CITE US!
%% Enrico Amico & Joaquin Goni, Purdue University
% If you are using this code for your research, please kindly cite us:
% Mapping the functional connectome traits of levels of consciousness.
% Amico, E. ... & Goñi, J. (2017). NeuroImage, 148, 201-211.
% http://www.sciencedirect.com/science/article/pii/S1053811917300204
% and
% Version used in Wirsich, J., Amico, E., Giraud A.L. Goñi, J, Sadaghiani S.,2020 
% Multi-timescale hybrid components of the functional brain connectome: A bimodal EEG-fMRI decomposition, Netw Neurosci (just accepted), https://doi.org/10.1162/netn_a_00135
% for an updated version try: https://engineering.purdue.edu/ConnplexityLab/publications/connICA_toolbox_v1.1.tar.gz
function [comp_match_run1, freq] = run_robust_connICA_fast(A,icasig,configs)

% reshape into 2D vecs 
A_2D = single(reshape(A,size(A,1),configs.numOfIC*configs.numRuns));
icasig_2D = single(reshape(icasig,size(icasig,1),configs.numOfIC*configs.numRuns));
% Compute correlation all vs all
disp ('Computing corr btw comps (might take some time)..')
C_A = single(corr(A_2D,A_2D));
C_icasig = single(corr(icasig_2D,icasig_2D));
mask_freq = (abs(C_A) > configs.corrMin_A) & (abs(C_icasig) > configs.corrMin); % if weights AND icasig are > min Corr value
% reshape freq back to (numComp, numRuns)
freq = reshape(sum(mask_freq),configs.numOfIC,configs.numRuns);
freq = freq./configs.numRuns;

robust_comps_runs = sum(freq>configs.minFreq);
reference_run = find(robust_comps_runs==max(robust_comps_runs),1);
if reference_run ~=1 % change reference run
    Aaux = A(:,:,1);
    A(:,:,1) = A(:,:,reference_run);
    A(:,:,reference_run) = Aaux;
    icasig_aux = icasig(:,:,1);
    icasig(:,:,1) = icasig(:,:,reference_run);
    icasig(:,:,reference_run) = icasig_aux;
    
    % Update corr mats
    disp ('New ref run: updating (might take some time)..')
    A_2D = single(reshape(A,size(A,1),configs.numOfIC*configs.numRuns));
    icasig_2D = single(reshape(icasig,size(icasig,1),configs.numOfIC*configs.numRuns));
    C_A = single(corr(A_2D,A_2D));
    C_icasig = single(corr(icasig_2D,icasig_2D));
    % update freq mat
    mask_freq = (abs(C_A) > configs.corrMin_A) & (abs(C_icasig) > configs.corrMin); % if weights AND icasig are > min Corr value
    % reshape freq back to (numComp, numRuns)
    freq = reshape(sum(mask_freq),configs.numOfIC,configs.numRuns);
    freq = freq./configs.numRuns;
    % compute comp_match_run1
    comp_match_run1 = zeros(configs.numOfIC,configs.numRuns);
    comp_match_run1(:,1) = 1:configs.numOfIC;
    for i=1:configs.numOfIC
        for k=2:configs.numRuns
            range = ((k-1)*configs.numOfIC)+1:((k-1)*configs.numOfIC)+configs.numOfIC;
            c = abs(C_A(i,range));
            ctrait = abs(C_icasig(i,range));
            if max(c) > configs.corrMin_A && max(ctrait) > configs.corrMin
                pos = find(ctrait==max(ctrait));
                comp_match_run1(i,k) = pos;
            end
        end
    end    
else
    % compute comp_match_run1
    comp_match_run1 = zeros(configs.numOfIC,configs.numRuns);
    comp_match_run1(:,1) = 1:configs.numOfIC;
    for i=1:configs.numOfIC
        for k=2:configs.numRuns
            range = ((k-1)*configs.numOfIC)+1:((k-1)*configs.numOfIC)+configs.numOfIC;
            c = abs(C_A(i,range));
            ctrait = abs(C_icasig(i,range));
            if max(c) > configs.corrMin_A && max(ctrait) > configs.corrMin
                pos = find(ctrait==max(ctrait));
                comp_match_run1(i,k) = pos;
            end
        end
    end    
end
