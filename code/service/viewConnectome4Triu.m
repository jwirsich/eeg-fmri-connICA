function [fig, label] = viewConnectome4Triu(vec, regions, label, isSaturate)
% Visualize Connectome Matrix in triu order
%
% Jonathan Wirsich, Enrico Amico 2020
%
% Wirsich, J., Amico, E., Giraud A.L. Goñi, J, Sadaghiani S.,2020 
% Multi-timescale hybrid components of the functional brain connectome: A bimodal EEG-fMRI decomposition
% Network Neuroscience (2020) 4 (3): 658–677. https://doi.org/10.1162/netn_a_00135
%
    reflect_folder = fileparts(mfilename('fullpath'));
    tmp = load(fullfile(reflect_folder, '/../../data/aparc_a2009_yeoRS7_148reg_eeg_nosubc.mat'));
    
    mask_ut = triu(true(regions,regions),1);
    conn = zeros(regions);
    
    if isSaturate
        ub = prctile(vec,99);
        lb = prctile(vec,1);
    end
    
    conn(mask_ut) = vec;
    conn = conn + conn';
    fig = figure('name', label);
    imagesc(conn(yeoOrder, yeoOrder), [lb,ub]);
    colormap('parula')
    axis square
    
    set(gca,'Yticklabel',[]) 
    set(gca,'Xticklabel',[])
    
%     set( gca, 'Units', 'normalized', 'Position', [0 0 1 1] );
%     export_fig
    set(gca,'LooseInset',get(gca,'TightInset'))
end