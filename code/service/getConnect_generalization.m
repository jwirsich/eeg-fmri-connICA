function conn = getConnect_generalization(sess, modal)
% Get connectivity from specified modality
%
% Jonathan Wirsich, Enrico Amico 2020
%
% Wirsich, J., Amico, E., Giraud A.L. Goñi, J, Sadaghiani S.,2020 
% Multi-timescale hybrid components of the functional brain connectome: A bimodal EEG-fMRI decomposition
% Network Neuroscience (2020) 4 (3): 658–677. https://doi.org/10.1162/netn_a_00135
    
    if strcmp(modal, 'fMRI')
        conn = sess(2, :);
    elseif strcmp(modal, 'delta')
        conn = sess(3, :);
    elseif strcmp(modal, 'theta')
        conn = sess(4, :);
    elseif strcmp(modal, 'alpha')
        conn = sess(5, :);
    elseif strcmp(modal, 'beta')
        conn = sess(6, :);
    elseif strcmp(modal, 'gamma')
        conn = sess(7, :);
    end
        
end