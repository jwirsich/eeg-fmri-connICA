% cut yeo order with subcortical to order with no subc
% 
% Jonathan Wirsich, Enrico Amico 2020
%
% Wirsich, J., Amico, E., Giraud A.L. Goñi, J, Sadaghiani S.,2020 
% Multi-timescale hybrid components of the functional brain connectome: A bimodal EEG-fMRI decomposition
% Network Neuroscience (2020) 4 (3): 658–677. https://doi.org/10.1162/netn_a_00135
%
yeoOrder_tmp = yeoOrder;
yeoROIs_eeg = yeoROIs;
Subc_index = [35:49,84];
for i=1:length(Subc_index)
    yeoOrder_tmp(yeoOrder_tmp==Subc_index(i))=nan;
    yeoROIs_eeg(yeoOrder_tmp==Subc_index(i))=nan;    
end
for i=1:length(yeoOrder_tmp)
    if yeoOrder_tmp(i)>49 && ~isnan(yeoOrder_tmp(i))         
        yeoOrder_tmp(i) = yeoOrder_tmp(i)-15;
    else
        yeoOrder_tmp(i) = yeoOrder_tmp(i);
    end
end

yeoOrder_eeg=yeoOrder_tmp(~isnan(yeoOrder_tmp));
yeoROIs_eeg=yeoROIs_eeg(~isnan(yeoROIs_eeg));
