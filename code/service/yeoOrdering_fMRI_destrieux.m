% cut yeo order with subcortical to order with no subc
% 
% Jonathan Wirsich, Enrico Amico 2020
%
% Wirsich, J., Amico, E., Giraud A.L. Goñi, J, Sadaghiani S.,2020 
% Multi-timescale hybrid components of the functional brain connectome: A bimodal EEG-fMRI decomposition
% Network Neuroscience (2020) 4 (3): 658–677. https://doi.org/10.1162/netn_a_00135

yeoOrder_tmp = yeoOrder;
yeoROIs_fMRI = yeoROIs;
Subc_index2 = [75:89,164];
Subc_index = [75, 82, 89, 164];
%map 76-81 + 83-88 to 1-12
for i=1:length(Subc_index)
    yeoOrder_tmp(yeoOrder_tmp==Subc_index(i))=nan;
%     yeoROIs_fMRI(yeoOrder_tmp==Subc_index(i))=nan;        
    yeoROIs_fMRI(Subc_index(i))=nan;    
end
for i=1:length(yeoOrder_tmp)
    if yeoOrder_tmp(i)<75 && ~isnan(yeoOrder_tmp(i))         
        yeoOrder_tmp(i) = yeoOrder_tmp(i)+12;
    elseif yeoOrder_tmp(i)>75 && yeoOrder_tmp(i)<83 && ~isnan(yeoOrder_tmp(i))         
        yeoOrder_tmp(i) = yeoOrder_tmp(i)-75;
    elseif yeoOrder_tmp(i)>82 && yeoOrder_tmp(i)<90 && ~isnan(yeoOrder_tmp(i))
        yeoOrder_tmp(i) = yeoOrder_tmp(i)-76;
    elseif yeoOrder_tmp(i)>89 && ~isnan(yeoOrder_tmp(i))
        yeoOrder_tmp(i) = yeoOrder_tmp(i)-3;
    else
        yeoOrder_tmp(i) = yeoOrder_tmp(i);
    end
end

yeoOrder_fMRI=yeoOrder_tmp(~isnan(yeoOrder_tmp));
yeoROIs_fMRI=yeoROIs_fMRI(~isnan(yeoROIs_fMRI));

% yeoOrder_fMRI = yeoOrder_fMRI([75:86, 1:74, 87:160]);
yeoROIs_fMRI = yeoROIs_fMRI([75:86, 1:74, 87:160]);