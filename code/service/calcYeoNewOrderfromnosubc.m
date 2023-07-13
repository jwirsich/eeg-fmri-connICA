% Desikan atlas: Put central eyefields into attention network (from orignial yeo7mapping)
% Jonathan Wirsich, Enrico Amico 2020
%
% Wirsich, J., Amico, E., Giraud A.L. Goñi, J, Sadaghiani S.,2020 
% Multi-timescale hybrid components of the functional brain connectome: A bimodal EEG-fMRI decomposition
% Network Neuroscience (2020) 4 (3): 658–677. https://doi.org/10.1162/netn_a_00135

%left caudal middle frontal = position 3
%right caudal middle frontal = position 37

temp_yeoRois_eeg = yeoROIs_eeg;

temp_yeoRois_eeg(3) = 3;
temp_yeoRois_eeg(37) = 3;

%find max 3
tempOrder = yeoOrder_eeg;
yeoOrder_eeg();

%attent network = position 22+23 (28+62)

%put l caudal on 22, 22 to 23, put r cuadal on 24 and 22 to 25
%shift all the rest +1

tempOrder = tempOrder([1:21 49 22 44 23:43 45:48 50:68]);

test = temp_yeoRois_eeg(tempOrder);
test2 = yeoROIs_eeg(yeoOrder_eeg);

yeoROIs_eeg = temp_yeoRois_eeg;
yeoOrder_eeg = tempOrder;