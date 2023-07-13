function saveFigures(outputpath, figures, labels)
% SAVEFIGURES to PNG format
%
% Jonathan Wirsich, Enrico Amico 2020
%
% Wirsich, J., Amico, E., Giraud A.L. Goñi, J, Sadaghiani S.,2020 
% Multi-timescale hybrid components of the functional brain connectome: A bimodal EEG-fMRI decomposition
% Network Neuroscience (2020) 4 (3): 658–677. https://doi.org/10.1162/netn_a_00135
%
for i = 1:length(figures)
    filename = [outputpath filesep 'figure' num2str(i) '-' labels{i}];
    print(figures(i), '-dpng', filename)
end

end

