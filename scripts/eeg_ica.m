function [icaweights, icasphere, artcomps] = eeg_ica(EEG)

ALLEEG = [];
CURRENTSET = 0;

% apply average reference after adding initial reference as zeros
EEG.nbchan = EEG.nbchan+1;
EEG.data(end+1, :) = zeros(1, EEG.pnts);
EEG.chanlocs(1, EEG.nbchan).labels = 'initialReference';
EEG = pop_reref(EEG, []);
EEG = pop_select(EEG, ...
    'nochannel',{'initialReference'});

% run independent component analysis (ICA)
dataRank = rank(double(EEG.data'));
[EEG.icaweights, EEG.icasphere] = runica(EEG.data, ...
    'pca', dataRank, ...
    'extended', 0);
EEG = eeg_checkset(EEG, 'ica');
icaweights = EEG.icaweights;
icasphere = EEG.icasphere;

% classify artifactual components using MARA
options = [0, 0, 0, 0, 0];
[~, EEG, ~] = processMARA(ALLEEG, EEG, CURRENTSET, options);

artcomps = find(EEG.reject.gcompreject);

end