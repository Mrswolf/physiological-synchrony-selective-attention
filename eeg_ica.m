function [icaweights, icasphere, artcomps] = eeg_ica(EEG)

ALLEEG = [];
CURRENTSET = 0;

% remove bad channels using trimOutlier
EEG_PRE = EEG;
channelSdLowerBound = 1;    % remove channels whose SD lower bound < 2 uV
channelSdUpperBound = 50; %50  % remove channels whose SD upper bound > 250 uV
amplitudeThreshold = 100; %150  % remove 100 ms around aplitudes larger than 250 uV
pointSpreadWidth = 100;
EEG = trimOutlier(EEG, channelSdLowerBound, channelSdUpperBound, amplitudeThreshold, pointSpreadWidth);

% check how many channels are left
if EEG.nbchan < 27
   warning(['Only %i correct channels left'], EEG.nbchan) 
end

% interpolate bad channels
EEG = pop_interp(EEG, EEG_PRE.chanlocs, 'spherical');

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