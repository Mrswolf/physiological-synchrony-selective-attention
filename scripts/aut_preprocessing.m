function [ALLEEG] = aut_preprocessing( filepath, filename)

N = length(filename);
ALLEEG = [];

if exist([filepath, '/', 'aut_processed']) ~= 7
    mkdir([filepath, '/', 'aut_processed']);
end

for n = 1 : N
    
    % load data of participant n
    EEG = pop_loadset( filename{n}, filepath);
        
    % Select only EDA and ECG channels
    EEG = pop_select( EEG, ...
        'channel', {'EXG1', 'GSR1'});
    
    % ecg r-peak detection
    [~, qrs_i, delay] = ECG_pan_tompkin( EEG.data(1,:), EEG.srate, 0);
    
    % extract interbeat interval
    ibi = diff(qrs_i / EEG.srate);
    time = (qrs_i(1) - delay)/EEG.srate + cumsum(ibi);
    time = [EEG.xmin, time, EEG.xmax];
    ibi = [nan, ibi, nan];
	[ibi, time] = resample(ibi, time, 32);
    
    % resample to 32 Hz
    EEG = pop_resample( EEG, 32);
       
    % compute phasic response using Ledalab for Matlab
    [phasic] = preproeda(EEG.data(2,:), EEG.srate);
    EEG.nbchan = EEG.nbchan+1;
    if length(phasic) < length(EEG.data)
        phasic(end+1) = 0;
    end
    EEG.data(end+1, :) = phasic;
    EEG.chanlocs(1, EEG.nbchan).labels = 'Phasic';
    
    % add ibi to data vector
    EEG.nbchan = EEG.nbchan+1;
    EEG.data(end+1, :) = ibi(1:EEG.pnts);
    EEG.chanlocs(1, EEG.nbchan).labels = 'IBI';
        
    % save processed dataset
    [~, name] = fileparts(filename{n});
    
    EEG = pop_saveset( EEG, ...
        'filename', [name, '.set'], ...
        'filepath', [filepath, '/aut_processed']);
    
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);
        
end

end