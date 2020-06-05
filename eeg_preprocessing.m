function [ALLEEG] = eeg_preprocessing( filepath, filename)

N = length(filename);
ALLEEG = [];

for n = 1 : N
    
    % load data of participant n
    EEG = pop_loadset( filename{n}, filepath);
        
    % Select only EEG channels
    EEG = pop_select( EEG, ...
        'channel', 1:32);
    
    % highpass filter data at 1 Hz
    EEG = pop_eegfiltnew( EEG, 1, 0);
    
    % notch filter data around 50 Hz
    EEG = pop_eegfiltnew( EEG, 45, 55, [], 1);
    
    % apply average reference after adding initial reference as zeros
    EEG.nbchan = EEG.nbchan+1;
    EEG.data(end+1, :) = zeros(1, EEG.pnts);
    EEG.chanlocs(1, EEG.nbchan).labels = 'initialReference';
    EEG = pop_reref(EEG, []);
    EEG = pop_select(EEG, ...
        'nochannel', {'initialReference'});
    
    % run independent component analysis (ICA)
    dataRank = rank(double(EEG.data'));
    [EEG.icaweights, EEG.icasphere] = runica(EEG.data, ...
        'pca', dataRank, ...
        'extended', 0);
    EEG = eeg_checkset(EEG, 'ica');
    
    % classify artifactual components using MARA
    options = [0, 0, 0, 0, 0];
    [~, EEG, ~] = processMARA(ALLEEG, EEG, CURRENTSET, options);
    artifact_comps = find(EEG.reject.gcompreject);
     
    % remove artifactual components

    EEG = eeg_checkset(EEG, 'ica');
    EEG = pop_subcomp(EEG, artifact_comps, 0);

    % remove bad samples based of which power is 4 std above or below the
    % mean in an iterative way with 4 repititions
    for i = 1 : EEG.nbchan
        for j = 1 : 4
            sqr_amp = EEG.data(i,:).^2;
            m_sqr_amp = mean(sqr_amp);
            std_sqr_amp = std(sqr_amp);

            idc = find(sqr_amp > (m_sqr_amp + 4*std_sqr_amp) | sqr_amp < (m_sqr_amp - 4*std_sqr_amp));

            EEG.data(i,idc) = nan;   
        end
    end

    % save processed dataset
    [~, name] = fileparts(filename{n});
    EEG = pop_saveset( EEG, ...
        'filename', [name, '.set'], ...
        'filepath', [filepath, '/', 'eeg_processed']);
    
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);

end

end