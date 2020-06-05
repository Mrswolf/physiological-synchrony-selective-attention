function [isc, p_isc, stats_isc, perf_isc, rand_perf, q_erp, stats_erp] = eeg_analysis( ALLEEG, Stimuli, flag_randclass)

% create indices of four different datasets: the whole experiment excluding
% baseline, beep tasks only, affective sounds only and the no-stimuli parts
stimulusIdc = cell(4,1);
stimulusIdc{1} = floor(ALLEEG(1).event(4).latency) + 60 * ALLEEG(1).srate : ceil(ALLEEG(1).event(end-2).latency);
for stim = 1 : 2
    stimIdc = find([Stimuli.id] == stim-1);
    extractIdc = ALLEEG(1).srate*[[Stimuli(stimIdc).timestamp]', [Stimuli(stimIdc).timestamp]' + [Stimuli(stimIdc).duration]'];
    stimulusIdc{stim+1} = [];
    for i = 1 : length(extractIdc)
        stimulusIdc{stim+1} = [stimulusIdc{stim+1}, round(extractIdc(i,1)) : round(extractIdc(i,2))];
    end    
end
stimulusIdc{4} = stimulusIdc{1}(~ismember(stimulusIdc{1}, stimulusIdc{2}) & ~ismember(stimulusIdc{1}, stimulusIdc{3}));

% compute event-related potentials with respect to the stimulid
[~, q_erp, stats_erp] = eeg_erpanalysis( ALLEEG, 'sepBeeps', 1);

% compute inter-subject correlations with respect to both groups
[isc, p_isc, stats_isc, perf_isc, rand_perf] = eeg_corrca( ALLEEG, stimulusIdc, flag_randclass);

end
