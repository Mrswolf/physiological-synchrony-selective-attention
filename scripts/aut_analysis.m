function [isc, p_isc, stats_isc, perf_isc, rand_perf, q_erp, stats_erp] = aut_analysis( ALLEEG, Stimuli, flag_randclass)

N = length(ALLEEG);

for n = 1 : N
	ALLEEG(n) = pop_select(ALLEEG(n), 'channel', {'Phasic', 'IBI'});
end

% find condition of all participants based on filename
conditionList = cellfun(@(x) strcmp(x(4), 's'), {ALLEEG.setname});

% create indices of four different datasets: the whole experiment excluding
% baseline, beep tasks only, affective sounds only and the no-stimuli parts
mwa_fs = 1; % moving window step size [s]  
stimulusIdc = cell(4,1);
stimulusIdc{1} = floor(ALLEEG(1).event(4).latency * mwa_fs/ALLEEG(1).srate) + 60*mwa_fs : ceil(ALLEEG(1).event(end-2).latency* mwa_fs/ALLEEG(1).srate);
for stim = 1 : 2
    stimIdc = find([Stimuli.id] == stim-1);
    extractIdc = mwa_fs*[[Stimuli(stimIdc).timestamp]', [Stimuli(stimIdc).timestamp]' + [Stimuli(stimIdc).duration]'];
    stimulusIdc{stim+1} = [];
    for i = 1 : length(extractIdc)
        stimulusIdc{stim+1} = [stimulusIdc{stim+1}, round(extractIdc(i,1)) : round(extractIdc(i,2))];
    end    
end
stimulusIdc{4} = stimulusIdc{1}(~ismember(stimulusIdc{1}, stimulusIdc{2}) & ~ismember(stimulusIdc{1}, stimulusIdc{3}));

for i = 1 : 2
       
    [~, q_erp, stats_erp] = eeg_erpanalysis( ALLEEG, 'epochTime', [-1, 30], ...
        'baselineTime', [-1000, 0], ...
        'channel', i);
end

[isc, perf_isc, rand_perf, p_isc, stats_isc] = iap2group(ALLEEG, conditionList + 1, stimulusIdc, flag_randclass);

% print some results
stimCondName = {'narrative & stimuli', 'beep counting task', 'affective sounds', 'narrative only'};
groupName = {'NA', 'SSA'};
for j = 1 : ALLEEG(1).nbchan
    for i = 1 : length(stimulusIdc)
        fprintf('Analysis considering the %s\n', stimCondName{i})
        fprintf(' Statistical difference between within-group and between-group ISC:\n')
        for g = 1 : 2
            fprintf('   %s: %s\n    t(%i) = %f, p = %f\n', groupName{g}, string(p_isc(i,g,j) < .05),  stats_isc(i,g,j).df, stats_isc(i,g,j).tstat, p_isc(i,g,j))
        end

        fprintf(' Classification performance:\n')
        fprintf('   %f\n', mean(perf_isc(:,i,j)))
        fprintf('\n')
    end
end

end