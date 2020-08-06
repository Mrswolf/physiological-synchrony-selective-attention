% CODE CORRESPONDING TO (STULDREHER ET AL., 2020)

% Author: Ivo Stuldreher, TNO

clearvars; close all;

flag_preprocess = 0;
flag_randclass = 0;

%% STARTUP
% load stimulus information ans participant's answers to questions
load('Stimuli'); load('answers');

% find .set datafiles in directory specified in 'filepath'
filepath = '../../data';
d = dir([filepath, '/*.set']);
isub = [d(:).isdir] ==  false;
filename = {d(isub).name}';
filename(ismember(filename,{'.','..'})) = [];
N = length(filename);

% find condition of all participants based on filename
conditionList = cellfun(@(x) strcmp(x(4), 's'), filename)';

% start eeglab
eeglab;

% let user know to check if EEGLAB is set to double precision
fprintf('Confirm EEGLAB is set to double precision. Uncheck:\n')
fprintf('''File'', -> ''Memory and other options'' -> ''If set use single precision under...''\n')
fprintf('Then, press any key to continue...\n')
pause;
fprintf('Continuing...\n')

%% PRE-PROCESSING
if flag_preprocess
    
    % pre-process eeg data (filter, remove artifactual components and remove
    % bad samples)
    [ALLEEG] = eeg_preprocessing( filepath, filename);
   
    % pre-process autonomic data (detect peaks in the ecg signal end
    % separate phasic and tonic eda components)
    [ALLAUT] = aut_preprocessing( filepath, filename);
    
end

%% ANALYSIS
if ~flag_preprocess
    ALLEEG = []; ALLAUT = [];
    
    for n = 1 : N

        % load processed data of participant n and store to ALLEEG
        EEG = pop_loadset( filename{n}, [filepath, '/eeg_processed']);       
        [ALLEEG, EEG, ~] = eeg_store(ALLEEG, EEG);
        
        % load processed data of participant n and store to ALLEEG
        AUT = pop_loadset( filename{n}, [filepath, '/aut_processed']);       
        [ALLAUT, AUT, ~] = eeg_store(ALLAUT, AUT);
    
    end
    
end

% perform analysis using eeg
[isc(:,:,:,1), p_isc(:,:,1), stats_isc(:,:,1), perf_isc(:,:,1), rand_perf(:,:,:,1), q_erp_eeg, stats_erp_eeg] = eeg_analysis( ALLEEG, Stimuli, flag_randclass);

% perform analysis using autonomic measures
[isc(:,:,:,2:3), p_isc(:,:,2:3), stats_isc(:,:,2:3), perf_isc(:,:,2:3), rand_perf(:,:,:,2:3), q_erp_aut, stats_erp_aut] = aut_analysis( ALLAUT, Stimuli, flag_randclass);

% correlate isc with answers on questions
[r_ans, p_ans, rank] = answer_correlation(squeeze(isc(:,:,1,:)), answers);

%% STATISTICAL TESTING

for mm = 1 : size(isc, 4)
  
    for gg = 1 : 2
        
        for ss = 1 : 4
            
            % test if data is normally distributed    
            [H(mm , gg, ss), pValue(mm, gg, ss), W(mm,gg, ss)] = swtest(isc(:,gg,ss,mm));

            % use t-test for normally distributed data and wilcoxon test
            % for non-normally distributed data
            if H(mm,gg,ss)
                [p_isc(ss,gg,mm), ~, stats_isc_alt(ss,gg,mm)] = signrank(isc(conditionList == gg-1, 1,ss,mm), isc(conditionList == gg-1, 2,ss,mm), 'method', 'approximate');
            else
                [~,p_isc(ss,gg,mm), ~, stats_isc(ss,gg,mm)] = ttest(isc(conditionList == gg-1, 1,ss,mm), isc(conditionList == gg-1, 2,ss,mm));
            end

        end

    end
    
end

%% VISUALIZATION

% figure 2
plot_isc(isc, p_isc, conditionList);

% figure 3
plot_erp;

%% PRINT RESULTS

% print some results
measureName = {'EEG', 'EDA', 'IBI'};
stimCondName = {'narrative & stimuli', 'beep counting task', 'affective sounds', 'narrative only'};
groupName = {'NA', 'SSA'};
for mm = 1 : size(isc, 4)
    for ss = 1 : size(isc, 3)
        fprintf('%s\n', measureName{mm})
        fprintf('Analysis considering the %s\n', stimCondName{ss})
        fprintf(' Statistical difference between within-group and between-group ISC:\n')
        for g = 1 : 2
            fprintf('   %s: %s\n    t(%i) = %f, p = %f\n', groupName{g}, string(p_isc(ss,g,mm) < .05),  stats_isc(ss,g,mm).df, stats_isc(ss,g,mm).tstat, p_isc(ss,g,mm))
        end

        fprintf(' Classification performance:\n')
        fprintf('   %f\n', mean(perf_isc(:,ss,mm)))
        fprintf('\n')
        
        % classification performance above chance
        [h(mm,ss), p(mm,ss), ~, stats(mm,ss)] = ttest2(mean(perf_isc(:,ss,mm)), mean(rand_perf(:,:,ss,mm), 2), ...
            'tail', 'right');
        
    end
end