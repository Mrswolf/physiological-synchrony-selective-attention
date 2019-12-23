
clearvars; close all;

load('Stimuli')

flag_preprocess = 0;
flag_randclass = 0;

% find .set datafiles in directory specified in 'filepath'
filepath = '../data';
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
[isc(:,:,:,1), p_isc(:,:,1), stats_isc(:,:,1), perf_isc(:,:,1), rand_perf(:,1,:), q_erp_eeg, stats_erp_eeg] = eeg_analysis( ALLEEG, Stimuli, flag_randclass);

% perform analysis using autonomic measures
[isc(:,:,:,2:3), p_isc(:,:,2:3), stats_isc(:,:,2:3), perf_isc(:,:,2:3), rand_perf(:,2:3,:), q_erp_aut, stats_erp_aut] = aut_analysis( ALLAUT, Stimuli, flag_randclass);

% correlate isc with answers on questions
[r_ans, p_ans, rank] = answer_correlation(squeeze(isc(:,:,1,:)), answers);

%% VISUALIZATION

% figure 1
plot_isc(isc, p_isc, conditionList);

% figure 2
plot_corr(squeeze(isc(:,:,1,:)), rank, r_ans, p_ans);

% figure 3
plot_erp;