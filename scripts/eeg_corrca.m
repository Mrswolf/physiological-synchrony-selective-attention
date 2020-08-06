
function [isc, p, stats, perf, rand_perf] =  eeg_corrca( ALLEEG, stimulusIdc, flag_randclass)

N = length(ALLEEG);
T = min([ALLEEG.pnts]);
D = ALLEEG(1).nbchan;
conditionList = cellfun(@(x) strcmp(x(4), 's'), {ALLEEG.setname});

X = zeros(T,D,N);

% number of stimulus conditions
C = length(stimulusIdc);

% add all eeg responses in one 3D-array
for n = 1 : N
    X(:,:,n) = ALLEEG(n).data(:,1:T)';
end
    
% find condition of all participants based on filename
nShuffle = 100;

% pre-assign variables for faster processing
p = zeros(C,2); perf = zeros(2,C); 
rand_perf = zeros(nShuffle,2,C);
isc = zeros(N,2,C);
isc_tg = cell(C,1);

% compute synchrony toward groups
for i = 1 : C
    x = nan(size(X));
    x(stimulusIdc{i},:,:) = X(stimulusIdc{i},:,:);
    
    [isc_tg{i}, perf(:,i)] = isc2group(X(stimulusIdc{i}, :, :), conditionList + 1);
    isc(:,:,i) = sum(isc_tg{i}(1:3,:,:));
    
    if flag_randclass
        for j = 1 : nShuffle
            conditionListAdapt = conditionList(randperm(length(conditionList)));

            [~, rand_perf(j,:,i)] = isc2group(X(stimulusIdc{i}, :, :), conditionListAdapt + 1);

        end
    end
    
end

%%
for c = 1 : C

    % statistical test
    for g = 1 : 2
        [~,p(c,g), ~, stats(c,g)] = ttest(squeeze(sum(isc_tg{c}(1:3,conditionList+1 == g, 1))), squeeze(sum(isc_tg{c}(1:3,conditionList+1 == g, 2))));
    end

end
