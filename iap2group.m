function [isc_tg, perf, rand_perf, p, stats, r_ans, p_ans, median_rand_perf, confInterv] = iap2group(data, G, stimulusIdc, flag_randperf)
% IAP2GROUP computes inter-subject correlations of interpersonal autonomic
% physiology (IAP) with respect to two groups with separate
% attentional focus.
%   [ISC_TG] = IAP2GROUP(DATA, G) computes the correlations of all
%   participants with the responses of two group of participants. Input
%   (data) should be an N-dimensional struct (where N is the number of
%   participants) with fieldnames corresponding to autonomic measures
%   (e.g., phasic and ibi) . Those substructs should contain fields
%   'samplerate', 'data' and 'time'. G should be an N-dimensional vector
%   with values 1 & 2 corresponding to the attenional group. This function
%   plots results in a bar plot, showing mean subject-to-group within and
%   between-groups as well as the relation between the separate scores of
%   an individual.
%
%   Dependencies: MWA
%
%   Copyright: Ivo Stuldreher, ivo.stuldreher@tno.nl
%
%   Last saved: 08-07-2019

N = length(data);
T = min([data.pnts]);
D = data(1).nbchan;
X = zeros(T,D,N);

% add all eeg responses in one 3D-array
for n = 1 : N
    X(:,:,n) = data(n).data(:,1:T)';
end


% compute isc for every possible combination of participants p1 and p2
isc_pp = cell(D,1);
isc_tg = zeros(N,2,4,D);

nShuffle = 100;
perf = zeros(2, length(stimulusIdc), length(D));
rand_perf = zeros(nShuffle, 2, length(stimulusIdc), length(D));

for c = 1 : D
    isc_pp{c} = zeros(N,N,length(stimulusIdc));
    for p1 = 1 : N
        for p2 = p1 : N
            [~, isc_tt{c}(:,p1,p2)] = mwa(X(:,c,p1), X(:,c,p2), data(1).srate, ...
                'CorWindow', 15, ...
                'CorStep', 1);
        end
    end
    
    for i = 1 : length(stimulusIdc)
        for p1 = 1 : N
            for p2 = 1 : N
                sampLarIdc = find(isc_tt{c}(:,p1,p2) > 0);
                sampLarIdc = sampLarIdc(ismember(sampLarIdc, stimulusIdc{i}));
                sampSmaIdc = find(isc_tt{c}(:,p1,p2) < 0);
                sampSmaIdc = sampSmaIdc(ismember(sampSmaIdc, stimulusIdc{i}));
                isc_pp{c}(p1,p2,i) = log(sum(isc_tt{c}(sampLarIdc,p1,p2)) / sum(abs(isc_tt{c}(sampSmaIdc,p1,p2))));
            end
        end
        
        isc_pp{c}(:,:,i) = triu(isc_pp{c}(:,:,i))'+triu(isc_pp{c}(:,:,i));    % mirror matrix along diagonal

        % compute participant to group PS isc for
        for p1 = 1 : N
            for g = 1 : 2
                isc_tg(p1,g,i,c) = nanmean(isc_pp{c}(p1, setdiff(find(G==g), p1), i));
            end
        end
        
        for g = 1 : 2
            perf(i,g,c) = (2/N) * length(find(isc_tg{c}(G == g, g,i) > isc_tg{c}(G == g, setdiff(1:2,g),i))); 
            
            [~,p(i,g,c), ~, stats(i,g,c)] = ttest(isc_tg{c}(G == g, 1,i), isc_tg{c}(G == g, 2,i));
            
        end
        
        if flag_randperf
            for s = 1 : nShuffle
                Gr = G(randperm(length(G)));
                isc_tg_rand = zeros(N, 2);
                for g = 1 : 2
                    for p1 = 1 : N
                        isc_tg_rand(p1, g) = nanmean(isc_pp{c}(p1, setdiff(find(G == g), p1), i));
                    end
                end

                for g = 1 : 2
                    rand_perf(s,g,i,c) = (2/N) * length(find(isc_tg_rand(Gr == g, g) > isc_tg_rand(Gr == g, setdiff(1:2,g))));

                end
            end
        end

    end

end

end