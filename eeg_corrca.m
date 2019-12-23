
function [isc_tg, p, stats, perf, rand_perf] =  eeg_corrca( ALLEEG, stimulusIdc, flag_randclass)

N = length(ALLEEG);
T = min([ALLEEG.pnts]);
D = ALLEEG(1).nbchan;
X = zeros(T,D,N);

% number of stimulus conditions
C = length(stimulusIdc);

% add all eeg responses in one 3D-array
for n = 1 : N
    X(:,:,n) = ALLEEG(n).data(:,1:T)';
end
    
% find condition of all participants based on filename
conditionList = cellfun(@(x) strcmp(x(4), 's'), {ALLEEG.setname});

% pre-assign variables for faster processing
p = cell(C, 1); stats = cell(C, 1); perf = cell(C,1); 
rand_perf = cell(C, 1);
isc = zeros(N,2,C);

% compute synchrony toward groups
nShuffle = 100;
for i = 1 : C
    [isc_tg, perf{i}] = isc2group(X(stimulusIdc{i}, :, :), conditionList + 1);
    isc(:,:,i) = sum(isc{i}(1:3,:,:));
    
    if flag_randclass
        for j = 1 : nShuffle
            conditionListAdapt = conditionList(randperm(length(conditionList)));

            [~, rand_perf{i}(j,:)] = isc2group(X(stimulusIdc{i}, :, :), conditionListAdapt + 1);

        end
    end
    
end

%%
for c = 1 : C
%     % compute mean isc to both groups averaged over both groups of
%     % participants
%     mean_isc = zeros(2,2);
%     for g = 1 : 2
%         mean_isc(g,:) = mean(sum(isc_tg{c}(1:3,conditionList+1 == g, :)));
%     end
% 
%     fig = figure('Position', [680 666 414 312]);
%     ax = axes(fig);
%     hold on;
%     b = bar(mean_isc, ...
%         'FaceColor', 'Flat', ...
%         'LineWidth', 1.5, ...
%         'LineStyle', 'none', ...
%         'BaseValue', -1);
%     b(1).CData = [.5 .5 .5];
% 
%     b(2).CData = [.8 .8 .8];
% 
%     % plot individual participants as colored dots
%     Color = {[1 .4 .4], [.4 .4 1]; [.4 .4 1],  [1 .4 .4]};
%     Linestyle = {':', '-'; '-', ':'};
%     nBars = 2;
%     nGroups = 2;
%     groupWidth = min(0.8, nBars/(nBars + 1.5));
%     for i = 1 : nGroups
%         x(i,:) = (1:nGroups) - groupWidth/2 + (2*i-1) * groupWidth / (2*nBars);    
%         for j = 1 : nBars
%             idc = find(conditionList+1 == j);
%             for k = 1 : length(idc)     
% 
%                 [~, idx] = max(squeeze(sum(isc_tg{c}(1:3,idc(k), :))));
% 
%                 plot(x(:,j), squeeze(sum(isc_tg{c}(1:3, idc(k), :))), ...
%                     'Color', Color{(idx==j)+1, j},...
%                     'Marker', '.', ...
%                     'LineWidth', 1.5, ...
%                     'LineStyle', Linestyle{(idx==j)+1, j}, ...
%                     'MarkerFaceColor', 'none', ...
%                     'MarkerEdgeColor', Color{(idx==j)+1, j}, ...
%                     'MarkerSize', 15);
% 
%             end          
%         end 
%     end
% 
%     % set axis specifications
%     set(ax, 'YLim', [-0.005, 0.04], ...
%         'FontSize', 14, ...
%         'FontName', 'Times New Roman', ...
%         'XTick', [1, 2], ...
%         'XTickLabel', {'NA', 'SSA'});

    % statistical test
    p{c} = zeros(1,2);
    for g = 1 : 2
        [~,p{c}(g), ~, stats{c}(g)] = ttest(squeeze(sum(isc_tg{c}(1:3,conditionList+1 == g, 1))), squeeze(sum(isc_tg{c}(1:3,conditionList+1 == g, 2))));
        sigstar({x(:,g)}, p{c}(g));
    end

%     % make legend
%     h(1) = bar(NaN,NaN, 'LineStyle', 'none', 'FaceColor', [.5 .5 .5]);
%     h(2) = bar(NaN,NaN, 'LineStyle', 'none', 'FaceColor', [.8 .8 .8]);
%     legend(h, {'ISC-NA', 'ISC-SSA'}, 'Location', 'nw')
%     legend('boxoff')
end
