
function [h, p, stats] = eeg_erpanalysis( ALLEEG, varargin)

% set detault options
epochTime = [-.1, 1];
baselineTime = [-100, 0];
channel = 13;
sepBeeps = 0;

% variable inputs
for i = 1 : length(varargin)/2
    switch varargin{2*i-1}
        case 'epochTime', epochTime = varargin{2*i};
        case 'baselineTime', baselineTime = varargin{2*i};
        case 'channel', channel = varargin{2*i};
        case 'sepBeeps', sepBeeps = varargin{2*i};
        otherwise, error(['Unknown option ' varargin{2*i-1}])
    end
end

if strcmp({ALLEEG(1).chanlocs(channel).labels}, 'Phasic')
    for i = 1 : length(ALLEEG)
        ALLEEG(i).data(channel,:) = zscore(ALLEEG(i).data(channel,:));
    end
end

% number of participants
N = length(ALLEEG);

% find all events with id 01 and add an event every second for the 14th
% seconds afterwards (i.e., corresponding to each beep onset)
if sepBeeps
    for n = 1 : N
        event_idc = find(strcmp({ALLEEG(n).event.type}, '01'));          
        for i = 1 : length(event_idc)
            for j = 1 : 14
                e = length(ALLEEG(n).event);
                ALLEEG(n).event(e+1).type = '019';
                ALLEEG(n).event(e+1).latency = ALLEEG(n).event(event_idc(i)).latency + (2*j*ALLEEG(n).srate);
                ALLEEG(n).event(e+1).urevent = n+1;
            end
        end
        % check for consistency and reorder the events chronologically...
        ALLEEG(n) = eeg_checkset(ALLEEG(n), 'eventconsistency');
    end
end

% find event types 
combinedEvents = [ALLEEG.event];
eventTypes = unique({combinedEvents.type});
del = zeros(1, length(eventTypes));
for e = 1 : length(eventTypes)
    del(e) = length(eventTypes{e}) < 2 | strcmp(eventTypes{e}, '31') ...
        | strcmp(eventTypes{e}, '41') | strcmp(eventTypes{e}, 'boundary') ...
        | strcmp(eventTypes{e}, '21');
end
eventTypes(find(del)) = [];
eventTypes = cellfun(@(x) x(1:2), eventTypes, 'UniformOutput', false);
eventTypes = unique(eventTypes);

% define lengths for loops
T = round(diff(epochTime)*ALLEEG(1).srate);
E = length(eventTypes);
D = ALLEEG(1).nbchan;

% preallocate variables for speed
ERP = zeros(D,T,N,E); cond = zeros(N,1);

for n = 1 : N
    
    % epoch EEG of participant n
    EEG = pop_epoch( ALLEEG(n), {}, epochTime);
    
    % find indices for baselining
    baselineIdc = EEG.times >= baselineTime(1) & EEG.times <= baselineTime(2);
    
    % compute response traces averaged over trials of each event type
    for e = 1 : E
        ERP(:,:,n,e) = nanmean(EEG.data(:,:,contains({EEG.epoch.eventtype}, eventTypes{e})), 3);
        ERP(:,:,n,e) = ERP(:,:,n,e) - nanmean(ERP(:,baselineIdc,n,e), 2);
    end
    
    cond(n) = strcmp(ALLEEG(n).setname(4), 's');

end

% preallocate variables for speed
h = zeros(T,E); p = zeros(T,E);
for e = 1 : E
    for t = 1 : T
        
        % test differences of response traces between groups
        [h(t,e), p(t,e), ~, stats(t,e)] = ttest2(ERP(channel,t,cond == 0,e), ERP(channel,t,cond == 1,e));
        
    end
    
    % correct for multiple comparisons
    [h(:,e), ~, ~, p(:,e)] = fdr_bh(p(:,e), 0.05);

end

% preallocate variables for speed
SEM = zeros(D,T,2,E);

pThresholds = [1, .05, .01, .001];
pText = {'*', '**', '***'};
Color = {'Blue', 'Red'};
linestyle = {'-', '--'};

fig = figure();

pos = [.1 .66125 .8 .2875; ...
        .1 .35 .8 .2875];

% create axes for plotting responses of channnel c
ax(2) = axes(fig, ...
    'Position', [.1, .1, .8 .2], ...
    'XColor', 'Black', ...
    'XLim', EEG.times([1,end]), ...
    'FontName', 'Times New Roman', ...
    'FontSize', 12);
xlabel('Time [ms]', 'Interpreter', 'Latex');
if diff(epochTime) > 2
    xlabel('Time [s]', 'Interpreter', 'Latex');
    xticklabels(cellfun(@(x) num2str(str2num(x) * 10), xticklabels, 'UniformOutput', false));
end
ylabel({'$$-\log_{10}(p)$$'}, 'Interpreter', 'Latex');
hold(ax(2), 'all');

color = {'Blue', 'Green'};
    
for e = 1 : E
    
    if all(isnan(h(:,e)))
        continue
    end
   
    % create axes for plotting responses of channnel c
    erpax(e) = axes(fig, ...
        'Position', pos(e,:), ...
        'XColor', 'none', ...
        'XLim', EEG.times([1,end]), ...
        'FontName', 'Times New Roman', ...
        'FontSize', 12);
    hold(erpax(e), 'all');
    
    % cycle through linestyles
%     set(groot,'DefaultAxesLineStyleOrder','-|--');
    
    % fill areas where there are significant beween-group differences
    i_sig = find(h(:,e) == 1);
    yLimit = [1.2*min(min(nanmean(ERP(channel,:,cond == 0,e), 3)), min(nanmean(ERP(channel,:,cond == 1,e), 3))), ...
        1.2*max(max(nanmean(ERP(channel,:,cond == 0,e), 3)), max(nanmean(ERP(channel,:,cond == 1,e), 3)))];
    for j = 1 : length(i_sig)
        time_sig = EEG.times(i_sig(j)) - 500*(1/EEG.srate) : 1 / EEG.srate : EEG.times(i_sig(j)) + 500*(1/EEG.srate);
        fill([time_sig, fliplr(time_sig)], [-10*ones(1, length(time_sig)), 10*ones(1, length(time_sig))], 'Black', ...
            'FaceAlpha', .2, ...
            'LineStyle', 'none');
    end
    
    set(erpax(e), 'YLim', yLimit);
    
    % text of subplot letter (A, B, C, ...)
    text(erpax(e), .01,1.01, char('A' + e - 1), 'Units', 'Normalized', 'HorizontalAlignment', 'lEFT','VerticalAlignment', 'Bottom', 'FontName', 'Times New Roman', 'FontWeight', 'Bold')
        
           
    for c = 1 : 2
                
        % standard area of the mean 
        SEM(:,:,c,e) = std(ERP(:,:,cond == c-1,e), [], 3) / sqrt(length(find((cond == c-1))));
        
        % plot ERP
        plo(c) = plot(erpax(e), EEG.times, nanmean(ERP(channel,:,cond == c-1,e), 3), ...
            'Color', Color{c}, 'LineStyle', linestyle{c});
        
        % plot and area of SEM around the mean response
        TIME = [EEG.times, fliplr(EEG.times)];
        AREA = [nanmean(ERP(channel,:,cond == c-1,e), 3) - SEM(channel,:,c,e), fliplr(nanmean(ERP(channel,:,cond == c-1,e), 3) + SEM(channel,:,c,e) )];
        TIME = TIME(~isnan(AREA));
        AREA = AREA(~isnan(AREA));
        fill(TIME, AREA, Color{c}, 'FaceAlpha', .1, 'LineStyle', 'none');
                
    end
    
    plot(erpax(e), [0, 0], [-10 10], ...
        'Color', [0,0, 0], ...
        'LineStyle', '-');
    
    if e == 1
        legend(erpax(e), plo, {'NA', 'SSA'}, 'Location', 'NorthWest');
    end
    
    if e == 1
        if strcmp(EEG.chanlocs(channel).labels, 'Phasic')
            yl = ylabel('Amplitude [-]', 'Interpreter', 'Latex', 'units', 'normalized');  
        elseif strcmp(EEG.chanlocs(channel).labels, 'IBI')
            yl = ylabel('$$\Delta$$IBI [ms]', 'Interpreter', 'Latex', 'units', 'normalized');
        else
            yl = ylabel('Amplitude [$$\mu$$V]', 'Interpreter', 'Latex', 'units', 'normalized');
        end
        yl.Position(2) = -.2;
        yl.Position(1) = -.06;
    end
    
    plot(ax(2), [0, 0], ax(2).YLim, ...
    'Color', [0,0, 0], ...
    'LineStyle', '-');
    
    % create horizontal lines at p-values .05, .01 and .001 and draw 1, 2
    % and 3 stars respectively
    for pv = 1 : 3
        if min(p(:,e)) < pThresholds(pv)
            plot(ax(2), EEG.times([1, end]), [1,1]*-log10(pThresholds(pv+1)), '--', ...
                'Color', 'Black')
            text(ax(2), .2, -log10(pThresholds(pv+1))+.1, pText{pv}, ...
                'FontSize', 14, ...
                'FontName', 'Times New Roman', ...
                'HorizontalAlignment', 'Center');
        end
    end
    
    % scatter plot p-values (as a Manhattan plot)
     sc(e) = scatter(ax(2), EEG.times, -log10(p(:,e)), -70*log(p(:,e)), color{e}, ...
         'Marker', '.');
    
end

lim = [min([erpax.YLim]), max([erpax.YLim])];
[~, icons] = legend(ax(2), sc, {'A', 'B'}, 'Location', 'Southeast');
for i = 3 : 4
    icons(i).Children.MarkerSize = 12;
end
% legend(ax(2), 'boxoff');

for e = 1 : E
    erpax(e).YLim = lim;

    if strcmp(EEG.chanlocs(channel).labels, 'IBI')
        erpax(e).YTickLabels = cellfun(@(x) num2str(str2num(x) * 1000), erpax(e).YTickLabels, 'UniformOutput', false);
    end
end

savefig(fig, ['tmperp_', ALLEEG(1).chanlocs(channel).labels]);

end