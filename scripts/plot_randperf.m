close all;

load('rand_perf_271219');
perf = [25/26, 19/26, 20/26; ...
    23/26, 18/26, 15/26; ...
    19/26, 19/26, 11/26; ...
    18/26, 16/26, 19/26];

ylab = {'EEG', 'EDA', 'IBI'};
stim = {'Whole stimulus', 'Beeps', 'Affective sounds', 'Narrative only'};
color = {[1, 1, 1], [1 1 1 ], [1 1 1 ], [1 1 1]};

% setup figure parameters
mar_hor = .1;
mar_ver = .1;
spa_hor = .05;
spa_ver = .05;
nhor = 4;
nver = 3;
horwidth = (1-2*mar_hor - (nhor-1)*spa_hor) / (nhor);
verheigth = (1-2*mar_ver - (nver-1)*spa_ver) / (nver);

fig = figure;

% create background axes
background = axes(fig, ...
    'Position', [0, 0, 1, 1], ...
    'XColor', 'None', 'YColor', 'None', ...
    'Xlim', [0, 1], 'YLim', [0,1]);
hold on;

for j = 1 : size(perf, 1)
    
    xpos = mar_hor + (j-1)*(horwidth + spa_hor);
   
    rectangle('Parent', background, ...
        'Position', [xpos-spa_hor/2 0 horwidth+(spa_hor) 1], ...
        'Linestyle', 'none', 'Facecolor', color{j});
    
end

% rectangle('position', [.125, .92, horwidth-2*.025, .005], 'linestyle', 'none', 'facecolor', 'black', 'parent', background);
% plot([.126, .157, .197 ,.23], [.920025 .920025 .920025 .920025], 'r*')

cnt = 1;
for pc = 1 : size(perf,2)
    
    for j = 1 : size(perf, 1)
        
        xpos = mar_hor + (j-1)*(horwidth + spa_hor);
        ypos = (1-mar_ver) - (pc-1)*(verheigth + spa_ver) - verheigth;
        pos = [xpos, ypos, horwidth, verheigth];
        
        % create subaxis
        ax = axes('units', 'normalized', 'position', pos, ...
                'XLim', [0, 26], ...
                'YLim', [0, 40], ...
                'fontsize', 10, ...
                'fontname', 'times new roman', 'fontsize', 10, ...
                'Parent', fig, ...
                'Color', color{j});
        hold on;
%         ax = subplot(size(perf,2), size(perf,1), cnt, ...
%             'XLim', [0, 26], ...
%             'YLim', [0, 40], ...
%             'fontsize', 10, ...
%             'fontname', 'times new roman', 'fontsize', 10, ...
%             'Parent', fig);
%         hold on;
        
        % plot histogram of random performance for measure pc considering
        % stimulus condition j
        h = histogram(mean(rand_perf(:,:,j,pc)*26, 2), 'numbins', 26, 'binedges', [-.5 : 1 : 25.5]);
        
        if ttest2(perf(j,pc), mean(rand_perf(:,:,j,pc), 2), 'tail', 'right')
            linecolor = [.2 .8 .2];
            linestyle = '-';
        else
            linecolor = [.8 .1 .1];
            linestyle = ':';
        end
        
        % plot obtained classification performance
        plot(perf(j,pc)*[26, 26], ylim, 'color', linecolor, 'linewidth', 1.5, 'linestyle', linestyle)
        
        set(ax, 'YTick', [0, 15, 30], ...
            'XTick', [5, 13, 21], ...
            'XTickLabels', {'20', '50', '80'});
        
        % remove y-ticks from plots except the leftmost
        if j~= 1
            set(ax, 'YTickLabels', {});
            if j == 4
                 text(ax, 1.02, .5, [ylab{pc}], 'units', 'normalized','horizontalalignment', 'left', 'verticalalignment', 'middle', 'fontname', 'times new roman', 'fontsize', 10);
            end
        else
           if pc == 2
                yl = ylabel(ax, {'Renditions [\%]'}, 'Interpreter', 'Latex', 'units', 'normalized', 'horizontalalignment', 'center', 'rotation', 90, 'fontsize', 10);
            end
        end
        
        % remove x-ticks from plots except the bottommost
        if pc~= 3
            set(ax, 'XTickLabels', {});
            if pc == 1
                text(ax, .5, 1.02, [stim{j}], 'units', 'normalized','horizontalalignment', 'center', 'verticalalignment', 'bottom', 'fontname', 'times new roman', 'fontsize', 10);
            end
        else
            if j == 2
                xlab = xlabel(ax, {'Accuracy [\%]'}, 'Interpreter', 'Latex', 'units', 'normalized');
                xlab.Position(1) = [1.25];
            end
        end
        
        % text of subplot letter (A, B, C, ...)
%         text(ax, -.02,1.02, char('A' + cnt - 1), 'Units', 'Normalized', 'HorizontalAlignment', 'Right','VerticalAlignment', 'Bottom', 'FontName', 'Times New Roman', 'FontWeight', 'Bold')
        
        cnt = cnt + 1;
        
    end
    
end        


% 
% axall = axes(fig, ...
%     'Position', [.5, .5, 1, 1]);
        
