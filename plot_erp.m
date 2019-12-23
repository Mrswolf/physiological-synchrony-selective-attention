function plot_erp

ylab = {'EEG [$$\mu$$V]', 'EDA [-]', 'IBI [ms]'};
measureName = {'-EEG', '-EDA', '-IBI'};
stimName = {'WMT', 'IADS'};
ylimit = {[-1.5, 1.8], [-.25, .5], [-.03 .04]};
xlimit = {[-100, 1050], [-3000, 31500], [-3000, 31500]};

xticks = {[0 250 500 750 1000], [0 10000 20000 30000], [0 10000 20000 30000]};
xticklabels = {{'0', '0.25', '0.5', '0.75', '1'}, {'0', '10', '20', '30'}, {'0', '10', '20', '30'}};

overviewFig = figure('Position', [680 314 827 664]);
cnt = 1;

% setup figure parameters
mar_hor = .1;
mar_ver = .1;
spa_hor = .05;
spa_ver = .05;
nhor = 2;
nver = 3;
horwidth = (1-2*mar_hor - (nhor-1)*spa_hor) / (nhor);
verheigth = (1-2*mar_ver - (nver-1)*spa_ver) / (nver);

% create background axes
background = axes(overviewFig, ...
    'Position', [0, 0, 1, 1], ...
    'XColor', 'None', 'YColor', 'None', ...
    'Xlim', [0, 1], 'YLim', [0,1]);
hold on;

for s = 1 : 2
    
    xpos = mar_hor + (s-1)*(horwidth + spa_hor);
   
    rectangle('Parent', background, ...
        'Position', [xpos .1 horwidth .8], ...
        'Linestyle', 'none', 'Facecolor', [.9 .9 .9]);
    
end

for m = 1 : 3
    
    filename = ['ERP', measureName{m}, '.fig'];
    fig = openfig(filename);
    axesPlot = findobj(fig.Children, 'Type', 'Axes');
    
    for s = 1 : 2
        
        xpos = mar_hor + (s-1)*(horwidth + spa_hor);
        ypos = (1-mar_ver) - (m-1)*(verheigth + spa_ver) - verheigth;
        pos = [xpos, ypos, horwidth, verheigth];
        
        ax_new(m,s) = axes('Position', pos, 'Units', 'normalized', 'Parent', overviewFig, ...
            'FontSize', 12, 'FontName', 'Times New Roman', ...
            'xtick', xticks{m}, 'xticklabels', xticklabels{m}, ...
            'xlim', xlimit{m}, ...
            'Box', 'Off', ...
            'Color', [.9 .9 .9]);
        hold on;
        
        ax = axesPlot(3-s);
        
        copyobj(get(ax, 'Children'), ax_new(m,s));
        
        if s == 1
            ylabel(ax_new(m,s), ylab{m}, 'Interpreter', 'Latex');
            if m == 3
                ax_new(m,s).YTickLabels = {'-20', '0', '20', '40'};
            end
        else
            ax_new(m,s).YTickLabels = {};
        end
        
        findobj(ax_new(m,s).Children, 'Type', 'Text');
        txt = findobj(ax_new(m,s).Children, 'Type', 'Text');
        if m ~= 1    
            txt.String = '';
        else
            txt.String = ['(', char('a' + s - 1), ')'];
            txt.Position = [-0.01, 1.01, 0];
            txt.HorizontalAlignment = 'Right';
            txt.FontSize = 16;
        end
                
        cnt = cnt + 1;
        
        ax_new(m,s).YLim = ylimit{m};
    end
    
    close(fig);
    
end

text(background, .5, .04, 'Time [s]', 'Units', 'Normalized', 'HorizontalAlignment', 'Center','VerticalAlignment', 'Middle', 'FontName', 'Times New Roman', 'Fontsize', 14)

end
        