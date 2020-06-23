function plot_corr(isc, rank, r_ans, p_ans)

% setup plot settings
markercolor = {[.2, .8, .2], 'Red', 'Blue'};
markertype = {'o', 'o', 'o'};
grouptocorr = [1, 2, 2, 1];

% pre-define plot labels
xlab = {'ISC-EEG', 'ISC-EDA', 'ISC-IBI'};
xlim = {[-.01, .03], [0, .4], [0, .4]};
xlimalt = {[-.025, .025], [-0.2, .2], [-.2, .2]};
ylab = {'Rank-N', 'Rank-I', 'Rank-B', '\delta Rank'};
ylimalt = [-20, 20];
leg = {'ISC-NA vs. narrative', 'ISC-SSA vs. IADS', 'ISC-SSA vs. beeps'};

% setup figure and axes
fig = figure('Position', [680 314 827 664]);

cnt = 1;
for j = 1 : size(rank, 2)
    
    for pc = 1 : size(isc, 3)   
        
    ax = subplot(size(rank, 2), size(isc, 3), cnt, ...
    'FontSize', 10, ...
    'FontName', 'Times New Roman', ...
    'YLim', [0,27], ...
    'XLim', xlim{pc});
    hold on;
    if j == 1
        set(ax, 'Color', [.9 .9 .9]);
    end
    if j == 4
        set(ax, 'XLim', xlimalt{pc}, ...
            'YLim', ylimalt, ...
            'Color', [.6 .6 .6]);
    end
        
        if p_ans(grouptocorr(j),j,pc) < .05
            linestyle = '-';
            marktype = 'o';
            markerfacecolor = 'black';
            linecolor = [.2 .8 .2];
        else
            linestyle = '-';
            marktype = 'o';
            markerfacecolor = 'none';
            linecolor = [.9 .1 .1];
        end

        % plot individuals participants as single dots
        if j == 4
            
            % fit a first order polynomial through the isc to group g for measure
            % pc vs. performance rank of answers j
            polfit = polyfit(diff(isc(~isnan(rank(:,j)),:,pc), [], 2), rank(~isnan(rank(:,j)),j), 1);
            polval = polyval(polfit, diff(isc(:,:,pc), [], 2));
            
            plot(diff(isc(:,:,pc), [], 2), rank(:,j), ...
                'Color', 'Black', ...
                'Marker', marktype, ...
                'MarkerFaceColor', markerfacecolor, ...
                'LineStyle', 'None', ...
                'MarkerSize', 5);

            % plot trend line
            pl(j) = plot(diff(isc(:,:,pc), [], 2), polval, 'linewidth', 1, 'color', linecolor, 'linestyle', linestyle, 'linewidth', 1.5);
            
        else
            
            % fit a first order polynomial through the isc to group g for measure
            % pc vs. performance rank of answers j
            polfit = polyfit(isc(~isnan(rank(:,j)),grouptocorr(j),pc), rank(~isnan(rank(:,j)),j), 1);
            polval = polyval(polfit, isc(:,grouptocorr(j),pc));
            
            plot(isc(:,grouptocorr(j),pc), rank(:,j), ...
                'Color', 'Black', ...
                'Marker', marktype, ...
                'MarkerFaceColor', markerfacecolor, ...
                'LineStyle', 'None', ...
                'MarkerSize', 5);

            % plot trend line
            pl(j) = plot(isc(:,grouptocorr(j),pc), polval, 'linewidth', 1, 'color', linecolor, 'linestyle', linestyle, 'linewidth', 1.5);
        
        end
        
        % remove y-ticks from plots except the leftmost
        if pc == 1
            ylabel(ax, ylab{j});
        else
            set(ax, 'YTickLabels', {});
        end
        
        if j == size(rank, 2)
            xlabel(ax, xlab{pc})
        else
%             set(ax, 'XTickLabels', {});
        end
        
        % text with r and p values
        strp = sprintf('%.3f', p_ans(grouptocorr(j),j,pc));
        str = sprintf('$$r =$$ %0.2f, $$p =$$ %s', r_ans(grouptocorr(j),j,pc), strp(2:end));
        text(ax, .04, .99, ['\bf{', str, '}'], 'units', 'normalized', 'verticalalignment', 'bottom', 'horizontalalignment', 'left', 'interpreter', 'latex', 'color', linecolor, 'fontweight', 'bold', 'fontsize', 10, 'fontname', 'Times New Roman'); 
        
        % text of subplot letter (A, B, C, ...)
%         text(ax, -.02,1.02, char('A' + cnt - 1), 'Units', 'Normalized', 'HorizontalAlignment', 'Right','VerticalAlignment', 'Bottom', 'FontName', 'Times New Roman', 'FontWeight', 'Bold')
                
        cnt = cnt + 1;
        
    end
       
end