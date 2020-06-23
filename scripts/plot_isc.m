function plot_isc(isc, p, conditionList)

[~, G, S, M] = size(isc);

% figure parameters
margin_h = .0875;
margin_v = 0;
spa_h = .025;
spa_v = .05;
width_h = (1-2*margin_h - (S-1)*spa_h) / (S);
heigth_v = (1-2*margin_v - (M-1)*spa_v) / (M);

fig = figure('Position', [680 249 827 729]);

% draw grey bars behind each figure column
ax_b = axes(fig, ...
    'Position', [0, 0, 1, 1], ...
    'XColor', 'None', 'YColor', 'None', ...
    'Xlim', [0, 1], 'YLim', [0,1], 'Color', [1 1 1]);
hold on;

for ss = 1 : 4
    
    xpos = margin_h + (ss-1)*(width_h + spa_h);
    rectangle('Parent', ax_b, ...
        'Position', [xpos margin_v width_h 1-2*margin_v], ...
        'Linestyle', 'none', 'Facecolor', [.9 .9 .9]);
    
end

% figure parameters
margin_h = .1;
margin_bottom = .1;
margin_top = .15;
spa_h = .05;
spa_v = .05;
width_h = (1-2*margin_h - (S-1)*spa_h) / (S);
heigth_v = (1-(margin_bottom + margin_top) - (M-1)*spa_v) / (M);

groupWidth = min(0.8, G/(G + 1.5));
color = {[1 .4 .4], [.4 .4 1]};
linestyle = {':', '-'};

ylab = {'ISC-EEG', 'ISC-EDA', 'ISC-IBI'};
yLim = {[-.01, .055], [-.1, .79], [-.1, .69]};

mean_isc = zeros(G,G,M,S);
for mm = 1 : M
    
    for ss = 1 : S
        
        % set axes parameters for figure with measure mm for condition ss
        xpos = margin_h + (ss-1)*(width_h + spa_h);
        ypos = (1-margin_top) - (mm-1)*(heigth_v + spa_v) - heigth_v;
        pos = [xpos, ypos, width_h, heigth_v];
        
        ax(mm,ss) = axes('units', 'normalized', 'position', pos, ...
            'FontSize', 10, 'FontName', 'Times New Roman', ...'XLim', [.69, 2.31], ...
            'XTickLabels', {}, 'XTick', [1 , 2], ...
            'YLim', yLim{mm}, ...
            'Box', 'Off', ...
            'Color', [.9 .9 .9]);
        hold on;
        
        for g1 = 1 : G
           
            mean_isc(:,g1,ss,mm) = nanmean(isc(conditionList == g1-1,:,ss,mm));
                         
            x(g1,:) = (1:G) - groupWidth/2 + (2*g1-1) * groupWidth / (2*G);  
            
        end
        
        b = bar(ax(mm,ss), mean_isc(:,:,ss,mm), ...
            'FaceColor', 'Flat', ...
            'LineStyle', 'none', ...
            'BaseValue', 0);    
        
        for g1 = 1 : G
           
            % color within-group correlations dark and between-group
            % correlations ligth grey
            b(g1).CData(g1,:) = [.5 .5 .5];
            b(g1).CData(setdiff(1:2, g1), :) = [.8 .8 .8];
            
            % remove baseline at y = 0
            b(g1).BaseLine.Visible = 'off';
            
                idc = find(conditionList+1 == g1);
                
                for ii = idc
                    
                    [~, g_ii] = max(isc(ii,:,ss,mm));
                    
                    plot(ax(mm,ss), x(:,g1), isc(ii,:,ss,mm), ...
                        'Color', color{(g_ii == g1)+1},...
                        'Marker', '.', ...
                        'LineWidth', 1.5, ...
                        'LineStyle', linestyle{(g_ii == g1)+1}, ...
                        'MarkerFaceColor', 'none', ...
                        'MarkerEdgeColor', color{(g_ii == g1)+1}, ...
                        'MarkerSize', 15);
                    
                end
            
            sigstar({x(:,g1)}, p(ss,g1,mm));
            
        end
        
        if ss == 1
            
            ylabel(ax(mm,ss), ylab{mm});
            
            if mm == 1
                
                legend(ax(mm,ss), {'within-group', 'between-group'}, 'Location', 'Northeast', 'FontSize', 10);
                legend(ax(mm,ss), 'boxoff');
                
            end
            
        else
           
            set(ax(mm,ss), 'YTickLabel', {});
            
        end
        
        if mm == 1
            
            % text of subplot letter (a, b, c, ...)
            text(ax(mm,ss), .5,1.55, ['(', char('a' + ss - 1), ')'], 'Units', 'Normalized', 'HorizontalAlignment', 'Center','VerticalAlignment', 'Bottom', 'FontName', 'Times New Roman', 'FontWeight', 'Bold', 'Fontsize', 12)
        
        elseif mm == 3
            
            set(ax(mm,ss), 'XTickLabels', {'NA', 'SSA'});
            
        end            
        
    end
    
end