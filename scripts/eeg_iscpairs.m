
% obtain number of channels (D) and number of participants (N)
[~,D,N] = size(X);

% compute cross covariance between all participant pairs
R = permute(reshape(nancov(X(:,:)),[D N  D N]),[1 3 2 4]);

v = corrca(X, 'shrinkage', .5);
%%
% pre-assign variables for faster processing
isc = zeros(N,N);
for n1 = 1 : N
    
    for n2 = n1+1 : N
        
        Rw = R(:,:,n1,n1) + R(:,:,n2,n2);
        Rb = R(:,:,n1,n2) + R(:,:,n2,n1);
        
        C = diag(v'*Rb*v)./diag(v'*Rw*v);
        isc(n1,n2) = sum(C(1:3));
        
    end
    
end

% mirror matrix along diagonal to obtain a symmetrical matrix
isc = (isc + isc') - eye(size(isc,1)).*diag(isc);

% perform pca
[coefs, score, ~, ~, explained] = pca(isc);

%%
% surf(isc, 'EdgeColor', 'none', 'FaceColor', [.5 .5 .5]);
hold on;
color = {'Blue', 'Red'};

myColors = zeros(N, 3);
rowsBlue = conditionList == 0;
rowsRed = conditionList == 1;

myColors(rowsBlue, 3) = 1;
myColors(rowsRed, 1) = 1;

for n1 = 1 : N
    
    subplot(6, 5, n1)
    
    plot(isc(:,n1), 'Color', color{conditionList(n1)+1});
    hold on;
    scatter(1:N, isc(:,n1), [], myColors)
    
%     plot3([1 : N], n1*ones(N,1), isc(:,n1), 'Color', color{conditionList(n1)+1}, 'LineWidth', 2);
%     plot3(n1*ones(N,1), [1 : N], isc(:,n1), 'Color', color{conditionList(n1)+1}, 'LineWidth', 2);
    
%     fill3([1 N N 1], [n1 n1 n1 n1], [-.01 -.01 .04 .04], color{conditionList(n1)+1}, 'FaceAlpha', .3, 'LineStyle', 'none');
%     fill3([n1 n1 n1 n1], [1 N N 1], [-.01 -.01 .04 .04], color{conditionList(n1)+1}, 'FaceAlpha', .3, 'LineStyle', 'none');
end


% Z = linkage(isc);
% T = cluster(Z,'maxclust',3);
% cutoff = median([Z(end-2,3) Z(end-1,3)]);
% dendrogram(Z,'ColorThreshold',cutoff);


% gmfit = fitgmdist(score(:,1:10), k);



[idx, C] = kmeans(isc, 2);
    

[isc_tg, perf] = isc2group(X, idx);

% bandWidth = .065;
% [clustCent,data2cluster,cluster2dataCell] = MeanShiftCluster(ISC,bandWidth);

%%
% [standard_data, mu, sigma] = zscore(isc);     % standardize data so that the mean is 0 and the variance is 1 for each variable
% [coeff, score, ~, ~, explained]  = pca(standard_data);     % perform PCA
% new_C = (C-mu)./sigma*coeff;     % apply the PCA transformation to the centroid data
% scatter(score(:, 1), score(:, 2), [], idx)     % plot 2 principal components of the cluster data (three clusters are shown in different colors)
% hold on
% plot(new_C(:, 1), new_C(:, 2), 'kx')     % plot 2 principal components of the centroid data