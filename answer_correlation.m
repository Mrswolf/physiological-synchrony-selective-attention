
function [r_ans, p_ans, rank] = answer_correlation( isc, answers)

[N, ~, ~, ~] = size(isc);

% initialize variables for faster processing
varname = fieldnames(answers); rank = zeros(N, length(varname));
r_ans = zeros(2, length(varname), size(isc, 3)); p_ans = zeros(2,length(varname), size(isc, 3));

% for beeps, performance is based on the absolute distance from the correct
% answers (a much lower estimation and much higher estimation are both
% incorrect). To align beeps with the other measures we take the inverse.
% Then for all measures holds: the higher, the better.
answers.beeps = abs(answers.beeps);
answers.beeps = 1./(answers.beeps + 1);

for ss = 1 : length(varname)

    for mm = 1 : size(isc, 3)

        % provide score for each participants based on rank. Best answer gets
        % rank 26, worst answer gets rank 1
        rank(:,ss) = tiedrank(answers.(varname{ss}));

        % compute correlations between isc with respect to group g and the rank
        % performance of answers j
        for g = 1 : 2

            [r_ans(g,ss,mm), p_ans(g,ss,mm)] = corr(isc(:,g,mm), rank(:,ss), ...
                'rows', 'complete');

        end
        
    end
    
end

% correlate difference isc with difference in question answers
score = zeros(N, 2);
score(:,1) = rank(:, strcmp(varname, 'narrative'));
score(:,2) = mean(rank(:, strcmp(varname, 'iads') | strcmp(varname, 'beeps')),   2);
rank(:,4) = diff(score, [], 2);

for mm = 1 : size(isc, 3)

    [r_ans(1,4,mm), p_ans(1,4,mm)] = corr(diff(isc(:,:,mm), [], 2), rank(:,4), 'rows', 'complete');
    
end

end