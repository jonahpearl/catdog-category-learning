% plot SVM timecourse results

%% Params and load SVM record
clearvars
close all

% Define paths to data
EXT_HD = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/';
CCL = '/Users/jonahpearl/Documents/MATLAB/catdog-category-learning';
svmRecordPath = 'XMA2/Monkey_structs/SVM_Records.mat';
pv_path = 'XMA2/Monkey_structs';
fullSVMPath = fullfile(EXT_HD, pv_path, 'SVM_results_%g.mat');
figureSavePath = '/Users/jonahpearl/Documents/BJR group/Catdog paper/';

% Load record
fullRecordPath = fullfile(EXT_HD, svmRecordPath);
load(fullRecordPath, 'Record');

% Load behavioral data
load(fullfile(EXT_HD, pv_path, 'MaxMarta_xma2_behav_and_metaNI.mat')) % behavior and neural summaries, but w/o spike times

% Statistics parameters
cluster_alpha = 0.05;

%% Find desired SVM data
% ID = 570049; % te and the three individual arrays, matching at 75%, with shuffle.
% ID = 844422; % copy of 570049 with unitinfo, cueinfo, etc.
% ID = 439280; % te, no matching, with 5x shuffle.
ID = 958736; % te, no matching, with 100x shuffle
% ID = 317849; % anterior, middle, posterior, no matching, with shuffle.

fNames = fields(Record);
row = find([Record.ID] == ID);
for f = 1:length(fNames)
    s = Record(row).(fNames{f});
    field = fNames{f};
    if strcmp(field, 'MatchInputParams')
        fprintf('%s: %s \n', field, mat2str(s))
    elseif strcmp(field, 'SessionsUsed')
        fprintf('%s: %s--%s, %s--%s \n', field, Monkeys(1).Name, mat2str(s{1}), Monkeys(2).Name, mat2str(s{2}))
    elseif isa(s, 'cell')
        fprintf('%s: %s \n', field, join(cellfun(@string ,s), ', '))
    else
        fprintf('%s: %s \n', field, string(s))
    end
    
end

%% Load SVM data
load(sprintf(fullSVMPath, ID), 'data');
Data = data;
clear data

%% Prepare list of intervals
first_time = -100;
last_time = 450;
width = 100;
step = 5;
starts = first_time:step:last_time;
rIntervals = cell(1, length(starts)); % "relevant Intervals"
for i = 1:length(starts)
    rIntervals{i} = [starts(i) (starts(i) + width)];
end

%% Choose what to plot

rSessionsByMonk = {[7 9] [6 7]};
% rArrayLocs = {'te', 'SHUFFLE_te'}; % just treat shuffle as a separate loc, will be easier.
rArrayLocs = {'te'};
% rArrayLocs = {'te', 'anterior', 'middle', 'posterior', 'SHUFFLE_te', 'SHUFFLE_anterior', 'SHUFFLE_middle', 'SHUFFLE_posterior'}; % just treat shuffle as a separate loc, will be easier.
% rArrayLocs = {'anterior', 'middle', 'posterior', 'SHUFFLE_anterior', 'SHUFFLE_middle', 'SHUFFLE_posterior'}; % just treat shuffle as a separate loc, will be easier.

%% Run cluster-based permutation statistics (shuffle vs real, pre vs post)

% This section of code automatically looks at shuffled data -- do not 
% include in list of locs.
rArrayLocs_stats = {'te'};

% 1. run ttests at each time point
% 2. Generate clusters and cluster statistics (real and shuffled)
% 3. Compare ranked statistics from shuffle to real to find pvals (mult
% comp correction)

for m = 1:length(Monkeys)
    rSessions = rSessionsByMonk{m};
    if numel(rSessions) ~= 2
        error('This code expects two sessions to compare')
    end
    
    for iLoc = 1:length(rArrayLocs_stats)
        loc = rArrayLocs_stats{iLoc};

        % Get list of real tstats
        tStats = getrealtstats(Data, m, rSessions, loc, rIntervals); % 1 x num intervals
        
        % Get lists of tstats from shuffled data
        tStats_SHUFFLED = getshuffledtstats(Data, m, rSessions, loc, rIntervals); % num shuffle permutations x num intervals
        
        % Find real clusters
        % clusterStats: 1 x num clusters vector (num clusters calculated in
        % function). Is sorted.
        % intsInCluster: of the I intervals in rIntervals, which ones are
        % in each cluster? Sorted to correspond to clusterStats.
        
        [clusterStats, intsInCluster] = calculateTClusterStats(tStats); 
        
        % Find shuffled clusters
        nShuff = size(tStats_SHUFFLED,1);
        nClust = size(clusterStats,2);
        clusterStats_SHUFFLED = zeros(nClust, nShuff);
        for iShuff = 1:nShuff
            clusterStats_SHUFFLED(:,iShuff) = calcShuffClusterStats(tStats_SHUFFLED(iShuff,:), nClust);
        end
        
        % Take largest from each shuffled cluster to use to create p-values.
        % (Can also rank-match, so it was worth getting the entire shuffled
        % cluster vector in case we want to change the pvalue method.)
        shuffled_vals = clusterStats_SHUFFLED(1,:);
        
        % Generate p-values. The transpose with the <= operator makes
        % MATLAB compare all to all. Love a good one-liner.
        signf_clusters = find(sum(clusterStats <= shuffled_vals') / numel(shuffled_vals) < cluster_alpha);
        
        % Of the intervals in rIntervals, which are significant?
        signf_ints = [intsInCluster{signf_clusters}];
        sigID = sprintf('ClustPermSignfInds_Sessions_%d_vs_%d', rSessions(1), rSessions(2));
        Monkeys(m).(sigID) = signf_ints;
        
    end
end

%% Plot the data

% Plotting params
plot_alpha = 0.4; % transparency of sem fill
mkYLims = {[0.45 0.85], [0.45 0.65]};


for m = 1:length(Monkeys)
    rSessions = rSessionsByMonk{m};
    
    % Pre-allocate vectors for plotting
    accMeans = zeros(length(rSessions), length(rIntervals), length(rArrayLocs));
    accSems = zeros(length(rSessions), length(rIntervals), length(rArrayLocs));
    
    % Prepare figure
    figure2
    hold on
    
    for iLoc = 1:length(rArrayLocs)
        loc = rArrayLocs{iLoc};
        
        % New subplot for each array location
%         subplot(2, length(rArrayLocs)/2, iLoc)
        hold on
        
        for i = 1:length(rSessions)
            sessn = rSessions(i);
            
            for iInt = 1:length(rIntervals)
                interval = rIntervals{iInt};
                
                % Get field names
                kflID = get_good_interval_name2(interval, loc, 'KFL');
%                 kflID_shuffle = get_good_interval_name2(interval, loc, 'KFL_SHUFFLE');
                
                % Get data
                kfls = Data(m).Sessions(sessn).(kflID);
                kfls = kfls(:); % reshape into 1 x numel
                accMeans(i, iInt, iLoc) = mean(1 - kfls);
                accSems(i, iInt, iLoc) = std(1 - kfls) / numel(kfls);
            end
            
            
            % Simplify variable names
            meanVector = accMeans(i, :, iLoc);
            semVector = accSems(i, iInt, iLoc);
            
            % Plot a line for each session with filled sem.
            % mlc() is a function to get the RGB vals for the default
            % MATLAB colors.
            plot(starts, meanVector, '-', 'LineWidth', 2, 'Color', mlc(i))
            fill([starts fliplr(starts)],...
                [meanVector + semVector, fliplr(meanVector - semVector)], mlc(i),...
                'FaceAlpha', plot_alpha,'linestyle','none');
            
            % Plot signf inds
            sigID = sprintf('ClustPermSignfInds_Sessions_%d_vs_%d', rSessions(1), rSessions(2));
            plot(starts(Monkeys(m).(sigID)), mkYLims{m}(2)-0.02, 'ko', 'MarkerFaceColor', 'k')
            
            % Add labels
            xlabel('Time from cue on (ms)')
            ylabel('SVM accuracy')
            
            % Make graphs have the same y-axes within each monkey
            ylim(mkYLims{m})
            
            
            
            % Detailed labels if desired
            %             ylabel('SVM abcat accuracy (mean +/- SEM)')
            %             title(sprintf('%s', loc), 'Interpreter', 'none')
            %             legend()
            
            % Make the plot look nice
            formatSVMPlot(gca, gcf)
        end
        
        % More detailed labels if desired.
        %         sgtitle(sprintf('%s', Monkeys(m).Name), 'Interpreter', 'none')
        
        % Save the plots
        saveas(gcf, fullfile(figureSavePath, sprintf('SVM_%s.svg', sigID)))
    end
end

%% Functions
function formatSVMPlot(ax, fig)
set(ax, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off',...
    'fontsize',28, 'YGrid', 'on',...
    'fontname', 'Helvetica',...
    'XColor', 'black', 'YColor', 'black')
set(fig, 'Color', 'white')
end

function tStats = getrealtstats(Data, m, rSessions, loc, rIntervals)
    tStats = zeros(1, length(rIntervals));
    for iInt = 1:length(rIntervals)
        interval = rIntervals{iInt};

        % Get field names
        kflID = get_good_interval_name2(interval, loc, 'KFL');

        % Run t-tests
        s1 = Data(m).Sessions(rSessions(1)).(kflID);
        s2 = Data(m).Sessions(rSessions(2)).(kflID);
        [~,~,~,stats] = ttest2(s1, s2, 'Tail', 'both');
        tStats(iInt) = stats.tstat;
    end
end

function tStats_SHUFFLED = getshuffledtstats(Data, m, rSessions, loc, rIntervals)

    % First calculate number of shuffles
    kflID = get_good_interval_name2(rIntervals{1}, loc, 'KFL_SHUFFLE');
    exampleKFLs = Data(m).Sessions(rSessions(1)).(kflID);
    nShuffles = size(exampleKFLs,2);
    
    tStats_SHUFFLED = zeros(nShuffles, length(rIntervals));
    for iInt = 1:length(rIntervals)
        interval = rIntervals{iInt};
        
        % Get field names
        kflID = get_good_interval_name2(interval, loc, 'KFL_SHUFFLE');

        % Get data
        s1 = Data(m).Sessions(rSessions(1)).(kflID);
        s2 = Data(m).Sessions(rSessions(2)).(kflID);
        
        for iShuff = 1:nShuffles
            [~,~,~,stats] = ttest2(s1(:,iShuff), s2(:, iShuff), 'Tail', 'both');
            tStats_SHUFFLED(iShuff, iInt) = stats.tstat;
        end
    end
end

function clusterStats = calcShuffClusterStats(tStats, nClust)
    
    % Run the same algorithm as the real data
    clusterStats = calculateTClusterStats(tStats);
    
    % make same length as real data
    if length(clusterStats) < nClust
        clusterStats(length(clusterStats):nClust) = 0;
    elseif length(clusterStats) > nClust
        clusterStats(nClust+1:end) = []; % sort is descending, so this removes smallest values
    end
end

function [clusterStats, intsInCluster] = calculateTClusterStats(tStats)
    starts = [1 find(diff(sign(tStats)))+1];
    ends = [starts(2:end)-1 length(tStats)];
    clusters = cell(1, length(starts));
    intsInCluster = cell(1, length(starts)); % inds in rIntervals that belong to each cluster, for plotting later
    for iClust = 1:length(clusters)
        clusters{iClust} = tStats(starts(iClust):ends(iClust));
        intsInCluster{iClust} = starts(iClust):ends(iClust);
    end
    [clusterStats, idx] = sort(abs(cellfun(@sum, clusters)), 'descend');
    intsInCluster = intsInCluster(idx); % sort these to correspond to sorted cluster list
end