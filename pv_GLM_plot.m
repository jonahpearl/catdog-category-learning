% plot SVM timecourse results

%% Params and load GLM record
clearvars
close all

% Define paths to data
EXT_HD = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/';
CCL = '/Users/jonahpearl/Documents/MATLAB/catdog-category-learning';
glmRecordPath = 'XMA2/Monkey_structs/GLM_Records.mat';
pv_path = 'XMA2/Monkey_structs';
fullGLMPath = fullfile(EXT_HD, pv_path, 'GLM_results_%g.mat');
figureSavePath = '/Users/jonahpearl/Documents/BJR group/Catdog paper/';

% Load record
fullRecordPath = fullfile(EXT_HD, glmRecordPath);
load(fullRecordPath, 'Record');

% Load behavioral data
load(fullfile(EXT_HD, pv_path, 'MaxMarta_xma2_behav_and_metaNI.mat')) % behavior and neural summaries, but w/o spike times

% Statistics parameters
cluster_alpha = 0.05;

% add Area labels
for m = 1:length(Monkeys)
    for i = 1:length(Monkeys(m).Sessions)
        for j = 1:length(Monkeys(m).Sessions(i).UnitInfo)
            switch Monkeys(m).Sessions(i).UnitInfo(j).Location
                case 'posterior'
                    Monkeys(m).Sessions(i).UnitInfo(j).Area = 'te';
                case 'middle'
                    Monkeys(m).Sessions(i).UnitInfo(j).Area = 'te';
                case 'anterior'
                    Monkeys(m).Sessions(i).UnitInfo(j).Area = 'te';
                case 'teo'
                    Monkeys(m).Sessions(i).UnitInfo(j).Area = 'teo';
            end
        end
    end
end

%% Find desired GLM data

ID = 20752;

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

%% Load GLM data
load(sprintf(fullGLMPath, ID), 'data');
Data = data;
clear data

%% Choose what to plot

% Choose sessions.
rSessionsByMonk = {[7 9] [6 7]};

% Choose arrays. Treat shuffle as a separate loc, will be easier.
% rArrayLocs = {'te', 'SHUFFLE_te'}; 
rArrayLocs = {'te'};
% rArrayLocs = {'te', 'anterior', 'middle', 'posterior', 'SHUFFLE_te', 'SHUFFLE_anterior', 'SHUFFLE_middle', 'SHUFFLE_posterior'};
% rArrayLocs = {'anterior', 'middle', 'posterior', 'SHUFFLE_anterior', 'SHUFFLE_middle', 'SHUFFLE_posterior'}; 
% rArrayLocs = {'anterior', 'middle', 'posterior'}; 

%% Plot "bars and stars" for a particular timepoint

mkYLims = {[0.45 0.85], [0.45 0.65]};
singleInterval = [175 275];
interval = singleInterval;
glm_alpha = 0.05;
            
% Prepare figure
figure2('Position', [1100 1000 250 630])
hold on
    
for m = 1:length(Monkeys)
    rSessions = rSessionsByMonk{m};
        
    % New subplot for each monkey
    subplot(length(Monkeys), 1, m)
    hold on
        
    for iLoc = 1:length(rArrayLocs)
        loc = rArrayLocs{iLoc};
        
        propSig = zeros(1, length(rSessions));
        for i = 1:length(rSessions)
            sessn = rSessions(i);
            
            % Get indices of units in the pValues matrix to use in calculating
            % proprotion of units with signf GLMs
            if strcmp(loc, 'te')
                units = strcmp({Monkeys(m).Sessions(sessn).UnitInfo.Area}, loc);
            else
                units = strcmp({Monkeys(m).Sessions(sessn).UnitInfo.Location}, loc);
            end
            
            % Get interval to use
            rIntervals = Data(m).Sessions(sessn).GLM_intervals;
            for iInt = 1:length(rIntervals)
                if all(rIntervals{iInt} == singleInterval)
                    int = iInt;
                    break
                else 
                    int = -1;
                end
            end
            if int == -1
                error('Requested interval not found in interval list')
            end
            
            % Get data
            pVals = Data(m).Sessions(sessn).GLM_Pvals;
            propSig(i) = sum(pVals(units, int) < glm_alpha) / sum(units);
        end
        
        % Plot the data
        plot(1:length(rSessions), propSig)
        
        % Stats testing if only two sessions (pre and post)
        if numel(rSessions)==2 && 0
            
        end
        
        
        % Add labels
        ylabel('SVM accuracy')
%         xticklabels({Data(m).Sessions(rSessions).ShortName})
%         xtickangle(45)
        xticklabels({'Pre', 'Post'})
        xticks(1:length(rSessions))
        
        xlim([0.5 0.5+length(rSessions)])
        % Make graphs have the same y-axes within each monkey
%         ylim(mkYLims{m})

        % Make the plot look nice
        formatSVMPlot(gca, gcf)
    end
end
sgtitle(sprintf('%d to %d', interval(1), interval(2)))
% Save the plots
% pause(0.5)
% saveas(gcf, fullfile(figureSavePath, sprintf('SVM_%s_%d_%s_%d_to_%d', sigID, ID, loc, interval(1), interval(2))), 'epsc')

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