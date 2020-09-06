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

ID = 633648; % 75-175, 175-275, 275-375

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

allYLims = [0.3 0.6];
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
        
        % Pre-allocate vectors (proportion = nSig/nTotal but collecting
        % this way makes stats tests easier)
        nSig = zeros(1, length(rSessions));
        nTotal = zeros(1, length(rSessions));
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
            idx = find(cellfun(@(a) all(a == interval), rIntervals));
            
            % Get data
            pVals = Data(m).Sessions(sessn).GLM_Pvals;
            nSig(i) = sum(pVals(units, idx) < glm_alpha);
            nTotal(i) = sum(units);
        end
        
        % Plot the data
        plot(1:length(rSessions), nSig ./ nTotal, 'o-', 'MarkerFaceColor', mlc(1))
        
        % Stats testing if only two sessions (pre and post)
        if numel(rSessions)==2 && prop_test([nSig(1) nSig(2)], [nTotal(1) nTotal(2)], false, glm_alpha)
            yval = max(nSig ./ nTotal);
            plot(1:length(rSessions), repelem(yval*1.02, 1, 2), 'k-')
            scatter(1.5, yval*1.05, 'k*')
        end
        
        
        % Add labels
        ylabel({'Fraction units', 'with signf. diff.'})
%         xticklabels({Data(m).Sessions(rSessions).ShortName})
%         xtickangle(45)
        xticklabels({'Pre', 'Post'})
        xticks(1:length(rSessions))
        
        xlim([0.5 0.5+length(rSessions)])
        % Make graphs have the same y-axes within each monkey
%         ylim(mkYLims{m})
        ylim(allYLims)

        % Make the plot look nice
        formatGLMPlot(gca, gcf)
    end
end
% sgtitle(sprintf('%d to %d', interval(1), interval(2)))
% Save the plots
% pause(0.5)
saveas(gcf, fullfile(figureSavePath, sprintf('GLM_%g_%s_%d_to_%d', ID, loc, interval(1), interval(2))), 'epsc')

%% Functions
function formatGLMPlot(ax, fig)
set(ax, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off',...
    'fontsize',28, 'YGrid', 'on',...
    'fontname', 'Helvetica',...
    'XColor', 'black', 'YColor', 'black')
set(fig, 'Color', 'white')
end

