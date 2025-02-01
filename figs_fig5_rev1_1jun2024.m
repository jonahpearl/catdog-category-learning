% script for Fig 5

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

%% Fig 5A

ID = 21110;  % pre / post, 175-275, matched trial nums. no shuffle.
load(sprintf(fullGLMPath, ID), 'data');
Data = data;
clear data

rSessionsByMonk = {[7 9] [6 7]}; % Fig 2!
rArrayLocs = {'te'};

%% Fig 5A: "bars and stars" for a particular timepoint, for all units

allYLims = [0.1 0.55]; % te
% allYLims = [0.1 0.8]; % arrays
singleInterval = [175 275];
interval = singleInterval;
glm_alpha = 0.05;
            
% Prepare figure
gray = [0.5 0.5 0.5];
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
        colors = zeros(length(rSessions), 3);
        
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
            
            % Set plotting color
            if contains(Monkeys(m).Sessions(sessn).ShortName, "Pre")
                colors(i,:) = mlc(1);
            elseif contains(Monkeys(m).Sessions(sessn).ShortName, "Post")
                colors(i,:) = mlc(2);
            else
                colors(i,:) = gray;
            end
        end
        
        % Plot the data
        % TE
        color = 'k';
        b = bar(1:length(rSessions), nSig ./ nTotal);
        b.FaceColor = 'flat';
        b.CData = colors;
        
        % Stats testing if only two sessions (pre and post)
        if numel(rSessions)==2 && prop_test([nSig(1) nSig(2)], [nTotal(1) nTotal(2)], false, glm_alpha)
            
            % te
            yval = max(nSig ./ nTotal);
            plot(1:length(rSessions), repelem(yval*1.02, 1, 2), 'k-') %
            scatter(1.5, yval*1.05, 'k*') % te
            
            [h, pval, chistat, df] = prop_test([nSig(1) nSig(2)], [nTotal(1) nTotal(2)], false, glm_alpha);
            fprintf('%s, session %s vs %s, signf. incr. in glm signf units (p = %d, chisq %0.5f, df %d \n',...
                Monkeys(m).Name,...
                Monkeys(m).Sessions(rSessions(1)).ShortName, Monkeys(m).Sessions(rSessions(2)).ShortName, ...
                pval, chistat, df)
            disp(pval)
        end
        
        
        % Add labels
        ylabel({'Fraction units', 'with signf. diff.'})
        xlim([0.5 0.5+length(rSessions)])
        xticks([])
%         ylim(allYLims)
        
        % Make the plot look nice
        formatGLMPlot(gca, gcf)
        
    end
end
sgtitle(sprintf('%d to %d', interval(1), interval(2)))

%% Fig 5B: same across arrays

allYLims = [0.1 0.8];
singleInterval = [175 275];
interval = singleInterval;
glm_alpha = 0.05;
            
rArrayLocs = {'anterior', 'middle', 'posterior'};

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
        plot(1:length(rSessions), nSig ./ nTotal, 'o-', ...
            'MarkerFaceColor', mlc(iLoc), 'Color', mlc(iLoc), ...
            'MarkerSize', 10)
  
        % Stats testing if only two sessions (pre and post)
        if numel(rSessions)==2
            
            [h, pval, chistat, df] = prop_test([nSig(1) nSig(2)], [nTotal(1) nTotal(2)], false, glm_alpha);
            fprintf('%s, session %s vs %s, array %s, signf. incr. in glm signf units (p = %.2g, chisq %0.5f, df %d \n',...
                Monkeys(m).Name,...
                Monkeys(m).Sessions(rSessions(1)).ShortName, ...
                Monkeys(m).Sessions(rSessions(2)).ShortName, ...
                loc, ...
                pval, chistat, df)
            
            if pval < glm_alpha
                yval = mean(nSig ./ nTotal);
                scatter(1.5, yval*1.1, 50, mlc(iLoc), '*') % arrays
            end
        end
        
        
        % Add labels
        ylabel({'Fraction units', 'with signf. diff.'})
        xlim([0.5 0.5+length(rSessions)])
        xticks([])
        ylim(allYLims)
        
        % Make the plot look nice
        formatGLMPlot(gca, gcf)
    end
end
sgtitle(sprintf('%d to %d', interval(1), interval(2)))

%% Fig 5D: Histogram of GLM coeffs for signf units 

interval = singleInterval;
glm_alpha = 0.05;
ranksum_alpha = 0.05;       
yval = 'probability';  % normalize hists

rArrayLocs = {"te"};

% Prepare figure
figure2('Position', [1100 1000 250 630])
    
for m = 1:length(Monkeys)
    rSessions = rSessionsByMonk{m};
        
    for iLoc = 1:length(rArrayLocs)
        loc = rArrayLocs{iLoc};
        
        nexttile
        hold on
        
        % Pre-allocate vectors (proportion = nSig/nTotal but collecting
        % this way makes stats tests easier)
        signf_coeffs = cell(length(rSessions), 1);
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
            signf_bools = (pVals(:, idx) < glm_alpha) & units';
            signf_coeffs{i} = Data(m).Sessions(sessn).GLM_coeffs(signf_bools);
            
            % Plot the data
            if strcmp(yval, 'probability')
%                 histogram(abs(signf_coeffs{i}), [0:0.05:1.2], 'Normalization', 'probability')
                histogram(signf_coeffs{i}, [-1.2:0.05:1.2], 'Normalization', 'probability', 'FaceColor', mlc(i))
            elseif strcmp(yval, 'count')
                histogram(signf_coeffs{i}, [-1.2:0.05:1.2], 'FaceColor', mlc(i)) % confirm absolute num is higher post -- yes.
%                 histogram(abs(signf_coeffs{i}), [0:0.05:1.2]) % confirm absolute num is higher post -- yes.
            end
            
            % Draw mean / std on plots
            mu = mean(signf_coeffs{i});
            sigma = std(signf_coeffs{i});
            yl = ylim;
            scatter(mu, yl(2) * 1.1, 100, mlc(i), 'filled', 'v',  'MarkerEdgeColor', 'k')
            plot([mu-sigma mu+sigma], repelem(yl(2) * 1.15, 2), '-', 'Color', mlc(i), 'LineWidth', 1.5)

        end
        
        % Stats testing if only two sessions (pre and post)
        [p, h] = ranksum(abs(signf_coeffs{1}), abs(signf_coeffs{2}));
        if numel(rSessions)==2 && (p < ranksum_alpha)
            fprintf('%s, median of abs GLM coeffs signf higher session %d > %d (p=%3g)\n',...
                Monkeys(m).Name, rSessions(1), rSessions(2), p)
        else
            fprintf('%s, median of abs GLM coeffs no signf diff session %d and %d (p=%3g)\n',...
                Monkeys(m).Name, rSessions(1), rSessions(2), p)
        end
        
        disp(median(signf_coeffs{1}))
        disp(median(signf_coeffs{2}))
        [p, h] = ranksum(signf_coeffs{1}, signf_coeffs{2});
        if numel(rSessions)==2 && (p < ranksum_alpha)
            fprintf('%s, median of SIGNED GLM coeffs different session %d > %d (p=%3g)\n',...
                Monkeys(m).Name, rSessions(1), rSessions(2), p)
        else
            fprintf('%s, median of SIGNED GLM coeffs no signf diff session %d and %d (p=%3g)\n',...
                Monkeys(m).Name, rSessions(1), rSessions(2), p)
        end
        
        
        % Add labels
%         ylabel(yval)
%         xlabel('Signf GLM coeffs.')
        if m == 1
            xlim([-1.2, 1.2])
%             xlim([0, 1.1])
        else
            xlim([-0.6, 0.6])
%             xlim([0, 0.5])
        end
        
        % Make the plot look nice
        formatGLMPlot(gca, gcf)
    end
end

%% Fig 5C: load data

ID = 49218;  % ditto but all sessions
load(sprintf(fullGLMPath, ID), 'data');
Data = data;
clear data

rSessionsByMonk = {[1 2 3 5 6 7 9], 1:7};
rArrayLocs = {'te'};

%% Fig 5C: "bars and stars" across baseline sessions


allYLims = [0.1 0.55]; % te
% allYLims = [0.1 0.8]; % arrays
singleInterval = [175 275];
interval = singleInterval;
glm_alpha = 0.05;
            
% Prepare figure
gray = [0.5 0.5 0.5];
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
        colors = zeros(length(rSessions), 3);
        
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
            
            % Set plotting color
            if contains(Monkeys(m).Sessions(sessn).ShortName, "Pre")
                colors(i,:) = mlc(1);
            elseif contains(Monkeys(m).Sessions(sessn).ShortName, "Post")
                colors(i,:) = mlc(2);
            else
                colors(i,:) = gray;
            end
        end
        
        % Plot the data
        % TE
        color = 'k';
        b = bar(1:length(rSessions), nSig ./ nTotal);
        b.FaceColor = 'flat';
        b.CData = colors;
        
        % Stats testing if only two sessions (pre and post)
        if numel(rSessions)==2 && prop_test([nSig(1) nSig(2)], [nTotal(1) nTotal(2)], false, glm_alpha)
            
            % te
            yval = max(nSig ./ nTotal);
            plot(1:length(rSessions), repelem(yval*1.02, 1, 2), 'k-') %
            scatter(1.5, yval*1.05, 'k*') % te
            
            [h, pval, chistat, df] = prop_test([nSig(1) nSig(2)], [nTotal(1) nTotal(2)], false, glm_alpha);
            fprintf('%s, session %s vs %s, signf. incr. in glm signf units (p = %d, chisq %0.5f, df %d \n',...
                Monkeys(m).Name,...
                Monkeys(m).Sessions(rSessions(1)).ShortName, Monkeys(m).Sessions(rSessions(2)).ShortName, ...
                pval, chistat, df)
            disp(pval)
        end
        
        
        % Add labels
        ylabel({'Fraction units', 'with signf. diff.'})
        xlim([0.5 0.5+length(rSessions)])
        xticks([])
%         ylim(allYLims)
        
        % Make the plot look nice
        formatGLMPlot(gca, gcf)
        
    end
end
sgtitle(sprintf('%d to %d', interval(1), interval(2)))

%% Functions
function formatGLMPlot(ax, fig)
set(ax, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off',...
    'fontsize',28, 'YGrid', 'on',...
    'fontname', 'Helvetica',...
    'XColor', 'black', 'YColor', 'black')
set(fig, 'Color', 'white')
end

