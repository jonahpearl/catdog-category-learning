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

% ID = 633648; % 75-175, 175-275, 275-375
% ID = 748804; % all days, same 3 ints
% ID = 224797; % same as above but with GLM coeffs

% === revisions ===
ID = 21110;  % pre / post, 175-275, matched trial nums. no shuffle.
% ID = 49218;  % ditto but all sessions

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
rSessionsByMonk = {[7 9] [6 7]}; % Fig 2!
% rSessionsByMonk = {[1 2 3 5 6 7 9], 1:7};
% rSessionsByMonk = {[1 6 7 9], [1 5 6 7]};

% Choose arrays. Treat shuffle as a separate loc, will be easier.
% rArrayLocs = {'te', 'SHUFFLE_te'}; 
rArrayLocs = {'te'};
% rArrayLocs = {'te', 'anterior', 'middle', 'posterior', 'SHUFFLE_te', 'SHUFFLE_anterior', 'SHUFFLE_middle', 'SHUFFLE_posterior'};
% rArrayLocs = {'anterior', 'middle', 'posterior', 'SHUFFLE_anterior', 'SHUFFLE_middle', 'SHUFFLE_posterior'}; 
% rArrayLocs = {'te', 'anterior', 'middle', 'posterior'}; 
% rArrayLocs = {'anterior', 'middle', 'posterior'}; 

%% Plot "bars and stars" for a particular timepoint, for all units

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


%% PLot same across arrays

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

%% Familiarity controls

% mkYLims = {[0.5 0.75], [0.49 0.56]};
% allYlims = [0.48 0.82];
singleInterval = [175 275];
interval = singleInterval;
% singleInterval = [175 350];

% Prepare figure
% figure2('Position', [2200 1300 250 630])
% figure2('Position', [2200 1300 400 1200])
figure
hold on

rSessionsByMonk = {[1 6 7 9], [1 5 6 7]}; % [baseFirst, baseLast, pre, post]
rArrayLocs = {'anterior', 'middle', 'posterior', 'te'};

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
            
            pVals = Data(m).Sessions(sessn).GLM_Pvals;
            nSig(i) = sum(pVals(units, idx) < glm_alpha);
            nTotal(i) = sum(units);
        end
        
        % Baseline days
        plot([1 2], nSig([1 2]) ./ nTotal([1 2]),... 
            'o-',...
            'LineWidth', 2,...
            'Color', mlc(iLoc),...
            'DisplayName', loc)

        % pre and post
        plot([3 4], nSig([3 4]) ./ nTotal([3 4]),...
            'o-',...
            'LineWidth', 2,...
            'Color', mlc(iLoc),...
            'HandleVisibility', 'off')
        
        % Add labels
        ylabel('Prop. signf. (GLM)')
        xticks(1:length(rSessions))
        xticklabels({Data(m).Sessions(rSessions).ShortName})
        xtickangle(45)
        xlim([0.5 0.5+length(rSessions)])
%         legend()
        title('Familiarity comparisons')

        % Make the plot look nice
        formatGLMPlot(gca, gcf)
    end
end

%% Histogram of GLM coeffs for signf units 

singleInterval = [175 275];
interval = singleInterval;
glm_alpha = 0.05;
ranksum_alpha = 0.05;       
rArrayLocs = {'te'};
yval = 'probability';  % normalize hists
% yval = 'count';  % don't normalize hists

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
% sgtitle(sprintf('%d to %d', interval(1), interval(2)))
% Save the plots
% pause(0.5)
% saveas(gcf, fullfile(figureSavePath, sprintf('GLM_%g_%s_%d_to_%d', ID, loc, interval(1), interval(2))), 'epsc')

%% Mean and std of GLM coeffs across days

% Params
singleInterval = [175 275];
interval = singleInterval;
glm_alpha = 0.05;
loc = 'te';

for m = 1:length(Monkeys)
    rSessions = rSessionsByMonk{m};
    
    for i = 1:length(rSessions)
        sessn = rSessions(i);
        
        % Get interval to use
        rIntervals = Data(m).Sessions(sessn).GLM_intervals;
        idx = find(cellfun(@(a) all(a == interval), rIntervals));
            
        % Get data
        units = strcmp({Monkeys(m).Sessions(sessn).UnitInfo.Area}, loc);
        pVals = Data(m).Sessions(sessn).GLM_Pvals;
        signf_bools = (pVals(:, idx) < glm_alpha) & units';
        signf_coeffs = Data(m).Sessions(sessn).GLM_coeffs(signf_bools);
        fprintf('%s, sessn %d, signf. GLM coeffs in %s — mean or std of (exp(abs(coeff))): mean %0.2f, std %0.2f\n', ...
            Monkeys(m).Name, sessn, loc, mean(exp(abs(signf_coeffs))), std(exp(abs(signf_coeffs))))
    end
end

%% Functions
function formatGLMPlot(ax, fig)
set(ax, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off',...
    'fontsize',28, 'YGrid', 'on',...
    'fontname', 'Helvetica',...
    'XColor', 'black', 'YColor', 'black')
set(fig, 'Color', 'white')
end

