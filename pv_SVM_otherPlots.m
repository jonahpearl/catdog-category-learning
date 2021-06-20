% other SVM plots


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

% Note: need 100x shuffle for permutation testing.

% ID = 570049; % te and the three individual arrays, matching at 75%, with shuffle.
% ID = 844422; % copy of 570049 with unitinfo, cueinfo, etc.
% ID = 439280; % te, no matching, with 5x shuffle.
% ID = 958736; % te, no matching, with 100x shuffle. FIG 2!
% ID = 317849; % anterior, middle, posterior, no matching, with 5x shuffle.
% ID = 886768; % all locs, no matching, 5x shuffle, all days, three main intervals.
% ID = 754886; % all locs, matching, no shuffle, all days, all ints.
    % This one isn't good for Max, b/c three of his baseline days have very
    % low trial counts, and it drags down the decoding accuracy for the
    % other four.
ID = 93591; % all locs, no matching, no shuffle, all days, all ints. SUPP FIG 6!
% ID = 467637; % Variable bin sizes! all locs, no matching, no shuffle, all days.
% ID = 696349; % copy of above
% ID = 983557; % Variable bin sizes as above but with matching.
% ID = 428279; % just pre and post, te, matched.

% ID = 178096; % image subset analysis: train and test on 240/240, just te
ID = 977195; % image subset analysis: train and test on 240/240, all locs / days
% ID = 772151; % image subset analysis: train on 240, test on 20, all locs / days
% ID = 339824; % image subset analysis: train on 20, test on 240, all locs / days

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

%% Image subset analyses

%% Familiarity comparisons

% mkYLims = {[0.5 0.75], [0.49 0.56]};
% allYlims = [0.48 0.82];
singleInterval = [175 275];
% singleInterval = [175 350];

% Prepare figure
% figure2('Position', [2200 1300 250 630])
figure2('Position', [2200 1300 400 1200])
% figure
hold on

rSessionsByMonk = {[1 6 7 9], [1 5 6 7]}; % [baseFirst, baseLast, pre, post]
rArrayLocs = {'te', 'anterior', 'middle', 'posterior'};

for m = 1:length(Monkeys)
    rSessions = rSessionsByMonk{m};
    
    % Pre-allocate vectors for plotting
    accMeans = zeros(length(rSessions), 1, length(rArrayLocs));
    accSems = zeros(length(rSessions), 1, length(rArrayLocs));
        
    % New subplot for each monkey
    subplot(length(Monkeys), 1, m)
    hold on
        
    for iLoc = 1:length(rArrayLocs)
        loc = rArrayLocs{iLoc};
        
        for i = 1:length(rSessions)
            sessn = rSessions(i);
            interval = singleInterval;

            % Get field names
            kflID = get_good_interval_name2(interval, loc, 'KFL');

            % Get data
            kfls = Data(m).Sessions(sessn).(kflID);
            kfls = kfls(:); % reshape into 1 x numel
            accMeans(i, 1, iLoc) = mean(1 - kfls);
            accSems(i, 1, iLoc) = std(1 - kfls) / numel(kfls);
        end
        
        % Baseline days
        errorbar([1 2], accMeans([1 2],1,iLoc), accSems([1 2],1,iLoc),...
            '-',...
            'LineWidth', 2,...
            'Color', mlc(iLoc),...
            'DisplayName', loc)

        % pre and post
        errorbar([3 4], accMeans([3 4],1,iLoc), accSems([3 4],1,iLoc),...
            '-',...
            'LineWidth', 2,...
            'Color', mlc(iLoc),...
            'HandleVisibility', 'off')
        
        % Add labels
        ylabel('SVM accuracy')
        xticks(1:length(rSessions))
        xticklabels({Data(m).Sessions(rSessions).ShortName})
        xtickangle(45)
        xlim([0.5 0.5+length(rSessions)])
%         legend
        title('Familiarity comparisons')

        % Make the plot look nice
        formatSVMPlot(gca, gcf)
    end
end
% sgtitle(sprintf('%d to %d', interval(1), interval(2)))
% Save the plots
% pause(0.5)
% saveas(gcf, fullfile(figureSavePath, sprintf('SVM_%s_%d_%s_%d_to_%d', sigID, ID, loc, interval(1), interval(2))), 'epsc')

%% Functions

function formatSVMPlot(ax, fig)
set(ax, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off',...
    'fontsize',22, 'YGrid', 'on',...
    'fontname', 'Helvetica',...
    'XColor', 'black', 'YColor', 'black')
set(fig, 'Color', 'white')
end