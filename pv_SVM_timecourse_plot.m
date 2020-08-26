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

% Load record
fullRecordPath = fullfile(EXT_HD, svmRecordPath);
load(fullRecordPath, 'Record');

% Load behavioral data
load(fullfile(EXT_HD, pv_path, 'MaxMarta_xma2_behav_and_metaNI.mat')) % behavior and neural summaries, but w/o spike times

%% Find desired SVM data
% ID = 570049; % matching at 75%, te and the three individual arrays, with shuffle.
% ID = 844422; % copy of 570049 with unitinfo, cueinfo, etc.
% ID = 439280; % no matching, te, with shuffle.
ID = 317849; % no matching, anterior, middle, posterior, with shuffle.
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

%% Plot the data

% Prepare list of intervals
first_time = -100;
last_time = 450;
width = 100;
step = 5;
starts = first_time:step:last_time;
rIntervals = cell(1, length(starts)); % "relevant Intervals"
for i = 1:length(starts)
    rIntervals{i} = [starts(i) (starts(i) + width)];
end

% Choose what to plot
rSessionsByMonk = {[7 9] [6 7]};
% rArrayLocs = {'SHUFFLE_te', 'te'}; % just treat shuffle as a separate loc, will be easier.
% rArrayLocs = {'te', 'anterior', 'middle', 'posterior', 'SHUFFLE_te', 'SHUFFLE_anterior', 'SHUFFLE_middle', 'SHUFFLE_posterior'}; % just treat shuffle as a separate loc, will be easier.
rArrayLocs = {'anterior', 'middle', 'posterior', 'SHUFFLE_anterior', 'SHUFFLE_middle', 'SHUFFLE_posterior'}; % just treat shuffle as a separate loc, will be easier.
for m = 1:length(Monkeys)
    rSessions = rSessionsByMonk{m};
    
    % Pre-allocate vectors for plotting
    accMeans = zeros(length(rSessions), length(rIntervals), length(rArrayLocs));
    accStds = zeros(length(rSessions), length(rIntervals), length(rArrayLocs));

    
    % Prepare figure
    figure2
    hold on
    
    for iLoc = 1:length(rArrayLocs)
        loc = rArrayLocs{iLoc};
        
        % New subplot for each array location
        subplot(2, length(rArrayLocs)/2, iLoc)
        hold on
        
        for i = 1:length(rSessions)
                sessn = rSessions(i);
                
            for iInt = 1:length(rIntervals)
                interval = rIntervals{iInt};

                % Get field name
                kflID = get_good_interval_name2(interval, loc, 'KFL');
            
                % Get data
                kfls = Data(m).Sessions(sessn).(kflID);
                accMeans(i, iInt, iLoc) = mean(1 - kfls);
                accStds(i, iInt, iLoc) = std(1 - kfls) / numel(kfls);
            end
            
            % Plot a line for each session
            errorbar(starts, accMeans(i, :, iLoc), accStds(i, :, iLoc), ...
                'DisplayName', Monkeys(m).Sessions(sessn).ShortName)
            xlabel('Time from cue on')
            ylabel('SVM abcat accuracy (mean +/- SEM)')
            if m == 1
                ylim([0.45 0.85])
            elseif m == 2
                ylim([0.45 0.65])
            end
            title(sprintf('%s', loc), 'Interpreter', 'none')
            legend()
        end
        
        % Format the plot more.
        sgtitle(sprintf('%s', Monkeys(m).Name), 'Interpreter', 'none')
        
    end
end