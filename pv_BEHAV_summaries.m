
%% Load in saved data 

% See jep_fix_cat_xma2_loading_raw_behav_analysis.m for loading in the raw
% data.

clearvars
close all


monkey_names = {'Marta_fix_cat_xma2', 'Max_fix_cat_xma2'};

% load data
EXT_HD = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/';
CCL = '/Users/jonahpearl/Documents/MATLAB/catdog-category-learning';
pv_path = 'XMA2/Monkey_structs';
figureSavePath = '/Users/jonahpearl/Documents/BJR group/Catdog paper/';

Monkeys = struct();
for m = 1:length(monkey_names)
    fname = strcat(monkey_names{m},'_raw_behav','.mat');
    tmp = load(fullfile(EXT_HD, pv_path, fname));
    Monkeys(m).Name =  tmp.Monkeys.Name;
    Monkeys(m).Sessions = tmp.Monkeys.Sessions;
end

%% Fixation errors wrt session
% conclusion: Marta makes 20 - 35 % fixation errors per session, no trend
% upwards or downwards across sessions

for m = 1:length(Monkeys)
    figure
    hold on
    xlabels = {};
    fe_mat = [];
    for i = 1:length(Monkeys(m).Sessions)
        exclude_interlude_fiter = [Monkeys(m).Sessions(i).TrialInfo.After_FE_interlude] == 0;
        total_trials = length(Monkeys(m).Sessions(i).TrialInfo(exclude_interlude_fiter));
%         sum([Monkeys(m).Sessions(i).TrialInfo(exclude_interlude_fiter).Fixn_err] == 1)
        perc_fe = sum([Monkeys(m).Sessions(i).TrialInfo(exclude_interlude_fiter).Fixn_err] == 1) / total_trials;
        bar(i,perc_fe, 'FaceColor', 'b');
        fe_mat(i) = perc_fe;
        xlabels{i} = sprintf('Day %d, %s', i, Monkeys(m).Sessions(i).DateStr);
    end
    xlabel('Session')
    ylabel('Fixation error rate')
    xticks(1:length(Monkeys(m).Sessions))
    xticklabels(xlabels)
    xtickangle(90)
    title(sprintf('%s, fixation errors by session', monkey_names{m}), 'Interpreter', 'none')
    
    % KS test vs uniform distribution
    fe_mat = fe_mat / sum(fe_mat);
    [h,p] = kstest(fe_mat, 'CDF', makedist('Uniform', 1, length(fe_mat)));
end

%% Fixation error timecourse
% conclusion: Marta's fixation error timecourse is very similar across
% sessions, and is more or less linear. Two sessions show a uptick at the
% end of the session, one is mild, 180214 is fairly severe, could consider
% looking at the end of that session more closely.
for m = 1:length(Monkeys)
    figure
    hold on
    for i = 1:length(Monkeys(m).Sessions)
        plot(Monkeys(m).Sessions(i).Fixn_err_TC, 1:length(Monkeys(m).Sessions(i).Fixn_err_TC));
        text(max(Monkeys(m).Sessions(i).Fixn_err_TC)+1,...
            length(Monkeys(m).Sessions(i).Fixn_err_TC)+1,...
            Monkeys(m).Sessions(i).DateStr,...
            'FontSize', 14);
    end
    xlabel('Time in sec')
    ylabel('Fixation errors')
%     title(sprintf('%s, fixation errors over time', monkey_names{m}), 'Interpreter', 'none')
    xticks(0:1000:5000)
    set(gca, 'FontSize', 24)
end

%% Completed trials timecourse
for m = 1:length(Monkeys)
    figure
    hold on
    for i = 1:length(Monkeys(m).Sessions)
        rwd_times = [Monkeys(m).Sessions(i).TrialInfo.Reward_on_time]/1000;
        plot(rwd_times, 1:length(rwd_times));
        text(max(rwd_times)+1, length(rwd_times)+1+i*3, Monkeys(m).Sessions(i).DateStr, 'FontSize', 14);
    end
    xlabel('Time in sec')
    ylabel('Completed trials')
%     title(sprintf('%s, completed trials over time', monkey_names{m}), 'Interpreter', 'none')
    xticks(0:1000:5000)
    set(gca, 'FontSize', 24)
end

%% Fixation errors wrt cue
for m = 1:length(Monkeys)
    concat_sessions = [];
    for i = 1:length(Monkeys(m).Sessions)
        fe_filter = [Monkeys(m).Sessions(i).TrialInfo.Fixn_err]==1;
        cues_of_fe = [Monkeys(m).Sessions(i).TrialInfo(fe_filter).Cue_before_FE]; % turns NaNs into zeros??
        cues_of_fe = cues_of_fe(cues_of_fe ~= 0); % remove the zeros
        figure
        histogram(cues_of_fe, 0:521) % bins include left edge, do not include right edge
        title(sprintf('%s, fixation errors by cue ID, session %d', monkey_names{m}, i))
        xlabel('Cue ID')
        ylabel('Num occurences')
        xticks(1.5:10:521.5)
        xticklabels(1:10:521)
       concat_sessions = [concat_sessions cues_of_fe];
    end
    figure
    histogram(concat_sessions, 0:521)
    ylabel('Num occurences')
    xlabel('Cue ID')
    xticks(1.5:10:521.5)
    xticklabels(1:10:521)
    title(sprintf('%s, fixation errors by cue ID, over all sessions', monkey_names{m}))
end

