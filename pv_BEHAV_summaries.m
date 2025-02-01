
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

%% If needed, get date strs and short names of sessions
marta_xls = fullfile(EXT_HD, 'RecordingMarta.xlsx');
max_xls = fullfile(EXT_HD, 'RecordingMax.xlsx');

for m = 1:length(Monkeys)
    date_strs = {Monkeys(m).Sessions.DateStr};
    if regexp(Monkeys(m).Name, 'Marta\w*')
%         shortNames = get_short_names(date_strs, marta_xls, '\w*cats\w*');
        shortNames = get_short_names(date_strs, marta_xls, '\w*cars\w*');
    else
%         shortNames = get_short_names(date_strs, max_xls, '\w*cats\w*');
        shortNames = get_short_names(date_strs, max_xls, '\w*cars\w*');
    end
    for i = 1:length(Monkeys(m).Sessions)
        Monkeys(m).Sessions(i).ShortName = shortNames{i};
    end
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
rSessionsByMonk = {[7 9], [6 7]};  % (Fig 2B / 2D)

for m = 1:length(Monkeys)
    figure
    hold on
    rSessions = rSessionsByMonk{m};
    for i = 1:length(rSessions)
        sessn = rSessions (i);
        rwd_times = [Monkeys(m).Sessions(sessn).TrialInfo.Reward_on_time]/1000;
        plot(rwd_times, 1:length(rwd_times), 'LineWidth', 2);
        disp(length(rwd_times) / rwd_times(end))
%         text(max(rwd_times)+1, length(rwd_times)+1+i*3, Monkeys(m).Sessions(sessn).DateStr, 'FontSize', 14);
    end
    xlabel('Time in sec')
    ylabel('Completed trials')
%     title(sprintf('%s, completed trials over time', monkey_names{m}), 'Interpreter', 'none')
%     xticks(0:1000:5000)
    set(gca, 'FontSize', 24)
    legend(["Pre", "Post"])
end

%% Fixation errors wrt cue
rSessionsByMonk = {[7 9], [6 7]};  % (Fig 2B / 2D)

for m = 1:length(Monkeys)
    concat_sessions = [];
    rSessions = rSessionsByMonk{m};
    
    for i = 1:length(rSessions)
        sessn = rSessions(i);
        fe_filter = [Monkeys(m).Sessions(sessn).TrialInfo.Fixn_err]==1;
        cues_of_fe = [Monkeys(m).Sessions(sessn).TrialInfo(fe_filter).Cue_before_FE]; % turns NaNs into zeros??
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

%% Fixation errors wrt category

rSessionsByMonk = {[7 9], [6 7]};  % (Fig 2B / 2D)
catg2_ind1 = 261;

for m = 1:length(Monkeys)
    concat_sessions = [];
    rSessions = rSessionsByMonk{m};
    
    fe_data = zeros(2, length(rSessions));  % rows = cat/dog, cols = sessions.
    n_trs = zeros(2, length(rSessions));
    for i = 1:length(rSessions)
        sessn = rSessions(i);
        fe_filter = [Monkeys(m).Sessions(sessn).TrialInfo.Fixn_err]==1;
        cues_of_fe = [Monkeys(m).Sessions(sessn).TrialInfo(fe_filter).Cue_before_FE]; % turns NaNs into zeros??
        cues_of_fe = cues_of_fe(cues_of_fe ~= 0); % remove the zeros
        cues_seen = [Monkeys(m).Sessions(sessn).TrialInfo.Cues_seen];
        cues_seen = cues_seen(cues_seen ~= 0);
        
        cat_fes = sum(cues_of_fe < catg2_ind1);
        dog_fes = sum(cues_of_fe >= catg2_ind1);
        fe_data(:,i) = [cat_fes; dog_fes];
        n_trs(:,i) = [sum(cues_seen < catg2_ind1); ...
            sum(cues_seen >= catg2_ind1)];
    end
    figure
    bar((fe_data./n_trs)', 'grouped')
    ylabel('Fraction of category stimuli')
    xlabel('Session')
    legend(["Cat", "Dog"])
    title(sprintf('%s, fixation errors by catg, over all sessions', monkey_names{m}))
    
    [h,p_cats,~,~] = prop_test(fe_data(1,:), n_trs(1,:), false, 0.05);
    [h,p_dogs,~,~] = prop_test(fe_data(2,:), n_trs(2,:), false, 0.05);
    
    fprintf("P-val for prop. cat FE differing pre vs. post: %g \n", p_cats)
    fprintf("P-val for prop. dog FE differing pre vs. post: %g \n", p_dogs)
    
    formatSVMPlot(gca, gcf)
    saveas(gcf, sprintf("/Users/jonahpearl/Documents/BJR group/Catdog_paper/20240423_revisions/behavior/%s_fixn_errs.eps", Monkeys(m).Name))
    
end


%% Functions
function formatSVMPlot(ax, fig, fontsize)
if nargin == 2
    fontsize=28;
end
set(ax, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off',...
    'fontsize', fontsize, 'YGrid', 'on',...
    'fontname', 'Helvetica',...
    'XColor', 'black', 'YColor', 'black')
set(fig, 'Color', 'white')
end
