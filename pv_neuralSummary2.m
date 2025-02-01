% look at neural firing rates across images

%% Load data
clearvars
close all

% Define paths to data
EXT_HD = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/';
CCL = '/Users/jonahpearl/Documents/MATLAB/catdog-category-learning';
pv_path = 'XMA2/Monkey_structs';
spikeCountPath = 'XMA2/Spike_count_mats';

% Load behavioral data
load(fullfile(EXT_HD, pv_path, 'MaxMarta_xma2_behav_and_metaNI.mat')) % behavior and neural summaries, but w/o spike times
TE_LOCS = {'anterior', 'middle', 'posterior'};

% Interval parameters for loading spike count mat
step = 5;
width = 100;

ignoreVal = 20;

%% Add 'area' label to unitinfo, get short session names

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

marta_xls = fullfile(EXT_HD, 'RecordingMarta.xlsx');
max_xls = fullfile(EXT_HD, 'RecordingMax.xlsx');
verification_string = 'cats';

for m = 1:length(Monkeys)
    date_strs = {Monkeys(m).Sessions.DateStr};
    if regexp(Monkeys(m).Name, 'Marta\w*')
        short_names = get_short_names(date_strs, marta_xls, sprintf('\\w*%s\\w*', verification_string)); % need double \\ to get single \ in output
    else
        short_names = get_short_names(date_strs, max_xls, sprintf('\\w*%s\\w*', verification_string));
    end
    for i = 1:length(Monkeys(m).Sessions)
        Monkeys(m).Sessions(i).ShortName = short_names{i};
    end
end

%% Collect Session Y (image id and catg id)
% very fast

catgs = {1:260, 261:520};
for m = 1:length(Monkeys)
    sessions_to_use = 1:length(Monkeys(m).Sessions);
        
    for i = 1:length(sessions_to_use)
        sessn = sessions_to_use(i);
        catg1 = Monkeys(m).Sessions(sessn).CueInfo(catgs{1});
        catg2 = Monkeys(m).Sessions(sessn).CueInfo(catgs{2});
        catg1_timeson = vertcat(catg1.Times_on);
        catg2_timeson = vertcat(catg2.Times_on);

        % make Y for catg and image ID
        Y = vertcat(repelem(1,length(catg1_timeson))' , repelem(2,length(catg2_timeson))');
        Monkeys(m).Sessions(sessn).Session_Y_catg = Y;
        
        Y = [];
        for j = 1:length(Monkeys(m).Sessions(sessn).CueInfo)
            Y = vertcat(Y, repelem(Monkeys(m).Sessions(sessn).CueInfo(j).CueID, Monkeys(m).Sessions(sessn).CueInfo(j).NumApp)');
        end
        if length(Y) ~= length(Monkeys(m).Sessions(sessn).Session_Y_catg)
            error('Error compiling session Y''s ')
        end
        Monkeys(m).Sessions(sessn).Session_Y_imageID = Y;
    end
end

%% Count num "ignored" units
rSessionsByMonk = {[1 2 3 5 6 7 9], [1 2 3 4 5 6 7]};  % baseline days with reasonable trial counts
manualIntervals = {[175 275]};
ignoreVal = 20;

total_ignored = 0;
for m = 1:length(Monkeys)
    rSessions = rSessionsByMonk{m};
    
    for i = 1:length(rSessions)
        sessn = rSessions(i);
        fprintf("Monkey %s, session %d \n", Monkeys(m).Name, sessn)
        
        viable_bool = strcmp({Monkeys(m).Sessions(sessn).UnitInfo.Area}, "te");
        ignored_bool = [Monkeys(m).Sessions(sessn).UnitInfo.SpikeNum] <= ignoreVal;
        n_ignored = sum(viable_bool & ignored_bool);
        total_ignored = total_ignored + n_ignored;
        
        all_scs = [Monkeys(m).Sessions(sessn).UnitInfo(viable_bool).SpikeNum];
        disp(nanmin(all_scs(all_scs > 20)))
        
%         disp(n_ignored)
        
    end
%     disp(total_ignored)
end

%% Summarize firing rates across imgs

% Params
% rSessionsByMonk = {[1 2 3 5 6 7 9], [1 2 6 7]};  % baseline days with reasonable trial counts
% manualIntervals = {[-100 -0], [75 175], [175 275], [275 375]};

rSessionsByMonk = {[7 9], [6 7]};  % (Fig 2B / 2D)
manualIntervals = {[175 275]};

ignoreVal = 20;
for iInt = 1:length(manualIntervals)
interval = manualIntervals{iInt};

% figure
    for m = 1:length(Monkeys)
        rSessions = rSessionsByMonk{m};
%         subplot(2,1,m)
%         hold on

        for i = 1:length(rSessions)
            sessn = rSessions(i);

            % Get data from appropriate storage place
            % X_full is matrix of spike counts of dims (trials) x (units) x
            % (intervals), rIntervals is list of each spike-count interval.
            MonkID = sprintf('%s_%s', Monkeys(m).Name, Monkeys(m).Sessions(sessn).ShortName);
            fileName = sprintf('%s_allNeurons_step%d_wd%d.mat', MonkID, step, width);
            [X_full, rIntervals] = load_interval_data(fullfile(EXT_HD, spikeCountPath, fileName)); 
            Y_imgIDs = Monkeys(m).Sessions(sessn).Session_Y_imageID; % list of image ids for each presentation
            Y = Monkeys(m).Sessions(sessn).Session_Y_catg;
            rIntervals_original = rIntervals; % for finding idx in 3rd dim of X_full

            % Pre-allocate the storage vectors
            max_avg_rates_per_img = zeros(size(X_full,2), length(rIntervals));  % units x intervals x 1
            all_avg_rates_per_img = zeros(size(X_full,2), length(rIntervals), 520);  % units x intervals x image

            for iUnit = 1:size(X_full,2)

                % Catch units with too few spikes
                if Monkeys(m).Sessions(sessn).UnitInfo(iUnit).SpikeNum < ignoreVal
                    max_avg_rates_per_img(iUnit, :) = NaN;
                    continue
                end

                % Run the analysis
                idx = find(cellfun(@(a) all(a == interval), rIntervals_original));
                spike_counts = X_full(:, iUnit, idx);
                stim_means = splitapply(@mean, spike_counts, Y_imgIDs); % god bless MATLAB becoming more like python
                
                interval_duration_sec = (interval(2) - interval(1)) / 1000;
                fr = stim_means / interval_duration_sec;  % convert to spikes/sec
                max_avg_rates_per_img(iUnit, idx) = max(fr);
                all_avg_rates_per_img(iUnit, idx, :) = fr;
            end
            
            % Plot the full data
%             all_spike_counts = reshape(X_full(:,:,idx), 1, []);
%             rates = all_spike_counts * (1000/diff(interval));
%             histogram(rates, 'FaceColor', mlc(i), 'BinEdges', 0:2:100) 
            
            Monkeys(m).Sessions(sessn).Max_avg_rates = max_avg_rates_per_img;
            Monkeys(m).Sessions(sessn).All_avg_rates = all_avg_rates_per_img;
            fprintf("Done with Monkey %s, session %d \n", Monkeys(m).Name, sessn)
            
        end
    end
end

%% Plot pre/post distributions

% rSessionsByMonk = {[1 2 3 5 6 7 9], [1 2 6 7]};
rSessionsByMonk = {[7 9], [6 7]};

for iInt = 1:length(manualIntervals)
   interval = manualIntervals{iInt};
   f1 = figure;
   f2 = figure;
   int_idx = find(cellfun(@(a) all(a == interval), rIntervals_original));
   
    for m = 1:length(Monkeys)
        rSessions = rSessionsByMonk{m};
        fprintf("\n Int %d to %d, monk %s \n", interval(1), interval(2), Monkeys(m).Name)
        figure(f1)
        hold on
        subplot(2,1,m)
        hold on
        
        figure(f2)
        hold on
    
        for i = 1:length(rSessions)
            sessn = rSessions(i);
            max_avg_rates_per_img = Monkeys(m).Sessions(sessn).Max_avg_rates;
            these_max_rates = max_avg_rates_per_img(:, int_idx);
            
            % Plot data
            figure(f1)
            histogram(these_max_rates, 'Facecolor', mlc(i), 'BinEdges', 0:2:100, 'Normalization', 'probability')
            
            % Draw mean / std on plots
            mu = nanmedian(these_max_rates);
            sigma = iqr(these_max_rates);
            yl = ylim;
            marker_yval = yl(2) * 1.01 + 0.002 * (i-1);
            scatter(mu, marker_yval , 100, mlc(i), 'filled', 'v',  'MarkerEdgeColor', 'k')
            plot([mu-sigma mu+sigma], repelem(marker_yval, 2), '-', 'Color', mlc(i), 'LineWidth', 1.5)
            
            % Ditto but for all avg rates, split up by category
            all_avg_rates_per_img = Monkeys(m).Sessions(sessn).All_avg_rates;
            these_avg_rates = all_avg_rates_per_img(:, int_idx, :);
            figure(f2)
            
            subplot(2, 2, 2*(m-1) + mod(1-1, 2) + 1)
            hold on
            histogram(these_avg_rates(:, 1:260), 'BinEdges', 0:2:100, 'Normalization', 'probability')
            title(sprintf("Cats, %s", Monkeys(m).Name))
            set(gca, "YScale", "log")
            
            subplot(2, 2, 2*(m-1) + mod(2-1, 2) + 1)
            hold on
            histogram(these_avg_rates(:, 261:520), 'BinEdges', 0:2:100, 'Normalization', 'probability')
            title(sprintf("Dogs, %s", Monkeys(m).Name))
            set(gca, "YScale", "log")
        end
        
        % Format fig f1
        figure(f1)
%         title(sprintf("%s", Monkeys(m).Name), "Interpreter", "none")
        xlabel("Max img-avg firing rate")
        ylabel("Probability")
        xl = xlim;
        xlim([0, xl(2)])
%         legend(["Pre", "Post"])

        % Stats for f1
        [p, h] = ranksum(Monkeys(m).Sessions(rSessions(1)).Max_avg_rates(:,int_idx), ...
                         Monkeys(m).Sessions(rSessions(2)).Max_avg_rates(:,int_idx));
        fprintf("Maximum avg rates ranksum, p = %0.4g \n", p)
        disp([nanmedian(Monkeys(m).Sessions(rSessions(1)).Max_avg_rates(:,int_idx)), ....
            nanmedian(Monkeys(m).Sessions(rSessions(2)).Max_avg_rates(:,int_idx))])
        formatSVMPlot(gca, gcf)

        % Format fig f2
        figure(f2)
        xlabel("All img-avg firing rates")
        ylabel("Probability")
        xl = xlim;
        xlim([0, xl(2)])
        
        % Stats for f2
        cat_pre = reshape(Monkeys(m).Sessions(rSessions(1)).All_avg_rates(:,int_idx, 1:260), 1, []);
        cat_post = reshape(Monkeys(m).Sessions(rSessions(2)).All_avg_rates(:,int_idx, 1:260), 1, []);
%         [p, h] = ranksum(cat_pre, cat_post);
%         fprintf("All avg rates for CATS, ranksum, p = %0.4g \n", p)
%         disp([nanmean(cat_pre), nanmean(cat_post)])
%         disp([nanstd(cat_pre), nanstd(cat_post)])
%         
        dog_pre = reshape(Monkeys(m).Sessions(rSessions(1)).All_avg_rates(:,int_idx, 261:520), 1, []);
        dog_post = reshape(Monkeys(m).Sessions(rSessions(2)).All_avg_rates(:,int_idx, 261:520), 1, []);
%         [p, h] = ranksum(dog_pre, dog_post);
%         fprintf("All avg rates for DOGS, ranksum, p = %0.4g \n", p)
%         disp([nanmean(dog_pre), nanmean(dog_post)])
%         disp([nanstd(dog_pre), nanstd(dog_post)])

        figure
        hold on
        histogram(cat_pre - dog_pre, 'BinEdges', -50:50, 'Normalization', 'Probability');
        histogram(cat_post - dog_post, 'BinEdges', -50:50, 'Normalization', 'Probability');
        disp([median(cat_pre - dog_pre) median(cat_post - dog_post)])
        
        
        
    end
%     sgtitle(sprintf("Interval %d - %d", interval(1), interval(2)))
end

%% Plot disributions across baseline days

% rSessionsByMonk = {[7 9], [6 7]};  % (Fig 2B / 2D)
% manualIntervals = {[175 275]};

for m = 1:length(Monkeys)
    rSessions = rSessionsByMonk{m};
    for iInt = 1:length(manualIntervals)
        interval = manualIntervals{iInt};
        int_idx = find(cellfun(@(a) all(a == interval), rIntervals_original));
        figure
        hold on
        for i = 1:length(rSessions)
            sessn = rSessions(i);
            max_avg_rates_per_img = Monkeys(m).Sessions(sessn).Max_avg_rates;
            [f,x] = ecdf(max_avg_rates_per_img(:, int_idx ));
            plot(x,f, 'LineWidth', 2, 'DisplayName', Monkeys(m).Sessions(sessn).ShortName)
            
        end
        legend
        title(sprintf("%s, interval %d to %d", Monkeys(m).Name, interval(1), interval(2)), "Interpreter", "none")
        xlabel("Firing rate")
        ylabel("ECDF")
    end
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

