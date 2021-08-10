% make KDEs for each unit / image

% %% Load in saved data
% clearvars
% close all
% 
% % load Max data
% monkey_names = {'Max_fix_cat_xma2'};
% path1 = '/Users/pearlje/Documents/MATLAB/matsumoto/XMA2/Monkey_structs';
% fname = strcat(monkey_names{1},'_ni','.mat');
% Data = load(fullfile(path1, fname)); 
% [status, Monkeys] = stitch_monkeyStruct_from_parts(Data);
% clear Data
% 
% 
% % load Marta data
% monkey_names = {'Marta_fix_cat_xma2'};
% path1 = '/Users/pearlje/Documents/MATLAB/matsumoto/XMA2/Monkey_structs';
% fname = strcat(monkey_names{1},'_ni','.mat');
% Data = load(fullfile(path1, fname));
% Monkeys(2).Name = Data.Monkeys.Name;
% Monkeys(2).Sessions = Data.Monkeys.Sessions;
% clear Data

%% Load data
clearvars
close all

monkey_names = {'MaxMarta_xma2'}; % MaxMartaMarta_fix_cat_xma2, Max_CARTRUCK_fix_cat_xma2
path1 = '/Users/pearlje/Documents/MATLAB/matsumoto/XMA2/Monkey_structs';
load_waveforms = false;
if load_waveforms
    fname = strcat(monkey_names{1},'_ni_wf','.mat'); 
else
    fname = strcat(monkey_names{1},'_ni','.mat');
end


Data = load(fullfile(path1, fname)); 
[status, Monkeys] = stitch_monkeyStruct_from_parts(Data);
clear Data

%% ANALYSIS -- calculate mean KDEs for each neuron

% We'll go from 200 ms before cue on to 500 ms after cue on.
% I plan to use this KDE analysis for a few things
% 1. We can normalize based on it, and begin to correlate neural responses
% 2. We can cross-correlate them, and find neurons with similar patterns of
% response (this could eventually lead to a clustering analysis, but CC
% seems less biased for now).

bw = 20; % arbitrary bandwidth for now
buffer = 2; % how many bw of buffer to give the KDE

% ms to use before and after cue
plot_window_back = 200;
plot_window_front = 500;

spike_bounds_back = plot_window_back + bw*buffer;
spike_bounds_front = plot_window_front + bw*buffer;
% same strategy as "Each unit aligned to cue on, all cues x all trials"
% code above.

categories = {1:260, 261:520};


KDE = struct();
for m = 1:length(Monkeys)
% for m = 1
    sessions_to_plot = 1:length(Monkeys(m).Sessions);
%     sessions_to_plot = [7 9];
    KDE(m).Name = Monkeys(m).Name;
    
    
    for i = 1:length(sessions_to_plot)
        sessn = sessions_to_plot(i);
        KDE(m).Sessions(sessn).SessionName = Monkeys(m).Sessions(sessn).SessionName;
        KDE(m).Sessions(sessn).DateStr = Monkeys(m).Sessions(sessn).DateStr;
        
        units_to_plot = 1:length(Monkeys(m).Sessions(sessn).UnitInfo);
        
        for j = 1:length(units_to_plot) % for each unit
            unit = units_to_plot(j);
            max_trials = max([Monkeys(m).Sessions(sessn).CueInfo.NumApp]);
            allspike_times = sort(Monkeys(m).Sessions(sessn).UnitInfo(unit).Spike_times); % sort once to save time with many binary searches later
            
            KDE(m).Sessions(sessn).UnitInfo(unit).Electrode = Monkeys(m).Sessions(sessn).UnitInfo(unit).Electrode;
            KDE(m).Sessions(sessn).UnitInfo(unit).Electrode_unit = Monkeys(m).Sessions(sessn).UnitInfo(unit).Electrode_unit;
            KDE(m).Sessions(sessn).UnitInfo(unit).Location = Monkeys(m).Sessions(sessn).UnitInfo(unit).Location;
            
            % preallocate rough size of arrays for plotting
            % say avg rate is 25 sp/s = 0.025 sp / ms
            avg_rate = 0.25*(spike_bounds_back + spike_bounds_front);
            x_vals = zeros(max_trials*260*avg_rate, 1); % (trials x cues) = row, row x num_spikes / row = num_spikes
            y_vals = zeros(max_trials*260*avg_rate, 1);
            catg_vals = zeros(max_trials*260*avg_rate, 1);
            xx = 1; % indexer
            
            % collect data
            for k = 1:length(Monkeys(m).Sessions(sessn).CueInfo) % for each cue
                
                % NB, this times_on field ONLY contains timestamps from
                % within completed trials (see behav_analysis code). If you want to include fixation
                % errors, either make a new field in the cueInfo struct
                % called err_times_on or something.
                % Or if you just want fixation errs, get time and cue from the
                % trialInfo struct for all trials that were fixation
                % errors.
                times_on = Monkeys(m).Sessions(sessn).CueInfo(k).Times_on;
                for q = 1:2
                    if ismember(Monkeys(m).Sessions(sessn).CueInfo(k).CueID, categories{q})
                        catg = q;
                        break
                    end
                end
                    
                bounds = [(times_on - spike_bounds_back) (times_on + spike_bounds_front)];
                for l = 1:length(bounds) % for each time that cue appeared
                    [lower_idx, upper_idx] = binarySearch_window(allspike_times, bounds(l, 1), bounds(l, 2)); % allspike_times must be sorted!
                    
                    % check edge cases
                    if lower_idx == -1 || upper_idx == -1 || upper_idx - lower_idx == -1
                        % no spikes in window
                        % lower_idx is -1 --> bounds are higher than
                        % all spike times
                        % uppder idx is -1 --> bounds are lower than spike
                        % times
                        % diff is -1 --> bounds are in between two spike
                        % times
                        continue
                    end
                    trial_spike_times = allspike_times(lower_idx:upper_idx);
                    norm_times = trial_spike_times - times_on(l);
                    x_vals(xx:xx + length(norm_times) - 1) = norm_times;
                    y_vals(xx:xx + length(norm_times) - 1) = repelem(k + 0.01*l, length(norm_times));
                    catg_vals(xx:xx + length(norm_times) - 1) = repelem(catg, length(norm_times));
                    xx = xx + length(norm_times);
                end
            end
            
            % remove extra zeros
            x_vals(xx:end) = [];
            y_vals(xx:end) = [];
            catg_vals(xx:end) = [];
            
            if isempty(x_vals)
                KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.RasterXVals = x_vals;
                KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.RasterYVals = y_vals;
                KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEXVals = kde_x_vals;
                KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEYVals = NaN;
                KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.BW = bw;
                KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.Buffer = buffer;
                continue
            end
            
            % run KDE
            kde_x_vals = -(plot_window_back):plot_window_front; % exclude bw buffers now
            f = ksdensity(x_vals, kde_x_vals, 'Bandwidth', bw);
            
            % run catg specific KDEs. use all points to create KDE
            if sum(catg_vals==1) == 0
                f_catg1 = zeros(1, numel(kde_x_vals));
                f_catg2 = ksdensity(x_vals(catg_vals==2), kde_x_vals, 'Bandwidth', bw);
            elseif sum(catg_vals==2) == 0
                f_catg1 = ksdensity(x_vals(catg_vals==1), kde_x_vals, 'Bandwidth', bw);
                f_catg2 = zeros(1, numel(kde_x_vals));
            else
                f_catg1 = ksdensity(x_vals(catg_vals==1), kde_x_vals, 'Bandwidth', bw);
                f_catg2 = ksdensity(x_vals(catg_vals==2), kde_x_vals, 'Bandwidth', bw);
            end


            % store
            KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.RasterXVals = x_vals;
            KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.RasterYVals = y_vals;
            KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.CatgVals = catg_vals;
            KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEXVals = kde_x_vals;
            KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEYVals_raw = f;
            KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEYVals_catg1_raw = f_catg1;
            KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEYVals_catg2_raw = f_catg2;
            KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.KDEYVals_SRs = [xx sum(catg_vals==1) sum(catg_vals==2)]; % [total, cats, dogs]
            KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.BW = bw;
            KDE(m).Sessions(sessn).UnitInfo(unit).CueOnAllCues.Buffer = buffer;
            
        end
    end
end

%% Save the analyzed data

path1 = '/Users/pearlje/Documents/MATLAB/matsumoto/XMA2/Monkey_structs';
fname = sprintf('MaxMarta_xma2_dualKDEs_bw%d.mat', bw);

t = whos('KDE');
if t.bytes > 2e9 % split large struct and save all the vars
    fprintf('Splitting monkey struct to save\n')
    [status, vnames] = split_monkeyStruct_in_parts(KDE);
    save(fullfile(path1, fname), vnames{:});
else
    fprintf('Saving whole struct (no splitting)\n')
    save(fullfile(path1, fname), 'KDE');
end
