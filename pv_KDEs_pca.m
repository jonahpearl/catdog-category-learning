% look at dimensionality reduced trajectories

% also see: /Volumes/Alex\'s\ Mac\ Backup/Documents/MATLAB/matsumoto/jep_fix_cat_xma2_rasters_KDEs.m
%% Load data
clearvars
close all

% Define paths to data
EXT_HD = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/';
CCL = '/Users/jonahpearl/Documents/MATLAB/catdog-category-learning';
pv_path = 'XMA2/Monkey_structs';

%% Params
rSessions = {[7 9], [6 7]};

bw = 20; % arbitrary bandwidth for now
buffer = 2; % how many bw of buffer to give the KDE

% ms to use before and after cue
plot_window_back = 200;
plot_window_front = 500;
kde_x_vals = -(plot_window_back):plot_window_front; % exclude bw buffers now
            
spike_bounds_back = plot_window_back + bw*buffer;
spike_bounds_front = plot_window_front + bw*buffer;

%% Load raw data
% Load raw neural data
fname = 'MaxMarta_xma2_ni.mat';
Data = load(fullfile(EXT_HD, pv_path, fname));
[status, Monkeys] = stitch_monkeyStruct_from_parts(Data);
clear Data

%% Get KDEs for all trials (slow)

% TODO: save this k/l cell array for each unit, in a subfolder of monkey structs.

KDE = struct();
for m = 1:length(Monkeys)
    sessions_to_use = rSessions{m};
    KDE(m).Name = Monkeys(m).Name;
    
    for i = 1:length(sessions_to_use)
        session = sessions_to_use(i);
        KDE(m).Sessions(session).SessionName = Monkeys(m).Sessions(session).SessionName;
        KDE(m).Sessions(session).DateStr = Monkeys(m).Sessions(session).DateStr;
        
        units_to_use= find(ismember({Monkeys(m).Sessions(session).UnitInfo.Location}, {'anterior', 'middle', 'posterior'}));
        max_trials = max([Monkeys(m).Sessions(session).CueInfo.NumApp]);
        
        for j = 1:length(units_to_use) % for each unit
            unit = units_to_use(j);
            allspike_times = sort(Monkeys(m).Sessions(session).UnitInfo(unit).Spike_times); % sort once to save time with many binary searches later
            
            % Transfer some unit info
            KDE(m).Sessions(session).UnitInfo(unit).Electrode = Monkeys(m).Sessions(session).UnitInfo(unit).Electrode;
            KDE(m).Sessions(session).UnitInfo(unit).Electrode_unit = Monkeys(m).Sessions(session).UnitInfo(unit).Electrode_unit;
            KDE(m).Sessions(session).UnitInfo(unit).Location = Monkeys(m).Sessions(session).UnitInfo(unit).Location;
            
            % pre-allocate matrix for results
            all_trial_kdes = cell(length(Monkeys(m).Sessions(session).CueInfo), max_trials);
            
            for k = 1:length(Monkeys(m).Sessions(session).CueInfo) % for each cue
                times_on = Monkeys(m).Sessions(session).CueInfo(k).Times_on;
                
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
                    all_trial_kdes{k,l} = ksdensity(norm_times, kde_x_vals, 'Bandwidth', bw);
                end
            end
            KDE(m).Sessions(session).UnitInfo(unit).CueOnTrialByTrial.KDE_array = all_trial_kdes;
            KDE(m).Sessions(session).UnitInfo(unit).CueOnTrialByTrial.BW = bw;
            KDE(m).Sessions(session).UnitInfo(unit).CueOnTrialByTrial.Buffer = buffer;
            fprintf('Done with Monkey %d session %d unit %d \n', m, session, unit)
        end
    end
end

%% Save output

fname = fullfile(EXT_HD, pv_path, sprintf('MaxMarta_xma2_trialByTrial_KDEs_bw%d.mat', bw));

t = whos('KDE');
if t.bytes > 2e9 % split large struct and save all the vars
    fprintf('Splitting monkey struct to save\n')
    [status, vnames] = split_monkeyStruct_in_parts(KDE);
    save(fname, vnames{:}, '-v7.3');
else
    fprintf('Saving whole struct (no splitting)\n')
    save(fname, 'KDE');
end


%% PCA on KDEs (ie each unit of time is a dimension, find temporal templates)

rSessions = {[7 9], [6 7]};

pct = 99;
diffid = sprintf('TrueDiff_BW%d', bw);
bootid = sprintf('BootstrappedDiffs_BW%d', bw);

for m = 1:length(KDE)
    sessions_to_use = rSessions{m};
    
    for i = 1:length(sessions_to_use)
        session = sessions_to_use(i);
        units_to_plot = find(ismember({KDE(m).Sessions(session).UnitInfo.Location}, {'anterior', 'middle', 'posterior'}));
        
        % Collect data
        X = zeros(length(units_to_plot), 701);
        good_unit_bool = zeros(length(units_to_plot),1);
        for j = 1:length(units_to_plot)
            unit = units_to_plot(j);
            if  ~ismember('KDEYVals_raw', fields(KDE(m).Sessions(session).UnitInfo(unit).CueOnAllCues))
                continue
            end
            good_unit_bool(j) = 1;
            kde_x_vals = KDE(m).Sessions(session).UnitInfo(unit).CueOnAllCues.KDEXVals;
            X(j,:) = KDE(m).Sessions(session).UnitInfo(unit).CueOnAllCues.KDEYVals_raw;
        end
        
        X = X(good_unit_bool==1, :);
        
        % Run pca
        X_standardized = (X - mean(X,1)) ./ (std(X,1));
        [components, X_new, ~, ~, explained_var] = pca(X_standardized);
        
        break
       
    end
    break
end

    
    
    
    