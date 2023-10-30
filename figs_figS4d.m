% vogels figs 

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

%% Load sparseness data

data = load(fullfile(EXT_HD, pv_path, 'MaxMarta_xma2_sparseness.mat'), 'Monkeys');
for m = 1:length(data.Monkeys)
    for i = 1:length(data.Monkeys(m).Sessions)
        Monkeys(m).Sessions(i).Unit_sparseness = data.Monkeys(m).Sessions(i).Unit_sparseness;
        Monkeys(m).Sessions(i).Unit_sparseness_dog = data.Monkeys(m).Sessions(i).Unit_sparseness_dog;
        Monkeys(m).Sessions(i).Unit_sparseness_cat = data.Monkeys(m).Sessions(i).Unit_sparseness_cat;
%         Monkeys(m).Sessions(i).Unit_sparseness_SHUFFLE = data.Monkeys(m).Sessions(i).Unit_sparseness_SHUFFLE;
    end
end
clear data

%% Load Nst data

data = load(fullfile(EXT_HD, pv_path, 'MaxMarta_xma2_Nst_sharedMax.mat'), 'Monkeys');
for m = 1:length(data.Monkeys)
    for i = 1:length(data.Monkeys(m).Sessions)
        Monkeys(m).Sessions(i).Nst = data.Monkeys(m).Sessions(i).Nst;
    end
end
clear data

%% Fig S4b: params
rSessionsByMonk = {[7 9] [6 7]};
ranksum_alpha = 0.05;
rArrayLocs = {'te'};

nst_iInt_to_plot = 2;  % intervals are [75 175], [175 275], [275 375]

sparseness_iInt_to_plot = 1;  %  % 1,2 are 75-175, 175-275

%% Plot Nst all (without catg)

figure

for m = 1:length(Monkeys)

rSessions = rSessionsByMonk{m};
        
    for iLoc = 1:length(rArrayLocs)
        loc = rArrayLocs{iLoc};
        
        for i = 1:length(rSessions)
            sessn = rSessions(i);
            
            % Get indices of units in the pValues matrix to use in calculating
            % proprotion of units with signf GLMs
            if strcmp(loc, 'te')
                units = strcmp({Monkeys(m).Sessions(sessn).UnitInfo.Area}, loc);
            else
                units = strcmp({Monkeys(m).Sessions(sessn).UnitInfo.Location}, loc);
            end
            
            
            % Get data, store for stats
            nst = Monkeys(m).Sessions(sessn).Nst.Nst_all(units, 2)/520;
            if regexp(Monkeys(m).Sessions(sessn).ShortName, 'Pre.*')
               pre =  nst; 
               color = mlc(1);
               lw = 1.5;
            elseif regexp(Monkeys(m).Sessions(sessn).ShortName, 'Post.*')
                post = nst;
                color = mlc(2);
                lw = 1.5;
            else
                color = [0.4 0.4 0.4];
                lw = 0.75;
            end
            
            % Plot data
            
            % One subplot per array, comparing sessions
            subplot(2, length(rArrayLocs), sub2ind([length(rArrayLocs), 2], iLoc, m))
            hold on
            [f,x] = ecdf(nst);
%             plot(x,f, 'DisplayName', Monkeys(m).Sessions(sessn).ShortName, 'Color', color, 'LineWidth', lw)
            plot(x,f, 'Color', color, 'LineWidth', lw)
            title(sprintf('%s, %s', Monkeys(m).Name, loc))
%             legend({'Pre', 'Post'})
            xlim([0 1])
        end

        formatPlot(gca, gcf)
%         legend
        
        % Stats test pre vs post: test for broadening of tuning curves
        [pval, h] = ranksum(pre, post, 'alpha', ranksum_alpha, 'tail', 'left');
        fprintf('%s, %s, Nst pre vs post (ranksum): h = %d, p = %0.2f \n',...
            Monkeys(m).Name, loc, h, pval)
          
    end
end

%% Nst cat vs Nst dog

rSessionsByMonk = {[7 9], [6 7]};
 
figure('Position', [400 400 400 400])
for m = 1:length(Monkeys)

rSessions = rSessionsByMonk{m};
        
    for iLoc = 1:length(rArrayLocs)
        loc = rArrayLocs{iLoc};
        
        for i = 1:length(rSessions)
            sessn = rSessions(i);
            
            % Get indices of units in the pValues matrix to use in calculating
            % proprotion of units with signf GLMs
            if strcmp(loc, 'te')
                units = strcmp({Monkeys(m).Sessions(sessn).UnitInfo.Area}, loc);
            else
                units = strcmp({Monkeys(m).Sessions(sessn).UnitInfo.Location}, loc);
            end
            
            if regexp(Monkeys(m).Sessions(sessn).ShortName, 'Pre.*')  
               lstyle = '--';
            elseif regexp(Monkeys(m).Sessions(sessn).ShortName, 'Post.*')
                lstyle = '-';
            else
                lstyle = '-';
            end
            
            % Get data
            cats = Monkeys(m).Sessions(sessn).Nst.Nst_catdog(units, nst_iInt_to_plot, 1)/260;
            dogs = Monkeys(m).Sessions(sessn).Nst.Nst_catdog(units, nst_iInt_to_plot, 2)/260;
            overall = Monkeys(m).Sessions(sessn).Nst.Nst_all(units, nst_iInt_to_plot)/520;
            
            
            subplot(2,1,m)
            hold on
            [f,x] = ecdf(cats);
            plot(x,f, 'Color', 'red', 'LineStyle', lstyle)
            [f,x] = ecdf(dogs);
            plot(x,f, 'Color', 'blue','LineStyle', lstyle)
            
            [f,x] = ecdf(overall);
            plot(x,f, 'Color', 'k', 'LineStyle', lstyle)
            
            
            title(sprintf('%s, %s', Monkeys(m).Name, loc));
%             legend({'Pre', 'Post'})
            xlim([0 1]);
            formatPlot(gca, gcf);
            axis square
            
            % Plot data
%             figure(3)
%             subplot(2,2, sub2ind([2,2], i, m))
%             hold on
%             title(sprintf('%s, session %s', Monkeys(m).Name, Monkeys(m).Sessions(sessn).ShortName))
%             scatter(dogs, cats, 'filled')
%             legend(rArrayLocs)
%             xlim([0 1])
%             ylim([0 1])
%             
%             figure(4)
%             subplot(2,3, sub2ind([3,2], iLoc, m))
%             hold on
%             scatter(dogs, cats, 'filled')
%             title(sprintf('%s, %s', Monkeys(m).Name, loc))
%             legend({'Pre', 'Post'})
%             xlim([0 1])
%             ylim([0 1])
        end
    end
end

%% Plot sparseness over sessions

% sparseness_iInt_to_plot = 2;  % 1,2,3 are 75-175, 175-275, 275-375
figure
show_shuffle = false;
for m = 1:length(Monkeys)

rSessions = rSessionsByMonk{m};
        
    for iLoc = 1:length(rArrayLocs)
        loc = rArrayLocs{iLoc};
        
        for i = 1:length(rSessions)
            sessn = rSessions(i);
            
            % Get indices of units in the pValues matrix to use in calculating
            % proprotion of units with signf GLMs
            if strcmp(loc, 'te')
                units = ismember({Monkeys(m).Sessions(sessn).UnitInfo.Location}, TE_LOCS);
            else
                units = ismember({Monkeys(m).Sessions(sessn).UnitInfo.Location}, loc);
            end
            
            sparsenesss_indexes = Monkeys(m).Sessions(sessn).Unit_sparseness(units, sparseness_iInt_to_plot); 
            if show_shuffle
                sparsenesss_indexes_shuff = Monkeys(m).Sessions(sessn).Unit_sparseness_SHUFFLE(units, sparseness_iInt_to_plot, :); 
            end
            if regexp(Monkeys(m).Sessions(sessn).ShortName, 'Pre.*')
               pre = sparsenesss_indexes; 
               color = mlc(1);
               lw = 1.5;
            elseif regexp(Monkeys(m).Sessions(sessn).ShortName, 'Post.*')
                post = sparsenesss_indexes;
                color = mlc(2);
                lw = 1.5;
            else
                color = [0.4 0.4 0.4];
                lw = 0.75;
            end
            
            
            % Prep the ecdf plot
            [f, x] = ecdf(sparsenesss_indexes);
            
            % Make the plot
            subplot(2, length(rArrayLocs), sub2ind([length(rArrayLocs), 2], iLoc, m))
            hold on
            
%             plot(x, f, 'Color', mlc(i), 'DisplayName', Monkeys(m).Sessions(sessn).ShortName)
            plot(x, f, 'Color', color, 'LineWidth', lw)
            title(sprintf('%s, %s', Monkeys(m).Name, loc))
            xlim([0 1])
%             xlabel('"Sparseness index" (inverse sparseness)')
            xlabel('s')
            
            % Show shuffle if requested
            if show_shuffle
                for iShuff=1:nShuffles
                    [f_shuff, x_shuff] = ecdf(sparsenesss_indexes_shuff(:,:,iShuff));
                    plot(x_shuff, f_shuff, 'Color', color, 'LineWidth', lw-1, 'LineStyle', '--')
                end
            end
            
        end
        % Stats test pre vs post: test for broadening of tuning curves (ie
        % greater sparseness index [i know, it's confusing that going up
        % means less sparseness...that's just how vogels defined it]).
        [pval, h] = ranksum(pre, post, 'alpha', ranksum_alpha, 'tail', 'left');
        fprintf('%s, %s, Nst pre vs post (ranksum): h = %d, p = %0.2f \n',...
        Monkeys(m).Name, loc, h, pval)
        formatPlot(gca, gcf)
    end
end

%% Plot sparseness within sessions, cat vs dog

sparseness_iInt_to_plot = 2;  % 1,2 are 75-175, 175-275
% sparseness_iInt_to_plot = 1;

figure

for m = 1:length(Monkeys)

rSessions = rSessionsByMonk{m};
        
    for iLoc = 1:length(rArrayLocs)
        loc = rArrayLocs{iLoc};
        
        for i = 1:length(rSessions)
            sessn = rSessions(i);
            
            % Get indices of units in the pValues matrix to use in calculating
            % proprotion of units with signf GLMs
            if strcmp(loc, 'te')
                units = ismember({Monkeys(m).Sessions(sessn).UnitInfo.Location}, TE_LOCS);
            else
                units = ismember({Monkeys(m).Sessions(sessn).UnitInfo.Location}, loc);
            end
            
            sparsenesss_indexes = Monkeys(m).Sessions(sessn).Unit_sparseness(units, sparseness_iInt_to_plot); 
            sparsenesss_indexes_cat = Monkeys(m).Sessions(sessn).Unit_sparseness_cat(units, sparseness_iInt_to_plot); 
            sparsenesss_indexes_dog = Monkeys(m).Sessions(sessn).Unit_sparseness_dog(units, sparseness_iInt_to_plot); 
            
            if regexp(Monkeys(m).Sessions(sessn).ShortName, 'Pre.*')
               pre = sparsenesss_indexes; 
               color = mlc(1);
               lw = 1.5;
               lstyle = '--';
            elseif regexp(Monkeys(m).Sessions(sessn).ShortName, 'Post.*')
                post = sparsenesss_indexes;
                color = mlc(2);
                lw = 1.5;
                lstyle = '-';
            else
                color = [0.4 0.4 0.4];
                lw = 0.75;
            end
            
            
            % Make the plot
%             subplot(2, length(rSessions), sub2ind([length(rSessions), 2], i, m))
            subplot(2,1,m)
            hold on
            
            % Plot overall popn sparseness
            [f, x] = ecdf(sparsenesss_indexes);
            plot(x, f, 'Color', 'k', 'LineWidth', lw, 'DisplayName', 'overall', 'LineStyle', lstyle)
            
            % Category-specific
            [f, x] = ecdf(sparsenesss_indexes_cat);
            plot(x, f, 'Color', 'red', 'LineWidth', lw, 'DisplayName', 'overall', 'LineStyle', lstyle)
            [f, x] = ecdf(sparsenesss_indexes_dog);
            plot(x, f, 'Color', 'blue', 'LineWidth', lw, 'DisplayName', 'overall', 'LineStyle', lstyle)
            
            % Test for catg differences
            [pval, h] = ranksum(sparsenesss_indexes_cat, sparsenesss_indexes_dog,...
                'alpha', ranksum_alpha, 'tail', 'left');  % test dogs less sparse --> dogs > cats --> Y > X --> left tail
            fprintf('%s, %s, %s, sparseness dog > cat (ranksum): h = %d, p = %0.2f \n',...
                Monkeys(m).Name, loc, Monkeys(m).Sessions(sessn).ShortName, h, pval)
            
            title(sprintf('%s, %s', Monkeys(m).Name, loc),...
                'Interpreter', 'none')
            xlim([0 1])
%             xlabel('"Sparseness index" (inverse sparseness)')
            xlabel('s')
            formatPlot(gca, gcf)
            
        end
        % Stats test pre vs post: test for broadening of tuning curves (ie
        % greater sparseness index [i know, it's confusing that going up
        % means less sparseness...that's just how vogels defined it]).
        [pval, h] = ranksum(pre, post, 'alpha', ranksum_alpha, 'tail', 'left');
        fprintf('%s, %s, overall sparseness pre vs post (ranksum): h = %d, p = %0.2f \n',...
            Monkeys(m).Name, loc, h, pval)
    end
end

%% Functions
function formatPlot(ax, fig)
set(ax, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off',...
    'fontsize',20, 'YGrid', 'on', 'XGrid', 'on',...
    'fontname', 'Helvetica',...
    'XColor', 'black', 'YColor', 'black')
set(fig, 'Color', 'white')
end