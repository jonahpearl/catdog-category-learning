% figs for the cat/dog manuscript, 18 Oct 2023

%% Fig 6 A, B: run pv_VR_scatter_plot.m

%% Fig 6C: load data

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

%% Fig 6C: plot sparseness within sessions, cat vs dog
rSessionsByMonk = {[7 9] [6 7]};
rArrayLocs = {'te'};
ranksum_alpha = 0.05;

idx_in_sparseness_mat = 2;  % 1,2 are 75-175, 175-275
% idx_in_sparseness_mat = 1;

close all
figure(1)

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
            
            sparsenesss_indexes = Monkeys(m).Sessions(sessn).Unit_sparseness(units, idx_in_sparseness_mat); 
            sparsenesss_indexes_cat = Monkeys(m).Sessions(sessn).Unit_sparseness_cat(units, idx_in_sparseness_mat); 
            sparsenesss_indexes_dog = Monkeys(m).Sessions(sessn).Unit_sparseness_dog(units, idx_in_sparseness_mat); 
            
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
            figure(1)
%             subplot(2, length(rSessions), sub2ind([length(rSessions), 2], i, m))
            subplot(2,1,m)
            hold on
            
            % Plot overall popn sparseness
%             [f, x] = ecdf(sparsenesss_indexes);
%             plot(x, f, 'Color', 'k', 'LineWidth', lw, 'DisplayName', 'overall', 'LineStyle', lstyle)
            
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
            
            [pval, h] = ranksum(sparsenesss_indexes_cat, sparsenesss_indexes_dog,...
                'alpha', ranksum_alpha);  % test two sided
            fprintf('%s, %s, %s, sparseness dog vs cat (ranksum): h = %d, p = %0.2f \n',...
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
        
        fprintf("\n")
    end
end

%% Fig 6D-F: run pv_population_geometry.m

%% Functions
function formatPlot(ax, fig)
set(ax, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off',...
    'fontsize',20, 'YGrid', 'on', 'XGrid', 'on',...
    'fontname', 'Helvetica',...
    'XColor', 'black', 'YColor', 'black')
set(fig, 'Color', 'white')
end
