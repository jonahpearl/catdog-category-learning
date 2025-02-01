% figs for the cat/dog manuscript, 18 Oct 2023

%% Load data
clearvars
close all
EXT_HD = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/';
CCL = '/Users/jonahpearl/Documents/MATLAB/catdog-category-learning';
skip_path = 'SXMA/Monkey_structs';
load(fullfile(EXT_HD, skip_path, 'MaxMarta_skip_preProcessed.mat'), 'Monkeys');

%% Fig 1 %%
%===========%

%% Fig 1C, i,ii: %% Percent correct by sextile and day
% Outcome numerical code:
    % TP = 1
    % TN = 2
    % FP = 3
    % FN = 4
    % Err = 5.
% NB: percent correct out of completed trials. Does not count errors.

cb = cbrewer('qual', 'Set1', 9);
purps = cbrewer('seq', 'Purples', 10);
roygbiv = cb([1 5 6 3 2 4],:);
roygbiv(3,:) = roygbiv(3,:) + [-0.1 -0.2 -0.1];% dull the yellow
colors = vertcat(roygbiv, purps([8 10],:), [0 0 0]);

figure('Position', [300 300 500 800])
hold on
for m = 1:length(Monkeys)
    sessions_to_plot = 1:length(Monkeys(m).Sessions)-1;
%     sessions_to_plot = length(Monkeys(m).Sessions);  % transfer test is the last session in each struct here. would be better to look for "{'SkipUNI-01'}" as session short name.
    
    % set up figure
    subplot(2,1,m)
    hold on
    
    % pre allocate
    labs = {};
    pc_mat = zeros(6,length(sessions_to_plot));
    
    % get percent correct data
    for i = 1:length(sessions_to_plot)
        sessn = sessions_to_plot(i);
        outcomes = [Monkeys(m).Sessions(sessn).TrialInfo.OutcomeCode];
        completed = [Monkeys(m).Sessions(sessn).TrialInfo.Completed];
        len = length(outcomes);
        sixth = floor(len/6);
        for k = 1:6
            if k == 6
                ind = sixth*(k-1)+1;
                num_corr = sum(outcomes(ind:end) == 1 | outcomes(ind:end) == 2); 
                pc_mat(k,i) = num_corr / sum(completed(ind:end));
            else
                ind = [sixth*(k-1)+1 sixth*k];
                num_corr = sum(outcomes(ind(1):ind(2)) == 1 | outcomes(ind(1):ind(2)) == 2); 
                pc_mat(k,i) = num_corr / sum(completed(ind(1):ind(2)));
            end
        end
        labs = [labs sprintf('Day %d', i)];
    end
    
    % plot data
    for i = 1:length(sessions_to_plot)
        sessn = sessions_to_plot(i);
        
        % get color, trial-unique or roygbiv
        if sessn == length(Monkeys(m).Sessions)
            col_ind = length(colors);
        else
            col_ind = sessn;
        end
        
        % plot
        h = plot(1:6, pc_mat(:,i), '-o', 'MarkerSize', 8,...
            'DisplayName', Monkeys(m).Sessions(sessn).ShortName, 'Color', colors(col_ind,:),...
            'LineWidth', 1.5, 'MarkerSize', 15);
        set(h, {'MarkerFaceColor'}, {get(h,'Color')}); 
    end
    set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
        'XMinorTick', 'off', 'YMinorTick', 'off',...
        'fontsize',42, 'YGrid', 'on',...
        'fontname', 'Helvetica',...
        'XColor', 'black', 'YColor', 'black')
    xlabel('Session sextile')
    xticks(1:6)
    yticks(0.5:0.25:1)
    line([0 6], [0.5 0.5], 'Color', 'k', 'LineStyle', '--', 'HandleVisibility', 'off')
    ylim([0.45 1])
    ylabel('Fraction correct')
end
set(gcf,'color','w');

%% Figure 1C, iii: Percent correct over days

colors = vertcat(roygbiv, purps([8 10],:), [0 0 0]);

figure('Position', [300 300 240 400])
hold on
for m = 1:length(Monkeys) 

    % set up figure
    subplot(length(Monkeys), 1, m)
    hold on
    
    % get data
    pc_mat = zeros(1, length(Monkeys(m).Sessions));
    for i = 1:length(Monkeys(m).Sessions) 
        outcomes = [Monkeys(m).Sessions(i).TrialInfo.OutcomeCode];
        
        % percent correct
        num_comp_tr = sum([Monkeys(m).Sessions(i).TrialInfo.Completed]); 
        num_corr = sum(outcomes == 1 | outcomes == 2);
        pc_mat(i) = num_corr / num_comp_tr;
        
        % print whether statistically different from chance
%         chance_corr = ceil(num_comp_tr / 2);
%         if prop_test([num_corr chance_corr], [num_comp_tr num_comp_tr], 0, 0.05)
%             [h, pval, chisq, df] = prop_test([num_corr chance_corr], [num_comp_tr num_comp_tr], 0, 0.05);
%             fprintf('%s, session %s, %0.0f percent correct, performance signf. above chance (p = %d, chisq %0.5f, df %d \n',...
%                 Monkeys(m).Name, Monkeys(m).Sessions(i).ShortName, pc_mat(i)*100, pval, chisq, df)
%         end
        
        % get color, trial-unique or roygbiv
        if i == length(Monkeys(m).Sessions)
            col_ind = length(colors);
        else
            col_ind = i;
        end
        
        % scatter
        scatter(i, pc_mat(i), 500, colors(col_ind,:), 'filled');
    end        
    
    % plot 20/20 days
    p1 = plot(1:length(Monkeys(m).Sessions)-1, pc_mat(1:end-1),...
        'k-', 'LineWidth', 1.5, 'DisplayName', 'Percent Correct');
    uistack(p1,'bottom')
    
    % format plot
    xlabel('Session')
    xlim([0.5 length(Monkeys(m).Sessions)+0.5])
    xticks(1:length(Monkeys(m).Sessions))
    xticklabels([string(1:length(Monkeys(m).Sessions)-1) 'Tr.'])
    yticks(0.5:0.25:1)
    line([0 length(Monkeys(m).Sessions)], [0.5 0.5], 'Color', 'k', 'LineStyle', '--', 'HandleVisibility', 'off')
    ylim([0.45 1])
    ylabel('Fraction correct')
    set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
        'XMinorTick', 'off', 'YMinorTick', 'off',...
        'fontsize',42, 'YGrid', 'on',...
        'fontname', 'Helvetica',...
        'XColor', 'black', 'YColor', 'black')
end
set(gcf,'color','w');


%% Fig 1C,iv: Plot RT for each outcome type over sessions
% TP = 1
% TN = 2
% FP = 3
% FN = 4

outcome_names = {'TP', 'TN', 'FP', 'FN'};
outcomes_to_plot = 2;
cb = cbrewer('qual', 'Set1', 9);
purps = cbrewer('seq', 'Purples', 10);
roygbiv = cb([1 5 6 3 2 4],:);
roygbiv(3,:) = roygbiv(3,:) + [-0.1 -0.2 -0.1];% dull the yellow
colors = vertcat(roygbiv, purps([8 10],:), [0 0 0]);
figure('Position', [300 300 500 800])
hold on
for m = 1:length(Monkeys) 
    for o = 1:length(outcomes_to_plot)
        outcome = outcomes_to_plot(o);
        
        % set up figure
        subplot(length(Monkeys), length(outcomes_to_plot), length(outcomes_to_plot)*(m-1) + mod(o-1, length(outcomes_to_plot)) + 1)
        hold on
        
        RTs_mat = nan(5000,length(Monkeys(m).Sessions)); % (over-estimate num trials) x (num sessions)
        for i = 1:length(Monkeys(m).Sessions) 
            outcomes = [Monkeys(m).Sessions(i).TrialInfo.OutcomeCode];
            inds = find(outcomes == outcome);
            rts = [Monkeys(m).Sessions(i).TrialInfo(inds).RT];
            RTs_mat(1:length(rts),i) = rts;
            
            % get color, trial-unique or roygbiv
            if i == length(Monkeys(m).Sessions)
                col_ind = length(colors);
            else
                col_ind = i;
            end

            % plot individual data
            scatter(i, median(rts, 'omitnan'), 500, colors(col_ind,:), 'filled');
            errorbar(i, median(rts, 'omitnan'), iqr(rts), 'k.')
            
        end 
        
        % Plot a line through the training points.
        plot(1:size(RTs_mat,2)-1, median(RTs_mat(:,1:end-1),'omitnan'),...
            '-k', 'LineWidth', 1.5, 'MarkerSize', 6)

        % Format the plot.
        xlabel('Session')
        xlim([0.5 length(Monkeys(m).Sessions)+0.5])
        xticks(1:length(Monkeys(m).Sessions))
        xticklabels([string(1:length(Monkeys(m).Sessions)-1) 'Tr.'])
        ylabel('ms')
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
        'XMinorTick', 'off', 'YMinorTick', 'off',...
        'fontsize',42, 'YGrid', 'on',...
        'fontname', 'Helvetica',...
        'XColor', 'black', 'YColor', 'black')
    end
end
set(gcf,'color','w');

%% Colored text as legend for above plots
figure
hold on
axis off
grid off
for i = 1:size(colors,1)
    if i < 9
        text(0.5, 1-0.1*i, sprintf('Day %d', i), 'Color', colors(i,:), 'FontSize', 32, 'FontName', 'Helvetica')
    else
        text(0.5, 1-0.1*i, 'Trial Unique', 'Color', colors(i,:), 'FontSize', 32, 'FontName', 'Helvetica')
    end
end


