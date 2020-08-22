% analyze category behavior, without any neural data


%% Load data
clearvars
close all
EXT_HD = '/Volumes/Alex''s Mac Backup/Documents/MATLAB/matsumoto/';
CCL = '/Users/jonahpearl/Documents/MATLAB/catdog-category-learning';
skip_path = 'SXMA/Monkey_structs';
load(fullfile(EXT_HD, skip_path, 'MaxMarta_skip_preProcessed.mat'), 'Monkeys');

%% Plot RT for each outcome type over sessions
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
        
        
        % Plot 20/20 days
%         errorbar(1:size(RTs_mat,2)-1, median(RTs_mat(:,1:end-1),'omitnan'), iqr(RTs_mat(:,1:end-1)),...
%             '-o', 'LineWidth', 1.5, 'Color', colors(outcome,:), 'MarkerSize', 6,...
%             'MarkerFaceColor', colors(outcome,:))

        % Plot a line through the training points.
        plot(1:size(RTs_mat,2)-1, median(RTs_mat(:,1:end-1),'omitnan'),...
            '-k', 'LineWidth', 1.5, 'MarkerSize', 6)
        
        
        
        % Plot trial unique day.
%         errorbar(size(RTs_mat,2), median(RTs_mat(:,end),'omitnan'), iqr(RTs_mat(:,end)),...
%             '-o', 'LineWidth', 1.5, 'Color', 'k', 'MarkerSize', 6,...
%             'MarkerFaceColor', 'k')



        % Format the plot.
%         title(outcome_names(outcomes_to_plot(o)))
        xlabel('Session')
        xlim([0.5 length(Monkeys(m).Sessions)+0.5])
        xticks(1:length(Monkeys(m).Sessions))
        xticklabels([string(1:length(Monkeys(m).Sessions)-1) 'Tr.'])
%         xtickangle(45)
        ylabel('ms')
%         title(sprintf('%s, %s', Monkeys(m).Name, outcome_names{outcome}), 'Interpreter', 'none')
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
        'XMinorTick', 'off', 'YMinorTick', 'off',...
        'fontsize',42, 'YGrid', 'on',...
        'fontname', 'Helvetica',...
        'XColor', 'black', 'YColor', 'black')
    end
end
set(gcf,'color','w');

%% Percent correct by sextile and day

% This should improve with learning.
% NHP should ROR for catg 1 and ROG for
% catg2. False positive, ie ROG for catg1, leads to timeout.
% True positive, ie ROG for catg 2, leads to reward. True and
% false negatives have no effect, just moves to next trial.
% Outcome numerical code:
    % TP = 1
    % TN = 2
    % FP = 3
    % FN = 4
    % Err = 5.
    
% NB: percent correct out of completed trials. Does not count errors.


% session by sextile, color pc
% figure
% pc_mat_flip = flipud(pc_mat); % make the first sextile be in the bottom left, ie row 6
% imagesc(pc_mat_flip)
% c = colorbar();
% c.Label.String = 'Percent correct';
% xticklabels(labs)
% yticklabels(6:-1:1)
% xlabel('Session date')
% ylabel('Session sextile')
% title('Marta SXMA, percent correct over session-time and days')
% set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
%     'XMinorTick', 'off', 'YMinorTick', 'off',...
%     'LineWidth', 1,  'fontsize',20, ...
%     'fontname', 'Helvetica','FontWeight','Bold', 'LineWidth', 2,...
%     'XColor', 'black', 'YColor', 'black')
% set(gcf,'color','w');


% ====== Build-up (slow roll) plots for a talk ======= %
% each session is a line
% figure
% hold on
% cb = cbrewer('qual', 'Set1', 9);
% purps = cbrewer('seq', 'Purples', 10);
% roygbiv = cb([1 5 6 3 2 4],:);
% roygbiv(3,:) = roygbiv(3,:) + [-0.1 -0.2 -0.1];% dull the yellow
% colors = vertcat(roygbiv, purps([8 10],:));
% xlabel('Sextile')
% ylabel('Percent correct')
% % title('Marta SXMA, percent correct over sextiles')
% 
% sessions_to_plot_arr = {1, 1:2, 1:3, 1:4, 1:5, 1:6, 1:8};
% for k = 1:length(sessions_to_plot_arr)
%     figure('Position', [1000  974 809 364])
%     hold on
%     sessions_to_plot = sessions_to_plot_arr{k};
%     for m = 1:length(Monkeys)
%         for i = 1:length(sessions_to_plot)
%             sessn = sessions_to_plot(i);
%             h = plot(1:6, pc_mat(:,sessn), '-o', 'MarkerSize', 15,...
%                 'DisplayName', labs{sessn}, 'Color', colors(sessn,:),...
%                 'LineWidth', 3);
%             set(h, {'MarkerFaceColor'}, {get(h,'Color')}); 
%         end
%     end
%     legend('Location', 'eastoutside')
%     set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
%         'XMinorTick', 'off', 'YMinorTick', 'off',...
%         'fontsize',30, 'YGrid', 'on',...
%         'fontname', 'Helvetica','FontWeight','Bold', 'LineWidth', 2,...
%         'XColor', 'black', 'YColor', 'black')
%     set(gcf,'color','w');
%     xlabel('Sextile')
%     xticks(1:6)
%     yticks(0.4:0.2:1)
%     line([0 6], [0.5 0.5], 'Color', 'k', 'HandleVisibility', 'off', 'LineWidth', 2)
%     yticklabels(0.4:0.2:1)
%     ylim([0.4 1])
%     ylabel('Percent correct')
%     
%     saveas(gcf, sprintf('sextiles_up_to_%d.png', k))
% end


% ======== Plot all data in one figure ==========%

cb = cbrewer('qual', 'Set1', 9);
purps = cbrewer('seq', 'Purples', 10);
roygbiv = cb([1 5 6 3 2 4],:);
roygbiv(3,:) = roygbiv(3,:) + [-0.1 -0.2 -0.1];% dull the yellow
colors = vertcat(roygbiv, purps([8 10],:), [0 0 0]);

% figure('Position', [300 300 240 400])
figure('Position', [300 300 500 800])
hold on
for m = 1:length(Monkeys)
    sessions_to_plot = 1:length(Monkeys(m).Sessions)-1;
%     sessions_to_plot = length(Monkeys(m).Sessions);
%     sessions_to_plot = 1:length(Monkeys(m).Sessions);
    
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
%         h = plot3(1:6, repelem(i, 1,6), pc_mat(:,sessn), '-o', 'MarkerSize', 8,...
%             'DisplayName', Monkeys(m).Sessions(sessn).ShortName, 'Color', colors(col_ind,:),...
%             'LineWidth', 1.5);
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
%     legend
end
set(gcf,'color','w');
% sgtitle('Cat/dog learning')

%% Percent correct over days


% colors =  [0    0.4470    0.7410; % the matlab colors!
%     0.8500    0.3250    0.0980;
%     0.9290    0.6940    0.1250;
%     0.4940    0.1840    0.5560;
%     0.4660    0.6740    0.1880;
%     0.3010    0.7450    0.9330;
%     0.6350    0.0780    0.1840;
%     0.3          0.3        0.3    ;
%     0          0        0    ];
colors = vertcat(roygbiv, purps([8 10],:), [0 0 0]);

figure('Position', [300 300 240 400])
% figure('Position', [300 300 500 800])
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
        chance_corr = ceil(num_comp_tr / 2);
        if prop_test([num_corr chance_corr], [num_comp_tr num_comp_tr], 0, 0.05)
            [h, pval, chisq, df] = prop_test([num_corr chance_corr], [num_comp_tr num_comp_tr], 0, 0.05);
            fprintf('%s, session %s, %0.0f percent correct, performance signf. above chance (p = %d, chisq %0.5f, df %d \n',...
                Monkeys(m).Name, Monkeys(m).Sessions(i).ShortName, pc_mat(i)*100, pval, chisq, df)
        end
        
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
    
    % plot UNI day
%     plot(length(Monkeys(m).Sessions), pc_mat(end),...
%         'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Percent Correct')
    
    % format plot
    xlabel('Session')
    xlim([0.5 length(Monkeys(m).Sessions)+0.5])
    xticks(1:length(Monkeys(m).Sessions))
    xticklabels([string(1:length(Monkeys(m).Sessions)-1) 'Tr.'])
    yticks(0.5:0.25:1)
    line([0 length(Monkeys(m).Sessions)], [0.5 0.5], 'Color', 'k', 'LineStyle', '--', 'HandleVisibility', 'off')
    ylim([0.45 1])
    ylabel('Fraction correct')
%     title(sprintf('%s', Monkeys(m).Name), 'Interpreter', 'none')
    set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
        'XMinorTick', 'off', 'YMinorTick', 'off',...
        'fontsize',42, 'YGrid', 'on',...
        'fontname', 'Helvetica',...
        'XColor', 'black', 'YColor', 'black')
%     legend
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


