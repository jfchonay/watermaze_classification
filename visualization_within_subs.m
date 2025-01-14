% config structure used for selection ROI
config_param.chanGroups(1).key           = 'FM';
config_param.chanGroups(1).full_name     = 'Frontal-midline';
config_param.chanGroups(1).chan_names    = {'y1','y2','y3','y25','y32'}; 
config_param.chanGroups(2).key           = 'PM';
config_param.chanGroups(2).full_name     = 'Parietal-midline';
config_param.chanGroups(2).chan_names    = {'r9', 'r10', 'r11', 'r27', 'r32'}; 
config_param.chanGroups(3).key           = 'LT';
config_param.chanGroups(3).full_name     = 'Left-temporal';
config_param.chanGroups(3).chan_names    = {'g1', 'y16', 'r15', 'r13'}; 
config_param.chanGroups(4).key           = 'RT';
config_param.chanGroups(4).full_name     = 'Right-temporal';
config_param.chanGroups(4).chan_names    = {'g24','y20', 'r18', 'r20'}; 

%% Load results structs
band = 'alpha';
run = 'end';
main_dir = ['P:\Jose_Chonay\frequency_sliding\' band];
% load populations and task
pt_s = load([main_dir '\ptnts_probe_stat_fsliding_' run '.mat']);
ct_s = load([main_dir '\ctrls_probe_stat_fsliding_' run '.mat']);

pt_m = load([main_dir '\ptnts_probe_mobi_fsliding_' run '.mat']);
ct_m = load([main_dir '\ctrls_probe_mobi_fsliding_' run '.mat']);

%% cut data
% cut the first and the las 200ms to avoid the noise created by the
% filtering windows also cut the time vector
to_cut = 0.2*ct_s.ctrl_struct.srate;
time = ct_s.ctrl_struct.ori_time{1,1}(1+to_cut:end-to_cut);

patients_s = pt_s.ptnts_struct.trials(:,1+to_cut:end-to_cut,:);
controls_s = ct_s.ctrl_struct.trials(:,1+to_cut:end-to_cut,:);

patients_m = pt_m.ptnts_struct.trials(:,1+to_cut:end-to_cut,:);
controls_m = ct_m.ctrl_struct.trials(:,1+to_cut:end-to_cut,:);
% cells of the data nested in condition then by population
control = {controls_s, controls_m};
patient = {patients_s, patients_m};
trials = {control, patient};
% palette of colors to be used in HEX
colors = {'#7ED857','#B157D8'};
% transform colors to rgb by using a function to convert to binary
for c = 1:size(colors,2)
    colors_rgb{c} = hex2rgb(colors{c});
end

%% create the figure

figure;
for pop = 1:2 %control or patient
    for Fi = 1:4 % ROI
        % find channel indices
        chan_indices = find(cellfun(@(x) any(strcmp(x, config_param.chanGroups(Fi).chan_names)), ct_m.ctrl_struct.chan_labels));
            
        switch pop
            case 1 % controls
                subplot(4, 2, Fi*2 - 1);
                % save legends as a cell so we can add empty spaces later
                legends_2b = {'';'stationary';'';'mobile'};
               
                for cond = 1:2 % condition ind 1 = stationary 2 = mobile
                    % Extract data for the channel of interest (e.g., channel Fi)
                    channel_data = squeeze(mean(trials{1,pop}{1,cond}(chan_indices, :, :),1)); % [1151x20 double]
                    
                    % Compute mean and confidence interval (95%) across
                    % participants for one band
                    mean_data = mean(channel_data, 2); % [1151×1]
                    std_data = std(channel_data, 0, 2); % [1151×1]
                    n_participants = size(channel_data, 2);
                    ci = 1.96 * (std_data / sqrt(n_participants)); % Confidence interval
                    
                    % Plot the mean with confidence interval as a ribbon
                    hold on;
                switch cond
                    case 1
                        fill([time, fliplr(time)], [mean_data + ci; flipud(mean_data - ci)].', colors_rgb{pop}, ...
                        'FaceAlpha', 0.4, 'EdgeColor', 'none'); % Confidence interval ribbon
                        plot(time, mean_data, '--', 'Color', colors_rgb{pop}, 'LineWidth', 1.5); % Mean
                    case 2
                        fill([time, fliplr(time)], [mean_data + ci; flipud(mean_data - ci)].', colors_rgb{pop}, ...
                        'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Confidence interval ribbon
                        plot(time, mean_data, '-', 'Color', colors_rgb{pop}, 'LineWidth', 1.5); % Mean
                end
                    % Labels and formatting
                    title(['Channel ' config_param.chanGroups(Fi).key ]);
                    ylim([9 11]);
                    xlabel('Time (s)');
                    ylabel('Frequency (Hz)');
                end
                % run a permutation between subjects
                x = squeeze(mean(trials{1,pop}{1,1}(chan_indices, :, :),1)); % stationary
                y = squeeze(mean(trials{1,pop}{1,2}(chan_indices, :, :),1)); % mobile
                % paramaters set to within subject analysis and p
                % threshold of 0.05
                [clusters, p_values, t_sums] = permutest(y,x,1,0.05,[],1);
                % when t_sums is empty there were no significant clusters
                if ~isempty(t_sums)
                    % for every significan cluster we will mark the
                    % position in the time axis
                    for i_c = 1:length(clusters)
                        hold on
                        y_value = 9.3; % Set the y-coordinate
                        plot(time(clusters{1,i_c}), y_value * ones(size(time(clusters{1,i_c}))),...
                            'k.', 'LineWidth', 0.1);
                        % add one empty space for every cluster line
                        % created
                        legends_2b{end+1,1} = '';
                    end
                end
                % use the legends created
                legend(legends_2b);
            case 2 % patients
                subplot(4, 2, Fi*2);
                % save legends as a cell so we can add empty spaces later
                legends_2b = {'';'stationary';'';'mobile'};
                
                for cond = 1:2 % group ind 1 = stationary 2 = mobile
                    % Extract data for the channel of interest (e.g., channel Fi)
                    channel_data = squeeze(mean(trials{1,pop}{1,cond}(chan_indices, :, :),1)); % [1151×10 double]
                    
                    % Compute mean and confidence interval (95%) across
                    % participants for one band
                    mean_data = mean(channel_data, 2); % [1151×1]
                    std_data = std(channel_data, 0, 2); % [1151×1]
                    n_participants = size(channel_data, 2);
                    ci = 1.96 * (std_data / sqrt(n_participants)); % Confidence interval
                    
                    % Plot the mean with confidence interval as a ribbon
                    hold on;
                switch cond
                    case 1
                        fill([time, fliplr(time)], [mean_data + ci; flipud(mean_data - ci)].', colors_rgb{pop}, ...
                        'FaceAlpha', 0.4, 'EdgeColor', 'none'); % Confidence interval ribbon
                        plot(time, mean_data, '--', 'Color', colors_rgb{pop}, 'LineWidth', 1.5); % Mean
                    case 2
                        fill([time, fliplr(time)], [mean_data + ci; flipud(mean_data - ci)].', colors_rgb{pop}, ...
                        'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Confidence interval ribbon
                        plot(time, mean_data, '-', 'Color', colors_rgb{pop}, 'LineWidth', 1.5); % Mean
                end                    
                    % Labels and formatting
                    title(['Channel ' config_param.chanGroups(Fi).key ]);
                    legend('','stationary','','mobile');
                    ylim([9 11]);
                    xlabel('Time (s)');
                    ylabel('Frequency (Hz)');
                end
                % run a permutation between subjects
                x = squeeze(mean(trials{1,pop}{1,1}(chan_indices, :, :),1)); % stationary
                y = squeeze(mean(trials{1,pop}{1,2}(chan_indices, :, :),1)); % mobile
                % paramaters set to within subject analysis and p
                % threshold of 0.05
                [clusters, p_values, t_sums] = permutest(y,x,1,0.05,[],1);
                % when t_sums is empty there were no significant clusters
                if ~isempty(t_sums)
                    % for every significan cluster we will mark the
                    % position in the time axis
                    for i_c = 1:length(clusters)
                        hold on
                        y_value = 9.3; % Set the y-coordinate
                        plot(time(clusters{1,i_c}), y_value * ones(size(time(clusters{1,i_c}))),...
                            'k.', 'LineWidth', 0.1);
                        % add one empty space for every cluster line
                        % created
                        legends_2b{end+1,1} = '';
                    end
                end
            legend(legends_2b);
            sgtitle(sprintf(['Instantenous frequency for all ' run ' trials filtered for ' ...
                band ' band \n controls vs patients']));
        end
    end
end
%print(gcf,[main_dir '\within_subs_' run '.png'],'-dpng','-r1500');