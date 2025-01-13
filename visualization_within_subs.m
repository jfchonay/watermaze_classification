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

pt_s = load([main_dir '\ptnts_probe_stat_fsliding_' run '.mat']);
ct_s = load([main_dir '\ctrls_probe_stat_fsliding_' run '.mat']);

pt_m = load([main_dir '\ptnts_probe_mobi_fsliding_' run '.mat']);
ct_m = load([main_dir '\ctrls_probe_mobi_fsliding_' run '.mat']);

%% cut data

to_cut = 0.2*ct_s.ctrl_struct.srate;
time = ct_s.ctrl_struct.ori_time{1,1}(1+to_cut:end-to_cut);

patients_s = pt_s.ptnts_struct.trials(:,1+to_cut:end-to_cut,:);
controls_s = ct_s.ctrl_struct.trials(:,1+to_cut:end-to_cut,:);

patients_m = pt_m.ptnts_struct.trials(:,1+to_cut:end-to_cut,:);
controls_m = ct_m.ctrl_struct.trials(:,1+to_cut:end-to_cut,:);

control = {controls_s, controls_m};
patient = {patients_s, patients_m};
trials = {control, patient};
colors = {'#7EC724', '#DE7F00', '#F6005C', '#9426C0'};

for c = 1:size(colors,2)
    colors_rgb{c} = hex2rgb(colors{c});
end

%% create the figure

figure;
for pop = 1:2 %control or patient
    for Fi = 1:4 % ROI
        switch pop
            case 1
                subplot(4, 2, Fi*2 - 1);
                
                % find channel indices
                chan_indices = find(cellfun(@(x) any(strcmp(x, config_param.chanGroups(Fi).chan_names)), ct_m.ctrl_struct.chan_labels));
                
                for cond = 1:2 % condition ind 1 = static 2 = mobile
                    % Extract data for the channel of interest (e.g., channel Fi)
                    channel_data = squeeze(mean(trials{1,pop}{1,cond}(chan_indices, :, :),1)); % [1251×10 double]
                    
                    % Compute mean and confidence interval (95%) across
                    % participants for one band
                    mean_data = mean(channel_data, 2); % [1251×1]
                    std_data = std(channel_data, 0, 2); % [1251×1]
                    n_participants = size(channel_data, 2);
                    ci = 1.96 * (std_data / sqrt(n_participants)); % Confidence interval
                    
                    % Plot the mean with confidence interval as a ribbon
                    hold on;
                    fill([time, fliplr(time)], [mean_data + ci; flipud(mean_data - ci)].', colors_rgb{cond}, ...
                        'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Confidence interval ribbon
                    plot(time, mean_data, 'Color', colors_rgb{cond}, 'LineWidth', 1.5); % Mean
                    
                    % Labels and formatting
                    title(['Channel ' config_param.chanGroups(Fi).key ]);
                    legend('','static','','mobile');
                    ylim([9 11]);
                    xlabel('Time (s)');
                    ylabel('Frequency (Hz)');
                end
            case 2
                 subplot(4, 2, Fi*2);
                
                % find channel indices
                chan_indices = find(cellfun(@(x) any(strcmp(x, config_param.chanGroups(Fi).chan_names)), ct_m.ctrl_struct.chan_labels));
                
                for cond = 1:2 % group ind 1 = patients 2 = controls
                    % Extract data for the channel of interest (e.g., channel Fi)
                    channel_data = squeeze(mean(trials{1,pop}{1,cond}(chan_indices, :, :),1)); % [1251×10 double]
                    
                    % Compute mean and confidence interval (95%) across
                    % participants for one band
                    mean_data = mean(channel_data, 2); % [1251×1]
                    std_data = std(channel_data, 0, 2); % [1251×1]
                    n_participants = size(channel_data, 2);
                    ci = 1.96 * (std_data / sqrt(n_participants)); % Confidence interval
                    
                    % Plot the mean with confidence interval as a ribbon
                    hold on;
                    fill([time, fliplr(time)], [mean_data + ci; flipud(mean_data - ci)].', colors_rgb{cond+2}, ...
                        'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Confidence interval ribbon
                    plot(time, mean_data, 'Color', colors_rgb{cond+2}, 'LineWidth', 1.5); % Mean
                    
                    % Labels and formatting
                    title(['Channel ' config_param.chanGroups(Fi).key ]);
                    legend('','static','','mobile');
                    ylim([9 11]);
                    xlabel('Time (s)');
                    ylabel('Frequency (Hz)');
                end

            sgtitle(sprintf(['Instantenous frequency for all ' run ' trials filtered for ' ...
                band ' band \n controls vs patients']));
        end
    end
end
%print(gcf,[main_dir '\within_subs_' run '.png'],'-dpng','-r1500');