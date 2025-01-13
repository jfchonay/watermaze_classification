clear;
clc;
%%
band = '3_8Hz';

main_dir = ['P:\Jose_Chonay\frequency_sliding\' band];

session = 'probe';
event = 'fsliding';
run = 'start';
ROI = 'FM';

ctrls_stat = load(fullfile(main_dir, ['ctrls_' session '_stat_' event '_' run '.mat']));
ctrls_mobi = load(fullfile(main_dir, ['ctrls_' session '_mobi_' event '_' run '.mat']));
ptnts_stat = load(fullfile(main_dir, ['ptnts_' session '_stat_' event '_' run '.mat']));
ptnts_mobi = load(fullfile(main_dir, ['ptnts_' session '_mobi_' event '_' run '.mat']));

%% cut 200 ms from the start and end

to_cut = 0.2*ctrls_mobi.ctrl_struct.srate;

ctrls_stat_trials = ctrls_stat.ctrl_struct.trials(:,1+to_cut:end-to_cut,:);
ctrls_mobi_trials = ctrls_mobi.ctrl_struct.trials(:,1+to_cut:end-to_cut,:);

ptnts_stat_trials = ptnts_stat.ptnts_struct.trials(:,1+to_cut:end-to_cut,:);
ptnts_mobi_trials = ptnts_mobi.ptnts_struct.trials(:,1+to_cut:end-to_cut,:);

ROI_ctrls_stat = squeeze(mean(ctrls_stat_trials(ctrls_stat.ctrl_struct.ROI.FM,:,:)));
ROI_ctrls_mobi = squeeze(mean(ctrls_mobi_trials(ctrls_mobi.ctrl_struct.ROI.FM,:,:)));

ROI_ptnts_stat = squeeze(mean(ptnts_stat_trials(ptnts_stat.ptnts_struct.ROI.FM,:,:)));
ROI_ptnts_mobi = squeeze(mean(ptnts_mobi_trials(ptnts_mobi.ptnts_struct.ROI.FM,:,:)));

%% permutest between subjects
% compare between subjects (different subject, same conditions)
x_b = ROI_ctrls_stat;
y_b = ROI_ptnts_stat;
% keep a p value threshold of 0.05
[clusters, p_values, t_sums, permutation_distribution ] = permutest(x_b,y_b,0);

%% plot permutest
new_time = ctrls_mobi.ctrl_struct.ori_time{1,1}(1+to_cut:end-to_cut);

figure
plot(new_time, mean(squeeze(ROI_ctrls_stat),2), 'c')
hold on
plot(new_time, mean(squeeze(ROI_ptnts_stat),2), 'b')
for i_c = 1:length(clusters)
    hold on
    fill([new_time(clusters{1,i_c}(1)),new_time(clusters{1,i_c}(end)), ...
        new_time(clusters{1,i_c}(end)),new_time(clusters{1,i_c}(1))], [3.9, 3.9, 4.1, 4.1],'r', 'FaceAlpha',0.3)
end
title(['Average instantenous frequency across all ' run ' trials for ' band ' band in ' ROI])
ylim([3 8])
xlabel('seconds')
ylabel('Hz')
legend('controls static','patients static')