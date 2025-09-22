clear;
clc;
%% set up the folder and directory
band = 'alpha';

main_dir = ['P:\Jose_Chonay\frequency_sliding\' band];

subs = 'ptnts';
session = 'probe';
task = 'mobi';
event = 'fsliding';
run = 'end';

%% select files

switch subs
    case 'ptnts'
        files = dir(fullfile(main_dir, ['sub-81*_' session '_' task '_' event '_' run '.mat']));
        labels = repmat({'patients'}, length(files)*24, 1);
        to_save = [subs '_'] ;
    case 'ctrls'
        files_82 = dir(fullfile(main_dir, ['sub-82*_' session '_' task '_' event '_' run '.mat']));
        files_83 = dir(fullfile(main_dir, ['sub-83*_' session '_' task '_' event '_' run '.mat']));
        files_84 = dir(fullfile(main_dir, ['sub-84*_' session '_' task '_' event '_' run '.mat']));
        files = [files_82; files_83; files_84];
        labels = repmat({'controls'}, length(files)*24, 1);
        to_save = [subs '_'] ;
end
%%
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

%% concatenate data and create new structre
for ROI = 1:4
    for i_sub = 1:length(files)
        fsliding_struct = load(fullfile(files(i_sub).folder, files(i_sub).name));
        freqsliding = fsliding_struct.freqsliding;
        chan_indices = find(cellfun(@(x) any(strcmp(x, config_param.chanGroups(ROI).chan_names)), freqsliding.chan_labels));
        to_cut = 0.2*freqsliding.srate; % 200ms in samples
        time = freqsliding.ori_time{1,1}(1+to_cut:end-to_cut);
        % here I have all the trials for one subject
        one_sub = cat(3, freqsliding.fslidemed{:});
        % get only for one ROI
        trials_ROI = squeeze(mean(one_sub(chan_indices, :, :), 1));
        % cut the data
        trials = trials_ROI(1+ to_cut:end-to_cut, :)'; % format (trials X features)
        all_trials(((i_sub-1)*freqsliding.nTrials)+1:freqsliding.nTrials*i_sub,:) = trials;
    end
    [labels{:,2}] = deal(task);
    [labels{:,3}] = deal(config_param.chanGroups(ROI).key);
    trialstrct = struct();
    trialstrct.trials = all_trials;
    trialstrct.labels = labels;
    trialstrct.time = time;
    trialstrct.ROI.key = config_param.chanGroups(ROI).key;
    trialstrct.ROI.name = config_param.chanGroups(ROI).full_name;
    trialstrct.ROI.chan_names = config_param.chanGroups(ROI).chan_names;
    save_name = [to_save session '_' task '_' event '_' config_param.chanGroups(ROI).key '_' run '.mat'];
    save(fullfile(main_dir, save_name), "trialstrct");
end
