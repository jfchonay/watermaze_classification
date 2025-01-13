clear;
clc;
%% set up the folder and directory
band = '3_8Hz';

main_dir = ['P:\Jose_Chonay\frequency_sliding\' band];

subs = 'ctrls';
session = 'probe';
task = 'mobi';
event = 'fsliding';
run = 'end';

%% select files

switch subs
    case 'ptnts'
        files = dir(fullfile(main_dir, ['sub-81*_' session '_' task '_' event '_' run '.mat']));
        all_sub_trials = zeros([129 1251 length(files)]);
    case 'ctrls'
        files_82 = dir(fullfile(main_dir, ['sub-82*_' session '_' task '_' event '_' run '.mat']));
        files_83 = dir(fullfile(main_dir, ['sub-83*_' session '_' task '_' event '_' run '.mat']));
        files_84 = dir(fullfile(main_dir, ['sub-84*_' session '_' task '_' event '_' run '.mat']));
        files = [files_82; files_83; files_84];
        all_sub_trials = zeros([129 1251 length(files)]);
end

%% concatenate data and create new structre
for i_sub = 1:length(files)
    fsliding_struct = load(fullfile(files(i_sub).folder, files(i_sub).name));
    freqsliding = fsliding_struct.freqsliding;
    % here I have all the trials for one subject
    all_trials = cat(3, freqsliding.fslidemed{:});
    avrg_trials = mean(all_trials, 3);
    all_sub_trials(:,:,i_sub) = avrg_trials;
end

%% save new file

switch subs
    case 'ptnts'
        ptnts_struct = rmfield(freqsliding, {'filtered','hilbert','phased','gam','fslidemed'});
        ptnts_struct.trials = all_sub_trials;
        to_save = [subs '_' session '_' task '_' event '_' run '.mat'];
        save(fullfile(main_dir, to_save),"ptnts_struct");
    case 'ctrls'
        ctrl_struct = rmfield(freqsliding, {'filtered','hilbert','phased','gam','fslidemed'});
        ctrl_struct.trials = all_sub_trials; 
        to_save = [subs '_' session '_' task '_' event '_' run '.mat'];   
        save(fullfile(main_dir, to_save),"ctrl_struct");
end

