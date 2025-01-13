clear;
clc;
%% 
band = '3_8Hz';
main_dir = dir('P:\Sein_Jeung\Project_Watermaze\WM_EEG_Data\7_epoched');
save_folder = ['P:\Jose_Chonay\frequency_sliding\' band];

if not(isfolder(save_folder))
    mkdir(save_folder);
end
%%
for i_sub = 3:size(main_dir)-1
    sub_name = main_dir(i_sub).name;
    sub_dir = dir(fullfile(main_dir(i_sub).folder, main_dir(i_sub).name));
    addpath(fullfile(main_dir(i_sub).folder, main_dir(i_sub).name));
    session = 'probe';
    task = ["mobi", "stat"];
    event = 'epoched';
    run = ["start","mid","end"];
    low_cut = 3;
    high_cut = 8;
    for i_task = 1:length(task)
        for i_run = 1:length(run)
            to_load = join([sub_name '_' session '_' task(i_task) '_' event '_' run(i_run) '.mat'], '');
            epoched = load(to_load);
            freqsliding = fsliding(epoched, low_cut, high_cut);
            to_save = join([save_folder '\' sub_name '_' session '_' task(i_task) '_fsliding_' run(i_run) '.mat'], '');
            save(to_save, "freqsliding");
        end
    end
end