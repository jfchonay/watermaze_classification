function bandData = util_WM_bandratios(data, ROI, fBandsIdx)
fBands      = {'theta', 'alpha', 'beta', 'gamma', 'highGamma'};
colNames    = {}; 
for Gi = 1:5
    colNames{end+1} = [ROI, fBands{Gi}];
end

% Separate the data into the fbands
for i_f = 1:5
    switch fBands{i_f}
        case 'theta'
            theta_s = find(fBandsIdx >= 4, 1);
            theta_e = find(fBandsIdx <= 8, 1, 'last');
            theta_trials = mean(data(:, theta_s:theta_e), 2);
        case 'alpha'
            alpha_s = find(fBandsIdx >= 8, 1);
            alpha_e = find(fBandsIdx <= 12, 1,'last');
            alpha_trials = mean(data(:,alpha_s:alpha_e), 2);
        case 'beta'
            beta_s = find(fBandsIdx >= 13, 1);
            beta_e = find(fBandsIdx <= 30, 1, 'last');
            beta_trials = mean(data(:,beta_s:beta_e), 2);
        case 'gamma'
            gamma_s = find(fBandsIdx >= 30, 1);
            gamma_e = find(fBandsIdx <= 80, 1, 'last');
            gamma_trials = mean(data(:, gamma_s:gamma_e), 2);
        case 'highGamma'
            high_s = find(fBandsIdx >= 80, 1);
            high_trials = mean(data(:,high_s:end), 2);
    end
end

dataTable = table(theta_trials, alpha_trials, beta_trials, gamma_trials, high_trials,...
    'VariableNames',colNames);

bandData = dataTable; 

end