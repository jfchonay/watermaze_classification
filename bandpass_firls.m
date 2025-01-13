function [filtered_data] = bandpass_firls(coefficients, data)
filtered_data = zeros(size(data));
    for n_c = 1:size(data,1)
        filtered_data(n_c,:) = filtfilt(coefficients,1,data(n_c,:));
    end
end