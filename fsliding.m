function [fsstruct] = fsliding(epochstruct, low_end, high_end)
% define the basic structure parameters of the new data
ftEEG = epochstruct.ftEEG;
fsEEG = struct();

fsEEG.srate    = ftEEG.hdr.Fs;
fsEEG.stim_dur = size(ftEEG.trial{1,1}, 2);
fsEEG.chan_labels = ftEEG.label;
fsEEG.ori_time = ftEEG.time;

fsEEG.nTrials = ftEEG.hdr.nTrials;
fsEEG.ROI.FM = [33, 34, 35, 57, 64];
fsEEG.ROI.PM = [73, 74, 75, 91, 96];
fsEEG.ROI.LT = [1, 48, 79, 77];
fsEEG.ROI.RT = [24, 52, 82, 84];
% define the parameters of the bandpass filter
f_order = 125; 

f_Nq = fsEEG.srate / 2;
f_low = low_end;
f_high = high_end;
t_bands = 0.15;
f_bands = [0, (1-t_bands)*f_low, f_low, f_high, f_high*(1+t_bands), f_Nq] / f_Nq;

f_amps = [0, 0, 1, 1, 0, 0];

f_coef = firls(f_order, f_bands, f_amps);

fsEEG.bandfilter.order = f_order;
fsEEG.bandfilter.low_end = low_end;
fsEEG.bandfilter.high_end = high_end;
fsEEG.bandfilter.transition_bandprctg = t_bands * 100;

% define the median filter paramters
mf_order = 10;
filtord = round(linspace(10,400,mf_order))/2;
filtord = round( filtord/(1000/fsEEG.srate) );
timeoi  = (fsEEG.ori_time{1,1}(1)*1000:(1000/fsEEG.srate):fsEEG.ori_time{1,1}(end)*1000)/1000;
timeoiidx = dsearchn(fsEEG.ori_time{1, 1}',timeoi');

fsEEG.medianfilter.order = mf_order;

% run the calculations of frequency sliding

fsEEG.filtered = cell(1, fsEEG.nTrials);
fsEEG.hilbert = cell(1, fsEEG.nTrials);
fsEEG.phased = cell(1, fsEEG.nTrials);
fsEEG.gam = cell(1, fsEEG.nTrials);
fslidemed = zeros(numel(fsEEG.chan_labels), fsEEG.stim_dur,mf_order);

for n_t = 1:fsEEG.nTrials
    fsEEG.filtered{1,n_t} =  bandpass_firls(f_coef, ftEEG.trial{1,n_t});
    fsEEG.hilbert{1, n_t} = transpose(hilbert(transpose(fsEEG.filtered{1,n_t})));
    fsEEG.phased{1, n_t} = angle(fsEEG.hilbert{1, n_t});
    % change in the order of the formula
    fsEEG.gam{1,n_t} = diff(fsEEG.srate*unwrap(fsEEG.phased{1,n_t},[],2),1,2)/(2*pi);
    % start of the median filter application
    frslide = zeros(numel(fsEEG.chan_labels), numel(fsEEG.gam{1,n_t}(1,:)));
    frslide(:,:) = fsEEG.gam{1,n_t}; 
    for oi=1:mf_order
        for ti=1:length(timeoi)
            temp = sort(frslide(:,max(timeoiidx(ti)-filtord(oi),1):min(timeoiidx(ti)+filtord(oi),numel(fsEEG.ori_time{1,n_t})-1)),2);
            fslidemed(:,ti,oi) = temp(:,floor(size(temp,2)/2)+1);
        end
    end
    fsEEG.fslidemed{1,n_t} = squeeze(median(fslidemed, 3));
end

fsstruct = fsEEG;
end