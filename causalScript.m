%%%% This script examines broad learning in the sum of sines experiment

%% Load Data
load('../../ProcessedData/DAT_5freq.mat')
DAT_5freq = DAT;

% load('../ProcessedData/DAT_7freq.mat')
% DAT_7freq = DAT;

% Add function path
addpath(genpath('../HelperFunctions'));

%% Extract basic learning curves
startInd = 400; % First trial of training
[n5freq, trainInd_5freq, ...
    aiming_5freq, handDir_5freq, impDir_5freq, rotDir_5freq, err_5freq]...
    = extractLearningCurves(DAT_5freq, tgt_all_N_5freq, startInd);
exptLen = length(rotDir_5freq) - 1;

rotSeq = DAT_5freq{1}.Cursor_dir - DAT_5freq{1}.Head_dir;
find(abs(rotSeq) ~= 0 & abs(rotSeq) < 20 & abs(rotSeq) > 1e-3, 1)

figure; plot()

%% Find driven frequencies and corresponding periods
tmpFFT = nanfft(rotDir_5freq(1:end - 1, 1));
fft_pert_5freq = tmpFFT(1:length(tmpFFT) / 2 + 1) / exptLen;

pertFreq_5freq = find(abs(fft_pert_5freq) > 1e-3);
% pertFreq_5freq = 1 + 2.^[0:4]; % Hard-coded version
pertPeriod_5freq = (size(handDir_5freq, 1) - 1) ./ (pertFreq_5freq - 1);

%%
figure; hold on;
plot(handDir_5freq)
plot(nanmean(handDir_5freq, 2), 'linewidth', 3, 'color', 'k')



    
[fft_hand_5freq, fft_aiming_5freq, ...
    fft_imp_5freq, fft_err_5freq, fft_pert_5freq]...
    = deal(NaN((size(handDir_5freq, 1) - 1) / 2 + 1, size(handDir_5freq, 2)));
for subj_i = 1:size(handDir_5freq, 2),
    tmpFFT = nanfft(aiming_5freq(1:end - 1, subj_i));
    fft_aiming_5freq(:, subj_i) = tmpFFT(1:length(tmpFFT) / 2 + 1) / exptLen;
    
    tmpFFT = nanfft(handDir_5freq(1:end - 1, subj_i));
    fft_hand_5freq(:, subj_i) = tmpFFT(1:length(tmpFFT) / 2 + 1) / exptLen;
    
    tmpFFT = nanfft(impDir_5freq(1:end - 1, subj_i));
    fft_imp_5freq(:, subj_i) = tmpFFT(1:length(tmpFFT) / 2 + 1) / exptLen;
    
    tmpFFT = nanfft(err_5freq(1:end - 1, subj_i));
    fft_err_5freq(:, subj_i) = tmpFFT(1:length(tmpFFT) / 2 + 1) / exptLen;
    

end
tmpFFT = nanfft(rotDir_5freq(1:end - 1));
fft_pert_5freq = tmpFFT(1:length(tmpFFT) / 2 + 1) / exptLen;