%%%% This script examines broad learning in the sum of sines experiment

%% Load Data
load('../../ProcessedData/DAT_5freq.mat')
DAT_5freq = DAT;

load('../../ProcessedData/DAT_7freq.mat')
DAT_7freq = DAT;

% Add function path
addpath(genpath('../HelperFunctions'));

%% Extract basic learning curves
% Find first trial of training
rotSeq = DAT_5freq{1}.Cursor_dir - DAT_5freq{1}.Head_dir;
startInd = find(abs(rotSeq) ~= 0 ...
    & abs(rotSeq) < 20 & abs(rotSeq) > 1e-3, 1) - 2;
% startInd = 400; % Hard-coded

[n5freq, trainInd_5freq, ...
    aiming_5freq, handDir_5freq, impDir_5freq, rotDir_5freq, err_5freq]...
    = extractLearningCurves(DAT_5freq, tgt_all_N_5freq, startInd);
exptLen = length(rotDir_5freq) - 1;

[n7freq, trainInd_7freq, ...
    aiming_7freq, handDir_7freq, impDir_7freq, rotDir_7freq, err_7freq]...
    = extractLearningCurves(DAT_7freq, tgt_all_N_7freq, startInd);

%% Outlier removal
for subj_i = 1:length(DAT_5freq),
    
    figure; hold on;
    plot(handDir_5freq, 'b.')
    plot(outlierRemove(handDir_5freq, 2, 3), 'r.')
    
    bsxfun(@minus, rotDir_5freq, outlierRemove(handDir_5freq, 2, 3))
    
    figure; hold on;
    plot(aiming_5freq, 'b.')
    plot(outlierRemove(aiming_5freq, 2, 3), 'r.')
    
    figure; plot(handDir_5freq)
    omitInd = find(abs(DAT_5freq{subj_i}.thetaT - DAT_5freq{subj_i}.thetaD) > 15);
    DAT_5freq{subj_i}.thetaT(omitInd) = NaN;
    DAT_5freq{subj_i}.thetaD(omitInd) = NaN;
end

%% Find driven frequencies and corresponding periods
tmpFFT = nanfft(rotDir_5freq(1:end - 1, 1));
fft_pert_5freq = tmpFFT(1:length(tmpFFT) / 2 + 1) / exptLen;

pertFreq_5freq = find(abs(fft_pert_5freq) > 1e-3);
% pertFreq_5freq = 1 + 2.^[0:4]; % Hard-coded version
pertPeriod_5freq = (size(handDir_5freq, 1) - 1) ./ (pertFreq_5freq - 1);

%% Perform FFT transform of learning curves

[fft_hand_5freq, fft_aiming_5freq, fft_imp_5freq, fft_err_5freq]...
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

%%

rectifyAngleVec

phaseDiff = NaN(size(fft_aiming_5freq));
for subj_i = 1:size(fft_aiming_5freq, 2)
    expPhase = 180 / pi * phase(fft_aiming_5freq(:, subj_i));
    impPhase = 180 / pi * phase(fft_imp_5freq(:, subj_i));
    phaseDiff(:, subj_i) = rectifyAngle(expPhase - impPhase);
end



figure; hold on;
plot((1:size(phaseDiff, 1)) - 1, circavg(phaseDiff, 2), '.')


figure; hold on;
plot([(size(handDir_5freq, 1) - 1) ./ (1:size(phaseDiff, 1) - 1), 0], circavg(phaseDiff, 2), '.')

plot(180 / pi * phase(fft_imp_5freq(:, 1)))



%% Calculate overall errors

err5freqOutrmv...
    = bsxfun(@minus, rotDir_5freq, outlierRemove(handDir_5freq, 2, 3));
perfMSE2 = nanmean(err5freqOutrmv.^2, 1);
perfMSE = nanmean(err_5freq.^2, 1);

perfMSE7freq = nanmean(err_7freq.^2, 1);

% aimOverall = nanmean(aiming_5freq.^2, 1);
aimOverall = nanmean(diff(aiming_5freq, [], 1).^2, 1);
% impOverall = nanmean(impDir_5freq.^2, 1);
impOverall = nanmean(diff(impDir_5freq, [], 1).^2, 1);

% aimOverall7freq = nanmean(aiming_7freq.^2, 1);
aimOverall7freq = nanmean(diff(aiming_7freq, [], 1).^2, 1);
% impOverall7freq = nanmean(impDir_7freq.^2, 1);
impOverall7freq = nanmean(diff(impDir_7freq, [], 1).^2, 1);

figure; 
subplot(2, 2, 1); hold on;
myRegress(perfMSE, aimOverall)
subplot(2, 2, 2); 
myRegress(perfMSE2, aimOverall)
subplot(2, 2, 3); hold on;
myRegress(perfMSE, impOverall)
subplot(2, 2, 4); 
myRegress(perfMSE2, impOverall)

figure; 
subplot(2, 2, 1); hold on;
myRegress(perfMSE7freq, aimOverall7freq)
subplot(2, 2, 2); 
myRegress(perfMSE7freq, impOverall7freq)


[bAim, ~, rAim, ~, statsAim] = regress(perfMSE', [ones(size(aimOverall')), aimOverall']);
[bImp, ~, rImp, ~, statsImp] = regress(perfMSE', [ones(size(impOverall')), impOverall']);
[bAimImp, ~, rAimImp, ~, statsAimImp] ...
    = regress(perfMSE', [ones(size(aimOverall')), aimOverall', impOverall']);

statsAim
statsImp
statsAimImp

statsAimImp(1) - statsAim(1)
statsAimImp(1) - statsImp(1)
help parR
1 - statsImp(1)
[p, partR, F] = parR(rAim, 2, rAimImp, 3)
[p, partR, F] = parR(rImp, 2, rAimImp, 3)


figure; 
subplot(1, 2, 1); hold on;
myRegress(perfMSE, aimOverall)

subplot(1, 2, 2); hold on;
myRegress(perfMSE2, aimOverall)
subplot(1, 2, 2); 
myRegress(perfMSE, aimOverall)
figure; plot(aimOverall, aimOverall2, 'o')
figure; hold on;
plot(aiming_5freq)
plot(nanmean(aiming_5freq, 2), 'linewidth', 3)


figure; imagesc(diff(aiming_5freq, [], 1))



figure; histogram()


figure; plot(err_5freq)

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