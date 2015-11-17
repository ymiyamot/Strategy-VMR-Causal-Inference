% This script checks out the Intertrial intervals in the 5 frequency sine
% Visuomotor rotation experiment

%% Load data
load('../../ProcessedData/DAT_5freq.mat')
DAT_5freq = DAT;

%% Extract and plot ITIs
ITIs = cell(1, length(DAT_5freq));
for subj_i = 1:length(DAT_5freq),
    ITIs{subj_i} = cell2mat(cellfun(@(x) diff(x(:, 1), [], 1), ...
        DAT_5freq{subj_i}.trial_start_and_end_times, ...
        'uniformoutput', false)');
end

allITIs = cell2matNaN(ITIs, 2, 'noCell');
allITIsTraining = allITIs(400:end, :);

myFigure([600, 300])
subplot(1, 2, 1); hold on;
plot(conv2(allITIs, 1 / 10 * ones(10, 1)), '.', 'markersize', 0.1)
h = plot(smooth(nanmean(allITIs, 2), 10), 'k', 'linewidth', 3);
xlim([0, 2000])
ylim([0, 10])
xlabel('Trial')
ylabel('ITI (sec)')
title('ITIs over course of expt')
legend(h, 'Mean across subjects')

subplot(1, 2, 2); hold on;
histogram(allITIsTraining(:), 'facecolor', 'w')
xlim([0, 12])
plot(nanmean(allITIsTraining(:)) * ones(1, 2), get(gca, 'ylim'), ...
    'k', 'linewidth', 3)
textinsert(0.4, 0.5, sprintf('Mean ITI = %.2f sec', nanmean(allITIsTraining(:))))
textinsert(0.4, 0.45, sprintf('Stdv ITI = %.2f sec', nanstd(allITIsTraining(:))))
xlabel('ITI (sec)')
ylabel('# Trials')
title('All ITIs (during training) aggregated')
myprintfig('Figures/IntertrialIntvs')