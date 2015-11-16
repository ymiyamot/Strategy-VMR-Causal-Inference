%%%% This script examines broad learning in the sum of sines experiment

%% Load Data
load('../ProcessedData/DAT_5freq.mat')
DAT_5freq = DAT;

load('../ProcessedData/DAT_7freq.mat')
DAT_7freq = DAT;

% Add function path
addpath(genpath('HelperFunctions'));