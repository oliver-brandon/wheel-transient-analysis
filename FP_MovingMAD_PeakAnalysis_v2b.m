clear all; clc; close all;
%removes polyfit warning%
warning('off','all');
VERSION = "2b";
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ________   __       _______ .______   .__   __.  __   __  ___     __          ___      .______  %
% |       /  |  |     |   ____||   _  \  |  \ |  | |  | |  |/  /    |  |        /   \     |   _  \ % 
% `---/  /   |  |     |  |__   |  |_)  | |   \|  | |  | |  '  /     |  |       /  ^  \    |  |_)  |% 
%    /  /    |  |     |   __|  |   _  <  |  . `  | |  | |    <      |  |      /  /_\  \   |   _  < % 
%   /  /----.|  `----.|  |____ |  |_)  | |  |\   | |  | |  .  \     |  `----./  _____  \  |  |_)  |% 
%  /________||_______||_______||______/  |__| \__| |__| |__|\__\    |_______/__/     \__\ |______/ % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%FP_MovingMAD_PeakAnalysis.m created by Brandon L. Oliver, M.A., adapted
%from Barker et al. (2017)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% input paramaters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 20; % time threshold below which we will discard
N = 100; %Downsample N times
sampleRate = 1017/N;
session_duration = 3600;
window_size_seconds = 15;
%Stream Stores%
DLS_ISOS = 'x405A'; % name of the 405A store
DLS_DA = 'x465A'; % name of the 465A store
NAc_ISOS = 'x405C'; % name of the 405C store
NAc_DA = 'x465C'; % name of the 465C store

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf("Version: %s\n",VERSION)
myDir = uigetdir('/Users/brandon/personal-drive/wheel-peak-test/tanks');
BLOCKPATH = myDir;
data = TDTbin2mat(BLOCKPATH, 'T2', session_duration, 'TYPE', {'streams'});



%time array used for all streams%
time = (1:length(data.streams.(DLS_DA).data))/data.streams.(DLS_DA).fs;
%removes the first (t) seconds where the data is wild due to turning on LEDs%

ind = find(time>t,1);% find first index of when time crosses threshold
time = time(ind:end); % reformat vector to only include allowed time
data.streams.(DLS_DA).data = data.streams.(DLS_DA).data(ind:end);
data.streams.(DLS_ISOS).data = data.streams.(DLS_ISOS).data(ind:end);
data.streams.(NAc_DA).data = data.streams.(NAc_DA).data(ind:end);
data.streams.(NAc_ISOS).data = data.streams.(NAc_ISOS).data(ind:end);

%downsample streams and time array by N times%
data.streams.(DLS_ISOS).data = downsample(data.streams.(DLS_ISOS).data, N);
data.streams.(DLS_DA).data = downsample(data.streams.(DLS_DA).data, N);
data.streams.(NAc_ISOS).data = downsample(data.streams.(NAc_ISOS).data, N);
data.streams.(NAc_DA).data = downsample(data.streams.(NAc_DA).data, N);
time = downsample(time, N);

%detrend & dFF%
%465A%
bls = polyfit(data.streams.(DLS_ISOS).data,data.streams.(DLS_DA).data,1);
Y_fit_all = bls(1) .* data.streams.(DLS_ISOS).data + bls(2);
Y_dF_all = data.streams.(DLS_DA).data - Y_fit_all; %dF (units mV) is not dFF
dFF = 100*(Y_dF_all)./Y_fit_all;
std_dFF = std(double(dFF));
detrend_465A = detrend(dFF);
detrend_465A = zscore(detrend_465A);
%465C%
bls2 = polyfit(data.streams.(NAc_ISOS).data,data.streams.(NAc_DA).data,1);
Y_fit_all2 = bls2(1) .* data.streams.(NAc_ISOS).data + bls2(2);
Y_dF_all2 = data.streams.(NAc_DA).data - Y_fit_all2; %dF (units mV) is not dFF
dFF2 = 100*(Y_dF_all2)./Y_fit_all2;
std_dFF2 = std(double(dFF2));
detrend_465C = detrend(dFF2);
detrend_465C = zscore(detrend_465C);

%set up moving MAD window%

peak_indicies_DLS = [];
peak_indicies_NAc = [];
signal_chunksA = [];
signal_chunksC = [];

% determines if window_size evenly divides into the streams and pads with 
% NaN if not
window_size = ceil(window_size_seconds * sampleRate);
remainderA = ceil((window_size - mod(length(detrend_465A),window_size)));
remainderC = ceil((window_size - mod(length(detrend_465C),window_size)));
if remainderA == window_size
    signalA = detrend_465A;
    signalA_ind = 1:length(signalA);
else
    padding = NaN(1,remainderA);
    signalA = [detrend_465A,padding];
    signalA_ind = 1:length(signalA);
end

if remainderC == window_size
    signalC = detrend_465C;
    signalC_ind = 1:length(signalC);
else
    padding = NaN(1,remainderC);
    signalC = [detrend_465C,padding];
    signalC_ind = 1:length(signalC);
end
time2 = (1:length(signalA))/sampleRate;

% fills signal_chunksA and signal_chunksC with the actual signal chunks
for ii = 1:window_size:length(signalA)
    x1 = signalA(1,ii:ii+window_size-1);
    signal_chunksA = [signal_chunksA; x1];
end
for ii = 1:window_size:length(signalC)
    x1 = signalC(1,ii:ii+window_size-1);
    signal_chunksC = [signal_chunksC; x1];
end

% finds peaks using previous window
signalA_pks = [];
signalA_pks(1,size(signal_chunksA,2)) = nan;
for jj = 2:height(signal_chunksA)
    signalA_median = median(signal_chunksA(jj-1,:));
    signalA_MAD = mad(signal_chunksA(jj-1,:),1);
    signalA_thresh = (signalA_median + (3 * signalA_MAD));
    pks = signal_chunksA(jj,:) > signalA_thresh;

    signalA_pks(jj,:) = pks;
end

% finds peaks using previous window
signalC_pks = [];
signalC_pks(1,size(signal_chunksC,2)) = nan;
for jj = 2:height(signal_chunksC)
    signalC_median = median(signal_chunksC(jj-1,:));
    signalC_MAD = mad(signal_chunksC(jj-1,:),1);
    signalC_thresh = (signalA_median + (3 * signalC_MAD));
    pks2 = signal_chunksC(jj,:) > signalC_thresh;

    signalC_pks(jj,:) = pks2;
end

% reshapes signal peak matrix into size of signal
signalA_pks = [signalA_pks(1, :) reshape(signalA_pks(2:end, :)', 1, [])];
signalC_pks = [signalC_pks(1, :) reshape(signalC_pks(2:end, :)', 1, [])];

% calculates number of peaks and frequency
totPks_A = length(signalA_pks(signalA_pks == 1));
transients_A = totPks_A/session_duration;
totPks_C = length(signalC_pks(signalC_pks == 1));
transients_C = totPks_C/session_duration;

% Plots the signal with markers for the peaks
plot(time2, signalA);
hold on;

plot(time2(signalA_pks == 1), signalA(signalA_pks == 1), 'ro', 'MarkerSize', 10, 'LineWidth', 2);

hold off;



