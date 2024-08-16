%% clear your workspace
clc; 
clear all; 

%%% OBJECTIVES %%%
% 1. Calculate Delta F/F0
% 2. Calculate peaks for frequency

% [file,path] = uigetfile('*.mat')
% currentFolder = split(path,"/");
% reverseFolder = flip(currentFolder);
% currentFile = reverseFolder(4);
% trialName = currentFile;
% load([path,file]);

% directoryPath = '/Users/kowteckfong/Desktop/MP InVivoCaImaging/IVCI Images & Videos/Raw dff zscore/';
% matFiles = dir(fullfile(directoryPath, '*raw dff zscore.mat'));

cd '/Volumes/DavidSSD/1 Stacked and suite2p data/';
path = '/Volumes/DavidSSD/1 Stacked and suite2p data/18 Oct 23/'; %% CHANGE folder
files = fullfile(path,'18Oct JNJ1.6 Idle','/suite2p/plane0/Fall.mat'); %% CHANGE trial
currentFolder = split(files,"/");
reverseFolder = flip(currentFolder);
currentFile = reverseFolder(4);
trialName = currentFile;
load(files);

%% collecting raw fluorescence, deconvoluted, and coordinates
rawF = [];
deconvSpike = [];
x_loc = [];
y_loc = [];
for a = 1:size(iscell,1)
    if iscell(a,1) == 1
        if F(a,:) == 0
            continue
        else
        rawF = [rawF;F(a,:)];
        deconvSpike = [deconvSpike; spks(a,:)];
        x_loc = [x_loc, stat{a}.med(2)];
        y_loc = [y_loc, stat{a}.med(1)];
        end
    end
end

%% delta F/F0
rawF_low20 = prctile(rawF,20);
for t = 1:size(rawF,2)
    window_size = ceil(3*14.6);
    window_start = max(1, t - window_size/2);
    window_end = min(size(rawF,2), t + window_size/2);
    rawF_baseline(t) = min(rawF_low20(window_start:window_end));
end
% rawF_baseline = movmean(rawF_baseline,50);
rawF_baseline = smoothdata(rawF_baseline,'movmean',window_size,'omitnan');

for i = 1:size(rawF,1)
ZRawFluo(i,:) = (rawF(i,:) - rawF_baseline) ./ rawF_baseline;
end

% figure();plot(ZRawFluo)
% figure();plot(mean(ZRawFluo,1))

%% z-score
for i = 1:size(rawF,1)
    mSpike = mean(rawF(i,:),'all');
    stdSpike = std(rawF(i,:),0,'all');
    ZDeconvFluo(i,:) = (rawF(i,:) - mSpike) ./ stdSpike;
end

%% find peak
for j = 1:size(ZRawFluo,1) % frequency
    % establish baseline
    mFluo = mean(ZRawFluo(j,:),'all');
    mStd = std(ZRawFluo(j,:),0,'all');
    baselinePeak = mFluo + 2*mStd;
    % delta f over f
    eachNeu_dff = ZRawFluo(j,:);
    Freq_dff(j,:) = eachNeu_dff > baselinePeak;
    belowB = find(eachNeu_dff <= baselinePeak);
    eachNeu_dff(belowB) = NaN;
    Amp_dff(j,:) = eachNeu_dff;

    % z score
    eachNeu_z = ZDeconvFluo(j,:);
    Freq_z_pos(j,:) = eachNeu_z > 2;
    Freq_z_neg(j,:) = eachNeu_z < -2;
end
Freq = {Freq_dff,Freq_z_pos,Freq_z_neg};

%% saving data
locs = {x_loc,y_loc};
allData = {
    'rawF', rawF;    
    'deconv fluo', deconvSpike;
    'delta F over F',ZRawFluo;
    'z score raw f',ZDeconvFluo;
    'freq: dff, z pos, z neg',Freq;
    'amp dff', Amp_dff;
    'locations',locs};

saveLocation = '/Users/kowteckfong/Desktop/MP InVivoCaImaging/IVCI Images & Videos/Raw dff zscore/';
save([saveLocation,trialName{1},' raw dff zscore.mat'],'allData');