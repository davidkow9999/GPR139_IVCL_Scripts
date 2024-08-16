clear;clc
directoryPath = '/Users/kowteckfong/Desktop/MP InVivoCaImaging/IVCI Images & Videos/Raw dff zscore/';
cd '/Users/kowteckfong/Desktop/MP InVivoCaImaging/IVCI Images & Videos/Raw dff zscore/'
matFiles = dir(fullfile(directoryPath, '*raw dff zscore.mat'));

%%% OBJECTIVES %%%
% 1. Further quantify frequency parameters
% 2. Calculate probability distribution
% 3. Statistical tests
% 4. Generating graphs

DMSO = [1,16,24,28,30,49];
JNJ = [3,26,32,39,51];
MK_JNJ = [5,13,18,34,41,53,61];
MK_JNJ_MK = [7,15,20,36,43,55,63];
TAK = [11,37,47,59,67];
MK_TAK = [8,21,44,56,64];
MK_TAK_MK = [10,23,46,58,66];
MK = [MK_JNJ_MK,MK_TAK_MK];

groups = {DMSO,JNJ,TAK,MK,MK_JNJ_MK,MK_JNJ,MK_TAK_MK,MK_TAK}; % 1v4, 1v5 1v7
% TAK [0.6 0.4 1]; DMSO [0 0 1]; JNJ [0 1 0.8]; MK [1 0.5 0.2]
visibility = 'on';
save_mode = 0; % 1 to save, 0 to not save

%% Amplitude and frequency
% COMMENT: Normalize by DMSO?
for k = 1:size(groups,2) % all groups
    GroupSumF = [];
    GroupSumA = [];
    dff_stack = [];
    dff_allsamples = [];
    for j = 1:size(groups{k},2) % all samples in a group
        load(matFiles(groups{k}(j)).name);
        dff_frame = allData{3,2}; % delta f over f0
        dffFreq_frame = allData{5,2}{1,1}; % z raw frequency
        dffFreq_min = dffFreq_frame;
%         dffFreq_min = reshape(dffFreq_frame(i,1:8100),15,[]);
        dff_allsamples = [dff_allsamples;dffFreq_min(:,1:8100)];
        zposFreq_frame = allData{5,2}{1,2};
        znegFreq_frame = allData{5,2}{1,3};
        dffAmp = allData{6,2}; % z raw amplitude
        for i = 1:size(dffFreq_frame,1)
            DFF_min = mean(reshape(dff_frame(i,1:8100),15,[]),1);
            Freq(i,:) = mean(reshape(dffFreq_frame(i,1:8100),15,[]),1,'omitnan'); % mean frequency by min
            zposFreq(i,:) = mean(reshape(zposFreq_frame(i,1:8100),15,[]),1,'omitnan');
            znegFreq(i,:) = mean(reshape(znegFreq_frame(i,1:8100),15,[]),1,'omitnan');
            Amp(i,:) = mean(reshape(dffAmp(i,1:8100),15,[]),1,'omitnan');
        end
        dff_stack = [dff_stack; DFF_min];
        Freq(Freq == 0) = NaN;
        zposFreq(zposFreq == 0) = NaN;
        znegFreq(znegFreq == 0) = NaN;
        sumFreq{k}{j} = mean(Freq,2,'omitnan'); % sum frequency
        sumAmp{k}{j} = mean(Amp,2,'omitnan'); % mean amplitude
        GroupSumF = [GroupSumF;sumFreq{k}{j}];
        GroupSumA = [GroupSumA;sumAmp{k}{j}];

        % for paired t test
        paired_F_A{k}(j,1) = mean(Freq,'all','omitnan');
        paired_F_A{k}(j,2) = mean(Amp,'all','omitnan');

    end
    dff_group{k} = dff_stack; % cells of dff for groups
    zpos_Freq{k} = sum(zposFreq,2,'omitnan'); 
    zneg_Freq{k} = sum(znegFreq,2,'omitnan'); 

    All_GSumF{k} = GroupSumF; % collect all groups
    All_GSumA{k} = GroupSumA;
    All_GSumF_MSS{k}(1,1) = mean(GroupSumF,1); % Mean, SD, SEM for frequency
    All_GSumF_MSS{k}(2,1) = std(GroupSumF,0,1);
    All_GSumF_MSS{k}(3,1) = All_GSumF_MSS{k}(2) ./ size(GroupSumF,1);

    All_GSumA_MSS{k}(1,1) = mean(All_GSumA{k},1,'omitnan'); % Mean, SD, SEM for amplitude
    All_GSumA_MSS{k}(2,1) = std(All_GSumA{k},0,1,'omitnan');
    All_GSumA_MSS{k}(3,1) = All_GSumA_MSS{k}(2) ./ size(All_GSumA,1);

    emptyNeu = sum(dff_allsamples,2);
    [row,~] = find(emptyNeu > 0);
    dff_allsamples_filtered = dff_allsamples(row,:);
    dff_allsamples_filtered_min = [];
    for f = 1:size(dff_allsamples_filtered,1)
        dff_allsamples_filtered_min(f,:) = mean(reshape(dff_allsamples_filtered(f,1:8100),15,[]),1);
    end
%     f_raster = figure();
%     [xPoints, yPoints] = plotSpikeRaster(logical(dff_allsamples_filtered_min));
%     f_raster.FontSize = 13;
%     XLabel('Time bin by second');
%     YLabel('Number of neurons');


end
emptyNeu = sum(dff_allsamples,2);
[row,~] = find(emptyNeu > 0);
dff_allsamples_filtered = dff_allsamples(row,:);
dff_allsamples_filtered_min = [];
for f = 1:size(dff_allsamples_filtered,1)
    dff_allsamples_filtered_min(f,:) = mean(reshape(dff_allsamples_filtered(f,1:8100),15,[]),1);
end
% [xPoints, yPoints] = plotSpikeRaster(logical(dff_allsamples_filtered_min)); hold on;
% LineFormat.Color = [0 0 0.8];


norm_DMSO = normalize(dff_group{1},'range',[0 1]);
norm_JNJ = normalize(dff_group{2},'range',[0 1]);
norm_TAK = normalize(dff_group{3},'range',[0 1]);
norm_MK = normalize(dff_group{4},'range',[0 1]);
norm_MKJNJ = normalize(dff_group{6},'range',[0 1]);
norm_MKTAK = normalize(dff_group{8},'range',[0 1]);

% Probability distribution
% TAK [0.6 0.4 1]; DMSO [0 0 1]; JNJ [0.4 0.8 0]; MK [1 0.5 0.2]
% TAK JNJ
f_hist1 = figure('Name','Probability distribution','Visible','off'); clf
h1 = histogram(norm_DMSO','EdgeColor',[0 0 1],'FaceAlpha',0.7,'Normalization','cdf','displaystyle','stairs','LineWidth',1.5); hold on
h2 = histogram(norm_JNJ','EdgeColor',[0.4 0.8 0],'FaceAlpha',0.7,'Normalization','cdf','displaystyle','stairs','LineWidth',1.5); hold on
h3 = histogram(norm_TAK','EdgeColor',[0.6 0.4 1],'FaceAlpha',0.7,'Normalization','cdf','displaystyle','stairs','LineWidth',1.5); hold on
h3 = gca;
h3.FontName = 'Arial';
h3.XAxis.FontSize = 13;
h3.YAxis.FontSize = 13;
h3.XAxis.LineWidth = 1.5;
h3.YAxis.LineWidth = 1.5;
h3.Box = "off";
h3.TickLength = [0 0];
legend({'DMSO','JNJ','TAK'},'LineWidth',1); hold off
ylim([0 1]);
ylabel('Cumulative Distribution of Neurons')
xlabel('Normalized Magnitude of Intensity')
pbaspect([0.45 0.9 0.9]);

% DMSO, MK, MKTAK, MKJNJ
% TAK [0.6 0.4 1]; DMSO [0 0 1]; JNJ [0 1 0.8]; MK [1 0.5 0.2]
f_hist3 = figure('Name','Probability distribution','Visible','off'); clf
h1 = histogram(norm_DMSO','EdgeColor',[0 0 1],'FaceAlpha',0.7,'Normalization','cdf','displaystyle','stairs','LineWidth',1.5); hold on
h2 = histogram(norm_MK','EdgeColor',[1 0.5 0.2],'FaceAlpha',0.7,'Normalization','cdf','displaystyle','stairs','LineWidth',1.5); hold on
h3 = histogram(norm_MKJNJ','EdgeColor',[0.4 0.8 0],'FaceAlpha',0.7,'Normalization','cdf','displaystyle','stairs','LineWidth',1.5); hold on
h4 = histogram(norm_MKTAK','EdgeColor',[0.6 0.4 1],'FaceAlpha',0.7,'Normalization','cdf','displaystyle','stairs','LineWidth',1.5); hold on
h3 = gca;
h3.FontName = 'Arial';
h3.XAxis.FontSize = 13;
h3.YAxis.FontSize = 13;
h3.XAxis.LineWidth = 1.5;
h3.YAxis.LineWidth = 1.5;
h3.Box = "off";
h3.TickLength = [0 0];
legend({'DMSO','MK-801','MKJNJ','MKTAK'},'LineWidth',1); hold off
ylim([0 1]);
ylabel('Cumulative Distribution of Neurons')
xlabel('Normalized Magnitude of Intensity')
pbaspect([0.45 0.9 0.9]);

% % MK MKJNJ
% f_hist2 = figure('Name','Probability distribution','Visible','off'); clf
% h1 = histogram(dff_group{1}','FaceColor',[0 0 1],'FaceAlpha',0.7,'BinWidth',0.5,'Normalization','cdf','displaystyle','stairs'); hold on
% h2 = histogram(dff_group{5}','FaceColor',[1 0.5 0.2],'FaceAlpha',0.7,'BinWidth',0.5,'Normalization','cdf','displaystyle','stairs'); hold on
% h3 = histogram(dff_group{6}','FaceColor',[0.4 0.8 0],'FaceAlpha',0.7,'BinWidth',0.5,'Normalization','cdf','displaystyle','stairs'); hold on
% h3 = gca;
% h3.FontName = 'Arial';
% h3.XAxis.FontSize = 13;
% h3.YAxis.FontSize = 13;
% h3.XAxis.LineWidth = 1.5;
% h3.YAxis.LineWidth = 1.5;
% h3.Box = "off";
% h3.TickLength = [0 0];
% legend({'DMSO','MK-801','MKJNJ'},'LineWidth',1); hold off
% xlim([-2 16]);
% ylim([0 1.1]);
% xlabel('Bins of z-scored fluorescence value');
% ylabel('Probability of neurons');
% pbaspect([0.6 0.7 0.7]);
% 
% % MK MKTAK
% f_hist3 = figure('Name','Probability distribution','Visible','off'); clf
% h1 = histogram(dff_group{1}','FaceColor',[0 0 1],'FaceAlpha',0.7,'BinWidth',0.5,'Normalization','cdf','displaystyle','stairs'); hold on
% h2 = histogram(dff_group{7}','FaceColor',[1 0.5 0.2],'FaceAlpha',0.7,'BinWidth',0.5,'Normalization','cdf','displaystyle','stairs'); hold on
% h3 = histogram(dff_group{8}','FaceColor',[0.6 0.4 1],'FaceAlpha',0.7,'BinWidth',0.5,'Normalization','cdf','displaystyle','stairs'); hold on
% h3 = gca;
% h3.FontName = 'Arial';
% h3.XAxis.FontSize = 13;
% h3.YAxis.FontSize = 13;
% h3.XAxis.LineWidth = 1.5;
% h3.YAxis.LineWidth = 1.5;
% h3.Box = "off";
% h3.TickLength = [0 0];
% legend({'DMSO','MK-801','MKTAK'},'LineWidth',1); hold off
% xlim([-2 16]);
% ylim([0 1.1]);
% xlabel('Bins of z-scored fluorescence value');
% ylabel('Probability of neurons');
% pbaspect([0.6 0.7 0.7]);


if save_mode == 1
hist_path = '/Users/kowteckfong/Desktop/MP InVivoCaImaging/IVCI Images & Videos/Graphs/';
hist_file = fullfile(hist_path,'hist zscore JNJ TAK.jpeg'); % CHANGE name
saveas(f_hist1, hist_file);
hist_path = '/Users/kowteckfong/Desktop/MP InVivoCaImaging/IVCI Images & Videos/Graphs/';
hist_file = fullfile(hist_path,'hist zscore MK MKJNJ.jpeg'); % CHANGE name
saveas(f_hist2, hist_file);
hist_path = '/Users/kowteckfong/Desktop/MP InVivoCaImaging/IVCI Images & Videos/Graphs/';
hist_file = fullfile(hist_path,'hist zscore MK MKTAK.jpeg'); % CHANGE name
saveas(f_hist3, hist_file);
end

%% Student's t test
% groups = {DMSO,JNJ,TAK,MK,MK_JNJ_MK,MK_JNJ,MK_TAK_MK,MK_TAK};
compPair = [1,2; % DMSO v JNJ 1
    1,3; % DMSO v TAK 2
    1,4; % DMSO v MK 3
    5,6; % MK vs MKJNJ 4 
    7,8; % MK vs MKTAK 5
    1,5; % DMSO v MKJNJMK 6
    1,7; % DMSO v MKTAKMK 7
    2,3; % JNJ v TAK 8
    6,8; % MKJNJ v MKTAK 9
    2,6; % JNJ v MKJNJ 10
    3,8]; % TAK v MKTAK 11

% 1v6(MKJNJMK) F: t(1617)=3.3955, p<0.001; A: t(1438)=6.6578, p<0.001
% 1v7(MKTAKMK) F: t(1740)=6.1914, p<0.001; A: t(1866)=-0.9725, p=0.3309
% 1v3(MK) F: t(1331)=4.8696, p<0.001; A: t(1374)=3.3281, p<0.001

for k = 1:size(compPair,1)
[hypF{k},proF{k},ciF{k},statsF{k}] = ttest2(All_GSumF{compPair(k,1)},All_GSumF{compPair(k,2)},'Vartype','unequal');
[hypA{k},proA{k},ciA{k},statsA{k}] = ttest2(All_GSumA{compPair(k,1)},All_GSumA{compPair(k,2)},'Vartype','unequal');

groupLength = max(length(All_GSumF{compPair(k,1)}), length(All_GSumF{compPair(k,2)}));
group1 = nan(groupLength,1);
group2 = nan(groupLength,1);
group1(1:length(All_GSumF{compPair(k,1)})) = All_GSumF{compPair(k,1)};
group2(1:length(All_GSumF{compPair(k,2)})) = All_GSumF{compPair(k,2)};
hedgeg_stats_ftt{k}=mes(group1,group2,'mdbysd','isDep',1,'missVal','listwise','nBoot',1000);
groupLength = max(length(All_GSumA{compPair(k,1)}), length(All_GSumA{compPair(k,2)}));
group1 = nan(groupLength,1);
group2 = nan(groupLength,1);
group1(1:length(All_GSumA{compPair(k,1)})) = All_GSumA{compPair(k,1)};
group2(1:length(All_GSumA{compPair(k,2)})) = All_GSumA{compPair(k,2)};
hedgeg_stats_att{k}=mes(group1,group2,'mdbysd','isDep',1,'missVal','listwise','nBoot',1000);
% groupA_freq = All_GSumF{compPair(k,1)}; % Mann Whitney
% groupB_freq = All_GSumF{compPair(k,2)};
% [MW_p{k},MW_h{k},MW_stats{k}] = ranksum(groupA_freq,groupB_freq);
% groupA_amp = All_GSumA{compPair(k,1)};
% groupB_amp = All_GSumA{compPair(k,2)};
% [MW_p_amp{k},MW_h_amp{k},MW_stats_amp{k}] = ranksum(groupA_amp,groupB_amp);
% groupA_zpos = zpos_Freq{compPair(k,1)};
% groupB_zpos = zpos_Freq{compPair(k,2)};
% [MW_p_zpos{k},MW_h_zpos{k},MW_stats_zpos{k}] = ranksum(groupA_zpos,groupB_zpos);
% groupA_zneg = zneg_Freq{compPair(k,1)};
% groupB_zneg = zneg_Freq{compPair(k,2)};
% [MW_p_zneg{k},MW_h_zneg{k},MW_stats_zneg{k}] = ranksum(groupA_zneg,groupB_zneg);

% kw_groups = [groupA;groupB]; % Kruskal Wallis
% kw_groups_by_num = [ones(1,size(groupA,1))*1,ones(1,size(groupB,1))*2];
% kw_labels = {'A','B'};
% [AB_p{k},AB_tbl{k},AB_stats{k}] = kruskalwallis(kw_groups,kw_labels(kw_groups_by_num),'off');
% % [c_MK,m_MK,h_MK,gnames_MK] = multcompare(stats_kw_MKTAKMKJNJ,'CriticalValueType','dunn-sidak');
end

%% paired; Permutation test
% frequency 
paired_groupA_freq = paired_F_A{compPair(4,1)}(:,1); % MK v MKJNJ
paired_groupB_freq = paired_F_A{compPair(4,2)}(:,1);
[pF{1}, odF{1}, esF{1}] = permutationTest(paired_groupA_freq, paired_groupB_freq, 100);
hedgeg_stats_fpaired{1}=mes(paired_groupA_freq,paired_groupB_freq,'mdbysd','isDep',1,'missVal','listwise','nBoot',1000);
paired_groupA_freq = paired_F_A{compPair(5,1)}(:,1); % MK v MKTAK
paired_groupB_freq = paired_F_A{compPair(5,2)}(:,1);
[pF{2}, odF{2}, esF{2}] = permutationTest(paired_groupA_freq, paired_groupB_freq, 100);
hedgeg_stats_fpaired{2}=mes(paired_groupA_freq,paired_groupB_freq,'mdbysd','isDep',1,'missVal','listwise','nBoot',1000);
% amplitude
paired_groupA_amp = paired_F_A{compPair(4,1)}(:,2);
paired_groupB_amp = paired_F_A{compPair(4,2)}(:,2);
[pA{1}, odA{1}, esA{1}] = permutationTest(paired_groupA_amp, paired_groupB_amp, 100);
hedgeg_stats_apaired{1}=mes(paired_groupA_amp,paired_groupB_amp,'mdbysd','isDep',1,'missVal','listwise','nBoot',1000);
paired_groupA_amp = paired_F_A{compPair(5,1)}(:,2);
paired_groupB_amp = paired_F_A{compPair(5,2)}(:,2);
[pA{2}, odA{2}, esA{2}] = permutationTest(paired_groupA_amp, paired_groupB_amp, 100);
hedgeg_stats_apaired{2}=mes(paired_groupA_amp,paired_groupB_amp,'mdbysd','isDep',1,'missVal','listwise','nBoot',1000);

%% Graphs
% TAK [0.6 0.4 1]; DMSO [0 0 1]; JNJ [0.4 0.8 0]; MK [1 0.5 0.2]

%%% TAK JNJ

% freq_matrix = [All_GSumF{1};All_GSumF{3};All_GSumF{2}]; % CHANGE group and row below
% empty_identifiers = [ones(length(All_GSumF{1}),1)*1; ones(length(All_GSumF{3}),1)*2; ones(length(All_GSumF{2}),1)*3];
% identifiers_name = {'DMSO','TAK','JNJ'}; % CHANGE name
% identifiers = identifiers_name(empty_identifiers);
% f_freq = figure('Name','Frequency by minute JNJ and TAK','Visible',visibility);clf
% vs = violinplot(freq_matrix,identifiers', ...
%     'ViolinColor',[1 0.5 0.2; 0.6 0.4 1; 0.4 0.8 0], ... % CHANGE color
%     'HalfViolin', 'full', ...
%     'ShowMean', true, ...
%     'EdgeColor',[0 0 0], ...
%     'LineWidth', 12, ...
%     'ViolinAlpha',0.1, ...
%     'BoxColor',[0 0 0], ...
%     'MedianMarkerSize', 100,...
%     'Width',0.25, ...
%     'DataStyle','none', ...
%     'LineWidth',1, ...
%     'BoxColor', [0 0 0], ...
%     'ShowMedian',true, ...
%     'QuartileStyle','shadow'); hold on
% ylabel('Mean Frequency (per min)', 'FontSize', 14);
% yticks([0:0.2:1]);
% xticks([1 2 3]);
% ylim([0 1.1]);
% set(gca,'FontSize',12);
% % add significance test result
% plot([1,2],[0.9 0.9],'LineWidth',1.5,'Color',[0 0 0]);
% text([1.5],[0.95],'* p=.0371','FontSize',14);
% 
% amp_matrix = [All_GSumA{1};All_GSumA{3};All_GSumA{2}]; % CHANGE group and row below
% empty_identifiers = [ones(length(All_GSumA{1}),1)*1; ones(length(All_GSumA{3}),1)*2; ones(length(All_GSumA{2}),1)*3];
% identifiers_name = {'DMSO','TAK','JNJ'}; % CHANGE name
% identifiers = identifiers_name(empty_identifiers);
% f_amp = figure('Name','Amplitude by minute JNJ and TAK','Visible',visibility);clf
% vs = violinplot(amp_matrix,identifiers', ...
%     'ViolinColor',[0 0 1; 0.6 0.4 1; 0.4 0.8 0], ... % CHANGE color
%     'HalfViolin', 'full', ...
%     'ShowMean', true, ...
%     'EdgeColor',[0 0 0], ...
%     'LineWidth', 12, ...
%     'ViolinAlpha',0.1, ...
%     'BoxColor',[0 0 0], ...
%     'MedianMarkerSize', 100,...
%     'Width',0.25, ...
%     'DataStyle','none', ...
%     'LineWidth',1, ...
%     'BoxColor', [0 0 0], ...
%     'ShowMedian',true, ...
%     'QuartileStyle','shadow'); hold on
% ylabel('Mean Amplitude (per min)', 'FontSize', 14);
% set(gca, 'XTickLabel', get(gca, 'XTickLabel'), 'FontSize', 14);
% set(gca, 'YTickLabel', get(gca, 'YTickLabel'), 'FontSize', 14);
% ylim([0 40]);
% xticks([1 2 3]);
% yticks([1:5:40]);
% set(gca,'FontSize',12);
% % add significance test result
% plot([1,2],[35 35],'LineWidth',1.5,'Color',[0 0 0]);
% text([1.5],[35.8],'* p=.0131','FontSize',14);

uniformColor = [0 0 0.8];

% JNJ v TAK
empty_identifiers = [ones(length(All_GSumF{1}),1)*1; ones(length(All_GSumF{2}),1)*2; ones(length(All_GSumF{3}),1)*3];
identifiers_name = {'DMSO','JNJ','TAK'}; % CHANGE name
identifiers = categorical(identifiers_name(empty_identifiers));
amp_matrix_MK = [All_GSumF{1};All_GSumF{2};All_GSumF{3}];
f_freq = figure('Name','Mean Frequency by minute JNJ TAK','Visible',visibility);clf
bc = boxchart(identifiers,amp_matrix_MK,'Notch','off','JitterOutliers','on','MarkerStyle','.'); hold on;
bc.BoxFaceColor = uniformColor;
bc.MarkerColor = uniformColor;
ylabel('Mean Frequency (per min)', 'FontSize', 14);
set(gca,'FontSize',12);
ax = gca;
ax.XAxis.Categories = categorical({'DMSO','JNJ','TAK'});
ylim([0 1]);
plot([1.1,1.9],[0.85 0.85],'LineWidth',1.5,'Color',[0 0 0]);
text(1.5,0.86,'*','FontSize',18,'FontWeight','bold');
plot([1.1,2.9],[0.91 0.91],'LineWidth',1.5,'Color',[0 0 0]);
text(1.9,0.94,'ns','FontSize',14);
ax.FontName = 'Arial';
ax.XAxis.FontSize = 13;
ax.YAxis.FontSize = 13;
ax.XAxis.LineWidth = 1.5;
ax.YAxis.LineWidth = 1.5;
ax.Box = "off";
ax.TickLength = [0 0];
pbaspect([0.5 0.7 0.7]);

empty_identifiers = [ones(length(All_GSumA{1}),1)*1; ones(length(All_GSumA{2}),1)*2; ones(length(All_GSumA{3}),1)*3];
identifiers_name = {'DMSO','JNJ','TAK'}; % CHANGE name
identifiers = categorical(identifiers_name(empty_identifiers));
amp_matrix_MK = [All_GSumA{1};All_GSumA{2};All_GSumA{3}];
f_amp = figure('Name','Mean Amplitude by minute JNJ and TAK','Visible','off');clf
bc = boxchart(identifiers,amp_matrix_MK,'Notch','off','JitterOutliers','on','MarkerStyle','.'); hold on;
bc.BoxFaceColor = uniformColor;
bc.MarkerColor = uniformColor;
ylabel('Mean Amplitude (per min)', 'FontSize', 14);
set(gca,'FontSize',12);
ax = gca;
ax.XAxis.Categories = categorical({'DMSO','JNJ','TAK'});
ylim([-3 35]);
plot([1.1,1.9],[30 30],'LineWidth',1.5,'Color',[0 0 0]);
text(1.5,30.9,'*','FontSize',18,'FontWeight','bold');
% plot([1,2,3],[mean(All_GSumF{1},'all','omitnan'),mean(All_GSumF{2},'all','omitnan'),mean(All_GSumF{3},'all','omitnan')],'-','Color',[0 0 1],'LineWidth',2);
ax.FontName = 'Arial';
ax.XAxis.FontSize = 13;
ax.YAxis.FontSize = 13;
ax.XAxis.LineWidth = 1.5;
ax.YAxis.LineWidth = 1.5;
ax.Box = "off";
ax.TickLength = [0 0];
pbaspect([0.5 0.7 0.7]);

%1v3(DMSOvMK) F: t(1331)=4.8696, p<0.001; A: t(1374)=3.3281, p<0.001

% DMSO v MK
empty_identifiers = [ones(length(All_GSumF{1}),1)*1; ones(length(All_GSumF{4}),1)*2];
identifiers_name = {'DMSO','MK'}; % CHANGE name
identifiers = categorical(identifiers_name(empty_identifiers));
amp_matrix_MK = [All_GSumF{1};All_GSumF{4}];
f_freq = figure('Name','Mean Frequency by minute DMSO MK','Visible',visibility);clf
bc = boxchart(identifiers,amp_matrix_MK,'Notch','off','JitterOutliers','on','MarkerStyle','.'); hold on;
bc.BoxFaceColor = uniformColor;
bc.MarkerColor = uniformColor;
ylabel('Mean Frequency (per min)', 'FontSize', 14);
set(gca,'FontSize',12);
ax = gca;
ax.XAxis.Categories = categorical({'DMSO','MK'});
plot([1.1,1.9],[0.85 0.85],'LineWidth',1.5,'Color',[0 0 0]);
text(1.4,0.87,'***','FontSize',18,'FontWeight','bold');
% plot([1,2,3],[mean(All_GSumF{1},'all','omitnan'),mean(All_GSumF{2},'all','omitnan'),mean(All_GSumF{3},'all','omitnan')],'-','Color',[0 0 1],'LineWidth',2);
ax.FontName = 'Arial';
ax.XAxis.FontSize = 13;
ax.YAxis.FontSize = 13;
ax.XAxis.LineWidth = 1.5;
ax.YAxis.LineWidth = 1.5;
ax.Box = "off";
ax.TickLength = [0 0];
pbaspect([0.5 0.7 0.7]);

empty_identifiers = [ones(length(All_GSumA{1}),1)*1; ones(length(All_GSumA{4}),1)*2];
identifiers_name = {'DMSO','MK'}; % CHANGE name
identifiers = categorical(identifiers_name(empty_identifiers));
amp_matrix_MK = [All_GSumA{1};All_GSumA{4}];
f_amp = figure('Name','Mean Amplitude by minute DMSO MK','Visible','off');clf
bc = boxchart(identifiers,amp_matrix_MK,'Notch','off','JitterOutliers','on','MarkerStyle','.'); hold on;
bc.BoxFaceColor = uniformColor;
bc.MarkerColor = uniformColor;
ylabel('Mean Amplitude (per min)', 'FontSize', 14);
set(gca,'FontSize',12);
ax = gca;
ax.XAxis.Categories = categorical({'DMSO','MK'});
ylim([-3 35]);
plot([1.1,1.9],[30 30],'LineWidth',1.5,'Color',[0 0 0]);
text(1.4,30.8,'***','FontSize',18,'FontWeight','bold');
% plot([1,2,3],[mean(All_GSumF{1},'all','omitnan'),mean(All_GSumF{2},'all','omitnan'),mean(All_GSumF{3},'all','omitnan')],'-','Color',[0 0 1],'LineWidth',2);
ax.FontName = 'Arial';
ax.XAxis.FontSize = 13;
ax.YAxis.FontSize = 13;
ax.XAxis.LineWidth = 1.5;
ax.YAxis.LineWidth = 1.5;
ax.Box = "off";
ax.TickLength = [0 0];
pbaspect([0.5 0.7 0.7]);

if save_mode == 1
    violin_path = '/Users/kowteckfong/Desktop/MP InVivoCaImaging/IVCI Images & Videos/Graphs/';
    violin_freq = fullfile(violin_path,'boxplot freq TAK JNJ.jpeg') % CHANGE name
    violin_amp = fullfile(violin_path,'boxplot amp TAK JNJ.jpeg') % CHANGE name
    saveas(f_freq, violin_freq);
    saveas(f_amp, violin_amp);
elseif save_mode == 0
    disp('Not saving');
end

% MK JNJ TAK
empty_identifiers = [ones(length(All_GSumF{4}),1)*1; ones(length(All_GSumF{6}),1)*2; ones(length(All_GSumF{8}),1)*3];
identifiers_name = {'MK','MKJNJ','MKTAK'}; % CHANGE name
identifiers = categorical(identifiers_name(empty_identifiers));
amp_matrix_MK = [All_GSumF{4};All_GSumF{6};All_GSumF{8}];
f_freq_MKJNJTAK = figure('Name','Mean Frequency by minute MK JNJ TAK','Visible',visibility);clf
bc = boxchart(identifiers,amp_matrix_MK,'Notch','off','JitterOutliers','on','MarkerStyle','.'); hold on;
bc.BoxFaceColor = uniformColor;
bc.MarkerColor = uniformColor;
ylabel('Mean Frequency (per min)', 'FontSize', 14);
set(gca,'FontSize',12);
ax = gca;
ax.XAxis.Categories = categorical({'MK','MKJNJ','MKTAK'});
ax.FontName = 'Arial';
ax.XAxis.FontSize = 13;
ax.YAxis.FontSize = 13;
ax.XAxis.LineWidth = 1.5;
ax.YAxis.LineWidth = 1.5;
ax.Box = "off";
ax.TickLength = [0 0];
pbaspect([0.5 0.7 0.7]);
plot([1.1,1.9],[0.9 0.9],'LineWidth',1.5,'Color',[0 0 0]);
text(1.5,0.91,'*','FontSize',18,'FontWeight','bold');
plot([1.1,2.9],[0.94 0.94],'LineWidth',1.5,'Color',[0 0 0]);
text(2,0.95,'*','FontSize',18,'FontWeight','bold'); hold off;
return


% MK JNJ
empty_identifiers = [ones(length(All_GSumF{5}),1)*1; ones(length(All_GSumF{6}),1)*2];
identifiers_name = {'MK','MKJNJ'}; % CHANGE name
identifiers = categorical(identifiers_name(empty_identifiers));
amp_matrix_MK = [All_GSumF{5};All_GSumF{6}];
f_freq_MKJNJTAK = figure('Name','Mean Frequency by minute MK JNJ','Visible',visibility);clf
bc = boxchart(identifiers,amp_matrix_MK,'Notch','off','JitterOutliers','on','MarkerStyle','.'); hold on;
bc.BoxFaceColor = uniformColor;
bc.MarkerColor = uniformColor;
ylabel('Mean Frequency (per min)', 'FontSize', 14);
set(gca,'FontSize',12);
ax = gca;
ax.XAxis.Categories = categorical({'MK','MKJNJ'});
plot([1.1,1.9],[0.85 0.85],'LineWidth',1.5,'Color',[0 0 0]);
text(1.5,0.88,'*','FontSize',18,'FontWeight','bold');
% plot([3.1,3.9],[0.85 0.85],'LineWidth',1.5,'Color',[0 0 0]);
% text(3.5,0.88,'* p<.0198','FontSize',14);
% plot([1,2,3],[mean(All_GSumF{1},'all','omitnan'),mean(All_GSumF{2},'all','omitnan'),mean(All_GSumF{3},'all','omitnan')],'-','Color',[0 0 1],'LineWidth',2);
ax.FontName = 'Arial';
ax.XAxis.FontSize = 13;
ax.YAxis.FontSize = 13;
ax.XAxis.LineWidth = 1.5;
ax.YAxis.LineWidth = 1.5;
ax.Box = "off";
ax.TickLength = [0 0];
pbaspect([0.5 0.7 0.7]);

empty_identifiers = [ones(length(All_GSumA{5}),1)*1; ones(length(All_GSumA{6}),1)*2];
identifiers_name = {'MK','MKJNJ'}; % CHANGE name
identifiers = categorical(identifiers_name(empty_identifiers));
amp_matrix_MK = [All_GSumA{5};All_GSumA{6}];
f_amp_MKJNJTAK = figure('Name','Mean Amplitude by minute MK JNJ','Visible','off');clf
bc = boxchart(identifiers,amp_matrix_MK,'Notch','off','JitterOutliers','on','MarkerStyle','.'); hold on;
bc.BoxFaceColor = uniformColor;
bc.MarkerColor = uniformColor;
ylabel('Mean Amplitude (per min)', 'FontSize', 14);
set(gca,'FontSize',12);
ax = gca;
ax.XAxis.Categories = categorical({'MK','MKJNJ'});
% plot([3.1,3.9],[25 25],'LineWidth',1.5,'Color',[0 0 0]);
% text(3.5,28,'* p=.0099','FontSize',14);
% plot([1,2,3],[mean(All_GSumF{1},'all','omitnan'),mean(All_GSumF{2},'all','omitnan'),mean(All_GSumF{3},'all','omitnan')],'-','Color',[0 0 1],'LineWidth',2);
ax.FontName = 'Arial';
ax.XAxis.FontSize = 13;
ax.YAxis.FontSize = 13;
ax.XAxis.LineWidth = 1.5;
ax.YAxis.LineWidth = 1.5;
ax.Box = "off";
ax.TickLength = [0 0];
pbaspect([0.5 0.7 0.7]);

% MK TAK
empty_identifiers = [ones(length(All_GSumF{7}),1)*1;ones(length(All_GSumF{8}),1)*2];
identifiers_name = {'MK','MKTAK'}; % CHANGE name
identifiers = categorical(identifiers_name(empty_identifiers));
amp_matrix_MK = [All_GSumF{7};All_GSumF{8}];
f_freq_MKJNJTAK = figure('Name','Mean Frequency by minute MK TAK','Visible',visibility);clf
bc = boxchart(identifiers,amp_matrix_MK,'Notch','off','JitterOutliers','on','MarkerStyle','.'); hold on;
bc.BoxFaceColor = uniformColor;
bc.MarkerColor = uniformColor;
ylabel('Mean Frequency (per min)', 'FontSize', 14);
set(gca,'FontSize',12);
ax = gca;
ax.XAxis.Categories = categorical({'MK','MKTAK'});
% plot([1.1,1.9],[0.85 0.85],'LineWidth',1.5,'Color',[0 0 0]);
% text(1.5,0.88,'* p=.0099','FontSize',14);
plot([1.1,1.9],[0.85 0.85],'LineWidth',1.5,'Color',[0 0 0]);
text(1.5,0.88,'*','FontSize',18,'FontWeight','bold');
% plot([1,2,3],[mean(All_GSumF{1},'all','omitnan'),mean(All_GSumF{2},'all','omitnan'),mean(All_GSumF{3},'all','omitnan')],'-','Color',[0 0 1],'LineWidth',2);
ax.FontName = 'Arial';
ax.XAxis.FontSize = 13;
ax.YAxis.FontSize = 13;
ax.XAxis.LineWidth = 1.5;
ax.YAxis.LineWidth = 1.5;
ax.Box = "off";
ax.TickLength = [0 0];
pbaspect([0.5 0.7 0.7]);

empty_identifiers = [ones(length(All_GSumA{7}),1)*1; ones(length(All_GSumA{8}),1)*2];
identifiers_name = {'MK','MKTAK'}; % CHANGE name
identifiers = categorical(identifiers_name(empty_identifiers));
amp_matrix_MK = [All_GSumA{7};All_GSumA{8}];
f_amp_MKJNJTAK = figure('Name','Mean Amplitude by minute MK TAK','Visible','off');clf
bc = boxchart(identifiers,amp_matrix_MK,'Notch','off','JitterOutliers','on','MarkerStyle','.'); hold on;
bc.BoxFaceColor = uniformColor;
bc.MarkerColor = uniformColor;
ylabel('Mean Amplitude (per min)', 'FontSize', 14);
set(gca,'FontSize',12);
ax = gca;
ax.XAxis.Categories = categorical({'MK','MKTAK'});
plot([1.1,1.9],[25 25],'LineWidth',1.5,'Color',[0 0 0]);
text(1.5,25.3,'*','FontSize',18,'FontWeight','bold');
% plot([1,2,3],[mean(All_GSumF{1},'all','omitnan'),mean(All_GSumF{2},'all','omitnan'),mean(All_GSumF{3},'all','omitnan')],'-','Color',[0 0 1],'LineWidth',2);
ax.FontName = 'Arial';
ax.XAxis.FontSize = 13;
ax.YAxis.FontSize = 13;
ax.XAxis.LineWidth = 1.5;
ax.YAxis.LineWidth = 1.5;
ax.Box = "off";
ax.TickLength = [0 0];
pbaspect([0.5 0.7 0.7]);


if save_mode == 1 
    violin_path = '/Users/kowteckfong/Desktop/MP InVivoCaImaging/IVCI Images & Videos/Graphs/';
    violin_freq = fullfile(violin_path,'box freq MK TAK JNJ.jpeg') % CHANGE name
    violin_amp = fullfile(violin_path,'box amp MK TAK JNJ.jpeg') % CHANGE name
    saveas(f_freq_MKJNJTAK, violin_freq);
    saveas(f_amp_MKJNJTAK, violin_amp);
elseif save_mode == 0
    disp('Not saving');
end

% freq_table = table(identifiers',freq_matrix, 'VariableNames',{'Identifiers','Values'});
% estimation_path = '/Users/kowteckfong/Desktop/MP InVivoCaImaging/IVCI Images & Videos/Graphs/';
% estimation_name = fullfile(estimation_path, 'freq_table.csv');
% writetable(freq_table,estimation_name);