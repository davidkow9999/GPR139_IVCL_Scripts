clear;clc
directoryPath = '/Users/kowteckfong/Desktop/MP InVivoCaImaging/IVCI Images & Videos/';
matFiles = dir(fullfile(directoryPath, 'Raw dff zscore/','*raw dff zscore.mat'));
matFiles_locations = dir(fullfile(directoryPath, 'all data/','* all data.mat'));
cd '/Users/kowteckfong/Desktop/MP InVivoCaImaging/IVCI Images & Videos/Raw dff zscore/';
addpath '/Users/kowteckfong/Desktop/MP InVivoCaImaging/IVCI Coding Scripts'

save_rng = '/Users/kowteckfong/Desktop/MP InVivoCaImaging/IVCI Images & Videos/graphpad_calciumimaging/this_rng.mat';
current_rng = rng;
save(save_rng,'current_rng');

DMSO = [1,16,24,28,30,49];
JNJ = [3,26,32,39,51];
MK_JNJ = [5,13,18,34,41,53,61];
MK_JNJ_MK = [7,15,20,36,43,55,63];
TAK = [11,37,47,59,67];
MK_TAK = [8,21,44,56,64]; %% 44 and 64 removed due to complete inactivity, ori: [8,21,44,56,64]
MK_TAK_MK = [10,23,46,58,66];
MK = [MK_JNJ_MK,MK_TAK_MK];

groups = {DMSO,JNJ,TAK,MK,MK_JNJ_MK,MK_JNJ,MK_TAK_MK,MK_TAK}; 
% TAK [0.6 0.4 1]; DMSO [0 0 1]; JNJ [0 1 0.8]; MK [1 0.5 0.2]
visibility = 'off';
save_mode = 0; % 1 to save, 0 to not save

%% Extract DFF, Frequency and (not using) Amplitude
% COMMENT: Normalize by DMSO?
for k = 1:size(groups,2) % all groups
    GroupSumF = [];
    GroupSumA = [];
    dff_stack = [];
    dff_allsamples = [];
    GroupSumF_bylarvae = [];
    GroupSumF_forscatter = [];
    Grouplocation_forscatter = [];
    for j = 1:size(groups{k},2) % all samples in a group
        Freq = [];
        DFF_min = [];
        ave_act = [];
        load(matFiles(groups{k}(j)).name); % drug state
        location_onelarvae{k}{j}(:,1) = normalize(allData{7,2}{1},'range',[1 1000]); % locations of all neurons for each larvae for each group
        location_onelarvae{k}{j}(:,2) = normalize(allData{7,2}{2},'range',[1 1000]);
        Grouplocation_forscatter = [Grouplocation_forscatter;location_onelarvae{k}{j}];
        dff_frame = allData{3,2}; % delta f over f0, from imported file
        dffFreq_frame = allData{5,2}{1,1}; % z raw frequency, from imported file
        trace_dff{k}{j} = dffFreq_frame; % for inspection: frequency of all larvae
        dffFreq_min = dffFreq_frame;
%         dffFreq_min = reshape(dffFreq_frame(i,1:8100),15,[]);
        dff_allsamples = [dff_allsamples;dffFreq_min(:,1:8100)]; % 9 minutes of all samples
        dffAmp = allData{6,2}; % z raw amplitude
        for i = 1:size(dffFreq_frame,1)
            DFF_min(i,:) = mean(reshape(dff_frame(i,1:8100),15,[]),1);
            ave_act(i,:) = trapz(DFF_min(i,:));
            Freq(i,:) = mean(reshape(dffFreq_frame(i,1:8100),15,[]),1,'omitnan'); % mean frequency by min
            peak_indx = find(Freq > 0);
            Freq(peak_indx) = 1;
            Amp(i,:) = mean(reshape(dffAmp(i,1:8100),15,[]),1,'omitnan');
        end
        ave_act_group{k}{j} = ave_act; 
        trace_dff_min{k}{j} = DFF_min;
        dff_stack = [dff_stack; DFF_min];
        sumFreq{k}{j} = mean(Freq,2,'omitnan'); % sum frequency, larvae were separately stored for each group
        sumAmp{k}{j} = mean(Amp,2,'omitnan'); % mean amplitude
        sumFreq_forscatter = sum(Freq,2); % not by larvae
        Freq_bylarvae{k}{j} = sum(Freq,2);
        GroupSumF = [GroupSumF;sumFreq{k}{j}]; % storing all larvae of a group together
        GroupSumF_forscatter = [GroupSumF_forscatter;sumFreq_forscatter];
        GroupSumF_bylarvae = [GroupSumF_bylarvae;Freq_bylarvae{k}{j}];
        GroupSumA = [GroupSumA;sumAmp{k}{j}];

        % for paired t test
        paired_F_A{k}(j,1) = mean(Freq,'all','omitnan');
        paired_F_A{k}(j,2) = mean(Amp,'all','omitnan');

    end
    dff_group{k} = dff_stack'; % cells of dff for groups

%     All_GSumF{k} = GroupSumF; % collect all groups
    All_Glocation_forscatter{k} = Grouplocation_forscatter;
    All_GSUMF_forscatter{k} = GroupSumF_forscatter;
    All_GSumF{k} = GroupSumF_bylarvae;
    All_GSumA{k} = GroupSumA;
    All_GSumF_MSS{k}(1,1) = mean(GroupSumF,1); % Mean, SD, SEM for frequency
    All_GSumF_MSS{k}(2,1) = std(GroupSumF,0,1);
    All_GSumF_MSS{k}(3,1) = All_GSumF_MSS{k}(2) ./ size(GroupSumF,1);

    All_GSumA_MSS{k}(1,1) = mean(All_GSumA{k},1,'omitnan'); % Mean, SD, SEM for amplitude
    All_GSumA_MSS{k}(2,1) = std(All_GSumA{k},0,1,'omitnan');
    All_GSumA_MSS{k}(3,1) = All_GSumA_MSS{k}(2) ./ size(All_GSumA,1);

%     emptyNeu = sum(dff_allsamples,2);
%     [row,~] = find(emptyNeu > 0);
%     dff_allsamples_filtered = dff_allsamples(row,:);
    for f = 1:size(dff_allsamples,1)
        dff_allsamples_filtered_min{k}(f,:) = mean(reshape(dff_allsamples(f,1:8100),15,[]),1); % time series all samples
        dff_meansample_filtered_min{k} = mean(dff_allsamples_filtered_min{k},1); % time series mean 
    end
end

%% raster plot
% figure();
% imagesc(dff_group{1}'); hold on; colorbar;
% figure();
% plot(ten_dff);

%% frequency
for k = 1:8
    % groups = {DMSO,JNJ,TAK,MK,MK_JNJ_MK,MK_JNJ,MK_TAK_MK,MK_TAK}; 
    % ANOVA
    % All_GSumF{k}

end

%% average activity
for k = 1:8
    mean_aotc_all = [];
    for j = 1:size(ave_act_group{k},2)
        mean_aotc{k}{j} = mean(ave_act_group{k}{j},'all');
        mean_aotc_all = [mean_aotc_all;mean_aotc{k}{j}];
    end
    mean_aotc_group{k} = mean_aotc_all;
end

% %% focality index (NOT NEEDED)
% % https://doi.org/10.1371/journal.pbio.0050178
% % DOI: 10.1126/sciadv.aaz3173
% for k = 1:8
%     focality_all = [];
%     for j = 1:size(sumFreq{k},2)
%         % left hb
%         left_neu_indx = find(location_onelarvae{k}{j}(:,1) < 450); 
%         left_neu_coor = location_onelarvae{k}{j}(left_neu_indx,:);
% %         norm_left_nc = normalize(left_neu_coor,'range');
%         left_euclidean{k}{j} = pdist(left_neu_coor,'seuclidean');
%         left_freq = sumFreq{k}{j}(left_neu_indx);
%         t_num_neu = size(left_freq);
%         [lar_neu_freq,lar_neu_idnx] = sort(left_freq,'descend');
%         left_top_ten_freq{k}{j} = lar_neu_freq(1:ceil(t_num_neu*0.1));
%         norm_freq = left_top_ten_freq{k}{j} / sum(left_top_ten_freq{k}{j});
%         top_ten_indx = lar_neu_idnx(1:ceil(t_num_neu*0.1));
%         selected_left_coor{k}{j} = left_neu_coor(top_ten_indx,:);
% %         norm_selected_lc = normalize(selected_left_coor{k}{j},'range');
%         left_ave_euclidean{k}{j} = pdist(selected_left_coor{k}{j},'seuclidean');
%         left_focality_indx = 1 - ((mean(left_ave_euclidean{k}{j},'all').*left_top_ten_freq{k}{j})./ mean(left_euclidean{k}{j},'all'));
% 
%         % right hb
%         right_neu_indx = find(location_onelarvae{k}{j}(:,1) > 550); 
%         right_neu_coor = location_onelarvae{k}{j}(right_neu_indx,:);
% %         norm_right_nc = normalize(right_neu_coor,'range');
%         right_euclidean = pdist(right_neu_coor,'seuclidean');
%         right_freq = sumFreq{k}{j}(right_neu_indx);
%         t_num_neu = size(right_freq);
%         [lar_neu_freq,lar_neu_idnx] = sort(right_freq,'descend');
%         right_top_ten_freq{k}{j} = lar_neu_freq(1:ceil(t_num_neu*0.1));
%         norm_freq = right_top_ten_freq{k}{j} / sum(right_top_ten_freq{k}{j});
%         top_ten_indx = lar_neu_idnx(1:ceil(t_num_neu*0.1));
%         selected_right_coor{k}{j} = right_neu_coor(top_ten_indx,:);
% %         norm_selected_rc = normalize(selected_right_coor{k}{j},'range');
%         right_ave_euclidean = pdist(selected_right_coor{k}{j},'seuclidean'); %% remove normalizd with intensity
%         
%         right_focality_indx = 1 - ((mean(right_ave_euclidean,'all').*norm_freq)./ mean(right_euclidean,'all'));
%         
%         focality_all = [focality_all, mean([right_focality_indx;left_focality_indx],'all')];
%     end
%     focality_grp{k} = focality_all;
% end

%% linear regression and corr based on dff - traced_dff_min or dff_group

% DONT BECAUSE U DIDN"T MEASURE REPEATEDLY
% two_groups_lr_matrix = [dff_group{2},dff_group{3}];
% two_groups_lr_x = [1:size(dff_group{2},2),1:size(dff_group{3},2)];
% group_labels = [zeros(1,size(dff_group{2},2)),ones(1,size(dff_group{3},2))];
% X_lr_x = [two_groups_lr_x;group_labels];
% 
% jnjtak_tbl_lr = table(two_groups_lr_x', mean(two_groups_lr_matrix,1)', group_labels', ...
%             'VariableNames', {'Neurons', 'TimeSeries', 'Group'});
% 
% jnj_tak_mdl = fitlm(jnjtak_tbl_lr,'TimeSeries ~ Neurons + Group');
% 
[rho_jnjtak_dff,pval_jnjtak_dff] = corr([dff_group{1},dff_group{2},dff_group{3}]);
[rho_mkjnjmktak_dff,pval_mkjnjmktak_dff] = corr([dff_group{4},dff_group{6},dff_group{8}]);
% 
% figure();h = plot(jnj_tak_mdl); hold on
% delete(h(1));delete(h(3));
% scatter(mean(dff_group{2},1),mean(dff_group{3},1));

% cd '/Users/kowteckfong/Desktop/MP InVivoCaImaging/IVCI Images & Videos/graphpad_calciumimaging/';
% f_corr_dmso_dff = figure('Position', [100, 100, 1000, 1000]); % [left, bottom, width, height]
% imagesc(rho_jnjtak_dff); hold on
% colorbar;
% colormap('hsv');
% xlabel('Neurons', 'FontName', 'Arial', 'FontSize', 30, 'FontWeight', 'bold');
% ylabel('Neurons', 'FontName', 'Arial', 'FontSize', 30, 'FontWeight', 'bold');
% set(gca, 'TickDir', 'none','FontSize', 20, 'FontName', 'Arial');
% group_boundary = size(dff_group{1},2) + 0.5; 
% total_end = size([dff_group{1},dff_group{2},dff_group{3}],2)
% line([group_boundary, group_boundary], [0, total_end+1], 'Color', 'k', 'LineWidth', 2); % Vertical line
% line([0, total_end+1], [group_boundary, group_boundary], 'Color', 'k', 'LineWidth', 2); % Horizontal line
% group_boundary = size(dff_group{1},2)+ size(dff_group{2},2) + 0.5;
% line([group_boundary, group_boundary], [0, total_end+1], 'Color', 'k', 'LineWidth', 2); % Vertical line
% line([0, total_end+1], [group_boundary, group_boundary], 'Color', 'k', 'LineWidth', 2); % Horizontal line
% exportgraphics(f_corr_dmso_dff, 'corr_dff_dmsojnjtak.pdf', 'Append', false);
% 
% f_corr_mk_dff = figure('Position', [100, 100, 1000, 1000]); % [left, bottom, width, height]
% imagesc(rho_mkjnjmktak_dff); hold on
% colorbar;
% colormap('hsv');
% xlabel('Neurons', 'FontName', 'Arial', 'FontSize', 30, 'FontWeight', 'bold');
% ylabel('Neurons', 'FontName', 'Arial', 'FontSize', 30, 'FontWeight', 'bold');
% set(gca, 'TickDir', 'none','FontSize', 20, 'FontName', 'Arial');
% group_boundary = size(dff_group{4},2) + 0.5; 
% total_end = size([dff_group{4},dff_group{6},dff_group{6}],2)
% line([group_boundary, group_boundary], [0, total_end+1], 'Color', 'k', 'LineWidth', 2); % Vertical line
% line([0, total_end+1], [group_boundary, group_boundary], 'Color', 'k', 'LineWidth', 2); % Horizontal line
% group_boundary = size(dff_group{4},2)+ size(dff_group{6},2) + 0.5;
% line([group_boundary, group_boundary], [0, total_end+1], 'Color', 'k', 'LineWidth', 2); % Vertical line
% line([0, total_end+1], [group_boundary, group_boundary], 'Color', 'k', 'LineWidth', 2); % Horizontal line
% exportgraphics(f_corr_mk_dff, 'corr_dff_mkmkjnjmktak.pdf', 'Append', false);


%% cluster topography, fidelity, selectivity
% DOI: 10.1016/j.cub.2014.01.015
% pairwise correlation between activity and distance of pair of neurons
for k = 1:8
    for j = 1:size(location_onelarvae{k},2)
%         left_hb_neuron_ind = find(location_onelarvae{k}{j}(:,1) < 450);
%         left_hb_neuron = location_onelarvae{k}{j}(left_hb_neuron_ind,1);
%         can use left_hb_neuron substitute location_onelarvae{k}{j} if
%         wanna check just one side of hb
        relevant_euclidean{k}{j} = double(pdist2(location_onelarvae{k}{j}, location_onelarvae{k}{j}, 'euclidean') / 1414.21 * 250.58); % achieved 250.33uM max
        [sorted_rev_euclidean{k}{j}, sorted_rev_euclidean_indx{k}{j}] = sort(relevant_euclidean{k}{j}(1,2:end)); % remove the first one, non relevant
        round_up{k}{j} = ceil(sorted_rev_euclidean{k}{j} / 20) * 20;
        [coef_pair_neuact{k}{j},prob_pair_neuact{k}{j}] = corrcoef(trace_dff_min{k}{j}');
    end
end

coef_disp = coef_pair_neuact{1}{5}(1,2:end);
% figure();plot(sorted_rev_euclidean{1}{5},coef_disp(sorted_rev_euclidean_indx{1}{5}));

x_interval = 10; % 10uM
for o = 1:size(sorted_rev_euclidean,2)
    total_unique_x_coor = 0;
    for p = 1:size(sorted_rev_euclidean{o},2) % sample
        x_coor_interval = ceil(sorted_rev_euclidean{o}{p} /x_interval) *x_interval; % previously * 5 / 5          
        total_unique_x_coor(p) = length(unique(x_coor_interval));
    end
    max_unique = max(total_unique_x_coor);
    x_length_sample = 0:x_interval:x_interval*max_unique; % change after changing x_coor_round_half!!
    int_by_xcoor_byplane = zeros(1,length(x_length_sample));
    for p = 1:size(sorted_rev_euclidean{o},2) % sample
        x_coor_forsample = ceil((sorted_rev_euclidean{o}{p} /x_interval)) *x_interval; % previously * 5 / 5
        int_by_point = 0;
        for i = 1:length(x_length_sample)
            ind_zero_one = x_coor_forsample == x_length_sample(i);
            ind_by_point = find(ind_zero_one == 1);

            if all(ind_by_point == 0)
                int_by_point(i) = NaN;
%                 int_by_point(i) = 0;
            else
                sam_pl_int = coef_pair_neuact{o}{p}(1,2:end);
                int_by_point(i) = mean(sam_pl_int(ind_by_point),'all');
            end
        end
      

        int_by_xcoor_byplane(p+1,:) = movmean(int_by_point,1,'omitnan');
%         int_by_xcoor_byplane(p+1,:) = int_by_point;
        
    end
    int_by_xcoor_bysample{o} = int_by_xcoor_byplane(2:end,:);
end


% {DMSO,JNJ,TAK,MK,MK_JNJ_MK,MK_JNJ,MK_TAK_MK,MK_TAK};
% figure();plot(mean(int_by_xcoor_bysample{1},1,'omitnan')); hold on
% plot(mean(int_by_xcoor_bysample{2},1,'omitnan'));
% plot(mean(int_by_xcoor_bysample{3},1,'omitnan'));
% xticks(x_length_sample);
% legend('DMSO','JNJ','TAK')
% 
% figure();plot(mean(int_by_xcoor_bysample{4},1,'omitnan')); hold on
% plot(mean(int_by_xcoor_bysample{6},1,'omitnan'));
% plot(mean(int_by_xcoor_bysample{8},1,'omitnan'));
% xticks(x_length_sample);
% legend('MK','MKJNJ','MKTAK')

%% Dimension reduction
%% all groups together - identify unique effect of compounds
group_size = []; data_reduced_all = []; data_noreduced_all = [];
for k = 1:8
data = dff_group{k};
coef{k} = corr(data);
high_corr{k} = abs(coef{k} > 0.65);
reduced_index = find(~any(triu(high_corr{k},1)) == 1);
data_reduced{k} = data(:, reduced_index);
group_size = [group_size;size(data_reduced{k},2)];
data_reduced_all = [data_reduced_all,data_reduced{k}];
data_noreduced_all = [data_noreduced_all,data];
end

[coeff_pca_group,score_pca_group,tsquared_pca_group,~,explained_pca_group] = pca(data_noreduced_all');
variance_percent = cumsum(explained_pca_group);
variance_idx = find(variance_percent > 99);

Y_tsne_group = tsne(score_pca_group);
% Y_tsne = tsne(score_pca_group(1:variance_idx(1),:));

% groups = {DMSO,JNJ,TAK, MK, MK_JNJ_MK,MK_JNJ, MK_TAK_MK,MK_TAK}; 
% 1-6; 2-5; 3-5; 4-5; 5-5; 6-5; 7-6; 8-5
% determine cluster number using elbow method
clust_num = [6,5,5,5,5,5,6,5];
clust_now = [];
    for f = 1:20
%         [clust_idx_postpca{f}, clust_C_postpca{f}, clust_sumd{f}] = kmeans(score_pca_group(1:variance_idx(1),:),f);
        [clust_idx_postpca_group{f}, clust_C_postpca_group{f}, clust_sumd_group{f}] = kmeans(score_pca_group,f);        
        clust_now = [clust_now,sum(clust_sumd_group{f})]; % for elbow method
    end
    clust_sumd_all_group = clust_now; % for elbow method
    % https://www.mathworks.com/matlabcentral/fileexchange/35094-knee-point
    [~, idx_of_result_group] = knee_pt(1:20,clust_sumd_all_group); % 6, but group is 8

% visualize cluster number using WCSS (sumd)
% figure();plot(1:10,clust_sumd_all_group{8},'-o')

% FIGURE: SAVE THIS
% cd '/Users/kowteckfong/Desktop/MP InVivoCaImaging/IVCI Images & Videos/Raw dff zscore/';
cd '/Users/kowteckfong/Desktop/MP InVivoCaImaging/IVCI Images & Videos/graphpad_calciumimaging/';
colorOrder = [
    0,    0.4470,    0.7410;
    0.8500,    0.3250,    0.0980;
    0.9290,    0.6940,    0.1250;
    0.4940,    0.1840,    0.5560;
    0.4660,    0.6740,    0.1880;
    0.3010,    0.7450,    0.9330;
    0.6350,    0.0780,    0.1840;
    1,         0,         1;  
    ];

% f = figure('Position', [100, 100, 1000, 1000]); % [left, bottom, width, height]
% % axes('Position', [0.1, 0.1, 1, 1]); % [left, bottom, width, height]
% gscatter(Y_tsne_group(:,1),Y_tsne_group(:,2),clust_idx_postpca_group{8},colorOrder,[],14)
% hold on;
% gscatter(Y_tsne_group(:,1),Y_tsne_group(:,2),clust_idx_postpca_group{8},'k','o',6);
% xlabel('tSNE 2', 'FontName', 'Arial', 'FontSize', 30, 'FontWeight', 'bold');
% ylabel('tSNE 1', 'FontName', 'Arial', 'FontSize', 30, 'FontWeight', 'bold');
% legend('Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4' , 'Cluster 5', 'Cluster 6', 'Cluster 7', 'Cluster 8' , ...
%     'FontName', 'Arial', 'FontSize', 20);
% set(gca, 'TickDir', 'none','FontSize', 20, 'FontName', 'Arial','ColorOrder',colorOrder);
% % exportgraphics(f, 'tSNE_all_groups.pdf', 'Append', false);

%% individual group - identify potential habenula subregions
for k = 1:8
    cluster_number = [];
    clust_euclidean_stack = [zeros(1,3)];
    for j = 1:size(trace_dff_min{k},2)
    [coeff_pca{k}{j},score_pca{k}{j},~,~,explained_pca{k}{j}] = pca(trace_dff_min{k}{j});
    variance_percent = cumsum(explained_pca{k}{j});
    variance_idx = find(variance_percent > 95);
    
    Y_tsne{k}{j} = tsne(score_pca{k}{j});
    % Y_tsne{k} = tsne(score_pca{k}(1:variance_idx(1),:));
    
    % groups = {DMSO,JNJ,TAK,MK,MK_JNJ_MK,MK_JNJ,MK_TAK_MK,MK_TAK}; 
    % 1-6; 2-5; 3-5; 4-5; 5-5; 6-8; 7-6; 8-5
    % determine cluster number using elbow method
    clust_num = [6,6,8,5,6,6,7,5];
    clust_now = [];
        for f = 1:20 % determine the number of cluster
    %         [clust_idx_postpca{k}{f}, clust_C_postpca{k}{f}, clust_sumd{k}{f}] = kmeans(score_pca{k}(1:variance_idx(1),:),f);
            [clust_idx_postpca{k}{j}{f}, clust_C_postpca{k}{j}{f}, clust_sumd{k}{j}{f}] = kmeans(score_pca{k}{j},f);
            [s_kclust{k}{j}{f},h_kclust{k}{j}{f}] = silhouette(score_pca{k}{j},clust_idx_postpca{k}{j}{f});
            clust_now = [clust_now,sum(clust_sumd{k}{j}{f})]; % for elbow method
        end
        clust_sumd_all{k}{j} = clust_now; % for elbow method
        % 6	6	8	5	6	6	7	5
        [~, idx_of_result{k}{j}] = knee_pt(1:20,clust_sumd_all{k}{j});
        cluster_number = [cluster_number, idx_of_result{k}{j}];
       
        for count = 1:idx_of_result{k}{j} % go thru all clust determined
            the_wanted_clust_num = clust_idx_postpca{k}{j}{idx_of_result{k}{j}};
            one_clust = find(the_wanted_clust_num == count); % f is always the num of determined clusters
            clust_dff_group_postpca{k}{j}{count} = trace_dff_min{k}{j}(one_clust,:); % can be used i

            clust_loc_group_postpca{k}{j}{count} = location_onelarvae{k}{j}(one_clust,:);
            clust_euclidean = double(pdist2(location_onelarvae{k}{j}(one_clust,:), location_onelarvae{k}{j}(one_clust,:), 'euclidean') / 1414.21 * 250.58);
            no_repeat = triu(clust_euclidean);
            no_repeat_vector = no_repeat(no_repeat ~= 0);
            clust_loc_euclidean{k}{j}{count} = no_repeat_vector;
            m_std_euclidean = [mean(no_repeat_vector,'all'),std(no_repeat_vector,[],'all'),size(no_repeat_vector,1)];
            clust_euclidean_stack = [clust_euclidean_stack;m_std_euclidean];
        end


%         for count = 1:idx_of_result{k}{j} % to store neu dff and loc according to cluster
%             one_clust = find(clust_idx_postpca{k}{j}{f} == count); % f is always the num of determined clusters
%             clust_dff_group_postpca{k}{j}{count} = trace_dff_min{k}{j}(one_clust,:); % can be used i
% %             clust_loc_group_postpca{k}{j}{count} = All_Glocation_forscatter{k}(one_clust,:);
%             clust_loc_group_postpca{k}{j}{count} = location_onelarvae{k}{j}(one_clust,:);
%             clust_euclidean = double(pdist2(location_onelarvae{k}{j}(one_clust,:), location_onelarvae{k}{j}(one_clust,:), 'euclidean') / 1414.21 * 250.58);
%             no_repeat = triu(clust_euclidean);
%             no_repeat_vector = no_repeat(no_repeat ~= 0);
%             clust_loc_euclidean{k}{j}{count} = no_repeat_vector;
%         end
    end
    total_cluster{k} = cluster_number;
    sort_eu_stack = clust_euclidean_stack(2:end,:);
    [~,ind_sorted_e_s] = sort(sort_eu_stack(:,1),'ascend');
    clust_loc_euclidean_allin{k} = sort_eu_stack(ind_sorted_e_s,:);
end
% visualize cluster number using WCSS (sumd)
% figure();plot(1:20,clust_sumd_all{8},'-o')

return

%% scatter graphs 
% JNJ TAK scatter
cd '/Users/kowteckfong/Desktop/MP InVivoCaImaging/IVCI Images & Videos/graphpad_calciumimaging/';
color_range = colormap('hsv');
sum_clust = floor(256/ sum(total_cluster{2}));
c_range = 1:sum_clust:256;
jnj_scatter = figure('Position', [100, 100, 1000, 1000]);

jnj_scatter_ind = [5,5,8,6,7];
ccount = 1;
for j = 1:5
    for i = 1:total_cluster{2}(j)
scatter(clust_loc_group_postpca{2}{j}{i}(:,1),clust_loc_group_postpca{2}{j}{i}(:,2),300, 'filled', ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', [color_range(c_range(ccount+1),:)],'LineWidth',2); hold on
    ccount = ccount + 1;
    end
end
set(gca, 'XTick', [], 'YTick', [])
% exportgraphics(jnj_scatter, 'scatter_jnj2.pdf', 'Append', false);

sum_clust = floor(256/ sum(total_cluster{3}));
c_range = 1:sum_clust:256;
tak_scatter = figure('Position', [100, 100, 1000, 1000]);
tak_scatter_ind = [6,7,6,7,6];
ccount = 1;
for j = 1:5
    for i = 1:total_cluster{3}(j)
scatter(clust_loc_group_postpca{3}{j}{i}(:,1),clust_loc_group_postpca{3}{j}{i}(:,2),300, 'filled', ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', [color_range(c_range(ccount+1),:)],'LineWidth',2); hold on
    ccount = ccount + 1;
    end
end
set(gca, 'XTick', [], 'YTick', [])
% exportgraphics(tak_scatter, 'scatter_tak2.pdf', 'Append', false);

% MKJNJ MKTAK scatter
color_range = colormap('hsv');
sum_clust = floor(256/ sum(total_cluster{6}));
c_range = 1:sum_clust:256;
mkjnj_scatter = figure('Position', [100, 100, 1000, 1000]);
mkjnj_scatter_ind = [9,6,6,6,6,7,5];
ccount = 1;
for j = 1:7
    for i = 1:total_cluster{6}(j)
scatter(clust_loc_group_postpca{6}{j}{i}(:,1),clust_loc_group_postpca{6}{j}{i}(:,2),300, 'filled', ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', [color_range(c_range(ccount+1),:)],'LineWidth',2); hold on
    ccount = ccount + 1;
    end
end
set(gca, 'XTick', [], 'YTick', [])
% exportgraphics(mkjnj_scatter, 'scatter_mkjnj2.pdf', 'Append', false);

sum_clust = floor(256/ sum(total_cluster{8}))-1;
c_range = 1:sum_clust:256;
mktak_scatter = figure('Position', [100, 100, 1000, 1000]);
mktak_scatter_ind = [7,5,8,8,5];
ccount = 1;
for j = 1:5
    for i = 1:total_cluster{8}(j)
scatter(clust_loc_group_postpca{8}{j}{i}(:,1),clust_loc_group_postpca{8}{j}{i}(:,2),300, 'filled', ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', [color_range(c_range(ccount+1),:)],'LineWidth',2); hold on
    ccount = ccount + 1;
    end
end
set(gca, 'XTick', [], 'YTick', [])
% exportgraphics(mktak_scatter, 'scatter_mktak2.pdf', 'Append', false);

return

%% find correlation of mean dff among clusters
% DMSO JNJ TAK
big_cluster_dmso = [];
for k = 1:3
    count = 0;
    for j = 1:size(clust_dff_group_postpca{k},2)
        for c = 1:size(clust_dff_group_postpca{k}{j},2)
            one_clus_mean_dff = mean(clust_dff_group_postpca{k}{j}{c},1);
            count = count + 1;
            big_cluster_dmso = [big_cluster_dmso,one_clus_mean_dff'];
        end
    end
    count_dmso{k} = count;
end
[rho_bigclust_dmso,pval_bigclust_dmso] = corr(big_cluster_dmso);

% MK MKJNJ MKTAK
big_cluster_mk = [];
ind_count = 1;
for k = [3,6,8]
    count = 0;
    for j = 1:size(clust_dff_group_postpca{k},2)
        for c = 1:size(clust_dff_group_postpca{k}{j},2)
            one_clus_mean_dff = mean(clust_dff_group_postpca{k}{j}{c},1);
            count = count + 1;
            big_cluster_mk = [big_cluster_mk,one_clus_mean_dff'];
        end
    end
    count_mk{ind_count} = count;
    ind_count = ind_count + 1;
end
[rho_bigclust_mk,pval_bigclust_mk] = corr(big_cluster_mk);

% % figure for dmso
% f_pc_dmso = figure('Position', [100, 100, 1000, 1000]);
% imagesc(rho_bigclust_dmso); hold on
% colormap('parula');
% group_boundary = count_dmso{1} + 0.5; 
% line([group_boundary, group_boundary], [0, size(big_cluster_dmso,2)+1], 'Color', 'k', 'LineWidth', 2); % Vertical line
% line([0, size(big_cluster_dmso,2)+1], [group_boundary, group_boundary], 'Color', 'k', 'LineWidth', 2); % Horizontal line
% group_boundary = count_dmso{1}+ count_dmso{2} + 0.5;
% line([group_boundary, group_boundary], [0, size(big_cluster_dmso,2)+1], 'Color', 'k', 'LineWidth', 2); % Vertical line
% line([0, size(big_cluster_dmso,2)+1], [group_boundary, group_boundary], 'Color', 'k', 'LineWidth', 2); % Horizontal line
% colorbar;
% % xlabel('Clusters', 'FontName', 'Arial', 'FontSize', 30, 'FontWeight', 'bold');
% % ylabel('Clusters', 'FontName', 'Arial', 'FontSize', 30, 'FontWeight', 'bold');
% set(gca, 'TickDir', 'none','FontSize', 20, 'FontName', 'Arial','ColorOrder',colorOrder);
% % exportgraphics(f_pc_dmso, 'pairwise_corr_dmso.pdf', 'Append', false);
% hold off;
% 
% % figure for mk
% f_pc_mk = figure('Position', [100, 100, 1000, 1000]);
% imagesc(rho_bigclust_mk); hold on
% colormap('parula');
% group_boundary = count_mk{1} + 0.5; 
% line([group_boundary, group_boundary], [0, size(big_cluster_mk,2)+1], 'Color', 'k', 'LineWidth', 2); % Vertical line
% line([0, size(big_cluster_mk,2)+1], [group_boundary, group_boundary], 'Color', 'k', 'LineWidth', 2); % Horizontal line
% group_boundary = count_mk{1}+ count_mk{2} + 0.5;
% line([group_boundary, group_boundary], [0, size(big_cluster_mk,2)+1], 'Color', 'k', 'LineWidth', 2); % Vertical line
% line([0, size(big_cluster_mk,2)+1], [group_boundary, group_boundary], 'Color', 'k', 'LineWidth', 2); % Horizontal line
% colorbar;
% set(gca, 'TickDir', 'none','FontSize', 20, 'FontName', 'Arial','ColorOrder',colorOrder);
% exportgraphics(f_pc_mk, 'pairwise_corr_mk.pdf', 'Append', false);
% hold off;

%% identify the clusters x and y coordinates
for k = 1:8
    stack_fish_statsc = []; lr_fish_coor_y = zeros(1,2); lr_fish_dff = zeros(1,2);
    for j = 1:size(clust_loc_group_postpca{k},2) % fish
        stack_statsc = zeros(1,10); together_stats =[]; lr_coor_y = zeros(1,2); lr_coor_y_single =[]; lr_dff = zeros(1,2);
        for h = 1:size(clust_loc_group_postpca{k}{j},2) % cluster
            left_c = find(clust_loc_group_postpca{k}{j}{h}(:,1) < 500);
            right_c = find(clust_loc_group_postpca{k}{j}{h}(:,1) > 500);
            % x coor
            left_coor = double(clust_loc_group_postpca{k}{j}{h}(left_c,1)); % identify x coor for left hb
            right_coor = double(clust_loc_group_postpca{k}{j}{h}(right_c,1));
            % y coor
            left_coor_y = double(clust_loc_group_postpca{k}{j}{h}(left_c,2)); % identify y coor for left hb
            right_coor_y = double(clust_loc_group_postpca{k}{j}{h}(right_c,2));
            % dff
            left_dff = double(clust_dff_group_postpca{k}{j}{h}(left_c,:)); % identify dff for left hb
            right_dff = double(clust_dff_group_postpca{k}{j}{h}(right_c,:));
            
            % mean, std, max, min for x-coor
            mean_leftc = mean(left_coor,'omitnan'); std_leftc = std(left_coor,'omitnan'); max_leftc = max(left_coor); min_leftc = min(left_coor);
            mean_rightc = mean(right_coor,'omitnan'); std_rightc = std(right_coor,'omitnan'); max_rightc = max(right_coor); min_rightc = min(right_coor);
            if isempty(max_leftc) == 1
                max_leftc = NaN;
            end
            if isempty(min_leftc) == 1
                min_leftc = NaN;
            end
            if isempty(max_rightc) == 1
                max_rightc = NaN;
            end
            if isempty(min_rightc) == 1
                min_rightc = NaN;
            end
            % x coor
            together_stats = [mean_leftc,std_leftc,max_leftc,min_leftc, mean_rightc,std_rightc,max_rightc,min_rightc, h,j];
            stack_statsc = [stack_statsc;together_stats];
            % y coor
            lr_coor_y_single = [mean(left_coor_y), mean(right_coor_y)];
            lr_coor_y = [lr_coor_y; lr_coor_y_single];
            % dff
            lr_dff_single = [mean(left_dff,'all'), mean(right_dff,'all')];
            lr_dff = [lr_dff; lr_dff_single];

            
        end
        % x coor
        stack_fish_statsc = [stack_fish_statsc; stack_statsc(2:end,:)];
        % y coor
        lr_fish_coor_y = [lr_fish_coor_y;lr_coor_y(2:end,:)];
        % dff
        lr_fish_dff = [lr_fish_dff;lr_dff(2:end,:)];
    end
    % x coor
    stack_fish_statsc(stack_fish_statsc == 0) = 1;
    clust_coor{k} = stack_fish_statsc;
    % y coor
    clust_y_coor{k} = lr_fish_coor_y(2:end,:);
    % dff
    clust_dff{k} = lr_fish_dff(2:end,:);
end

% % plot scatter
% group_num = 2;
% figure();
% for j = 1:size(clust_coor{group_num},1)
%     scatter(All_Glocation_forscatter{group_num}(:,1),All_Glocation_forscatter{group_num}(:,2),'o','k'); hold on
%     scatter(clust_coor{group_num}(j,1),clust_y_coor{group_num}(j,1),(clust_coor{group_num}(j,2)*5),'filled'); hold on
%     scatter(clust_coor{group_num}(j,5),clust_y_coor{group_num}(j,2),(clust_coor{group_num}(j,6)*5),'filled');
% end
% colormap("jet");
% % plot(200*ones(1,1000),1:1000,'Color','k');
% % plot(800*ones(1,1000),1:1000,'Color','k');
% xlabel('x coordinate'); ylabel('y coordinate'); zlabel('fishes');
% xlim([0 1000]); ylim([0 1000]);


return
%% NO USE: Pairwise, correlation coeff, and linear regression (try bootstrap?)
for k = 1:8
    for j = 1:size(location_onelarvae{k},2)
    %     above_80p = prctile(All_GSUMF_forscatter{k},60);
    %     ind_80p = find(All_GSUMF_forscatter{k} > above_80p);
    %    
    %     below_20p = prctile(All_GSUMF_forscatter{k},40);
    %     ind_20p = find(All_GSUMF_forscatter{k} < below_20p);
    
    %     figure();
    %     scatter(All_Glocation_forscatter{k}(ind_80p,1),All_Glocation_forscatter{k}(ind_80p,2),sumFreq_forscatter{k}(ind_80p)); hold on;
    %     scatter(All_Glocation_forscatter{k}(ind_20p,1),All_Glocation_forscatter{k}(ind_20p,2),sumFreq_forscatter{k}(ind_20p)); hold on;
        
        % sum frequency 
        value_left_hb = find(location_onelarvae{k}{j}(:,1) < 450);
        value_right_hb = find(location_onelarvae{k}{j}(:,1) > 550);
    
        [coef_dist_left{k}{j},prob_dist_left{k}{j}] = corr(Freq_bylarvae{k}{j}(value_left_hb),location_onelarvae{k}{j}(value_left_hb,1));
        [coef_dist_right{k}{j},prob_dist_right{k}{j}] = corr(Freq_bylarvae{k}{j}(value_right_hb),location_onelarvae{k}{j}(value_right_hb,1));

        left_mdl_lr_freq{k}{j} = fitlm(location_onelarvae{k}{j}(value_left_hb,:),Freq_bylarvae{k}{j}(value_left_hb));
        right_mdl_lr_freq{k}{j} = fitlm(location_onelarvae{k}{j}(value_right_hb,:),Freq_bylarvae{k}{j}(value_right_hb));
    
        % time series
        mat_loc = ones(size(trace_dff_min{k}{j},540));
        mat_left_loc = repmat(trace_dff_min{k}{j}(value_left_hb),1,540)';
        mat_right_loc = repmat(trace_dff_min{k}{j}(value_right_hb),1,540)';

        [coef_dist_left_timeseries{k}{j},prob_dist_left_timeseries{k}{j}] = corr(trace_dff_min{k}{j}(value_left_hb,:)',mat_left_loc);
        [coef_dist_right_timeseries{k}{j},prob_dist_right_timeseries{k}{j}] = corr(trace_dff_min{k}{j}(value_right_hb,:)',mat_right_loc);

        left_mdl_lr_dff{k}{j} = fitlm(location_onelarvae{k}{j}(value_left_hb,1),mean(trace_dff_min{k}{j}(value_left_hb,:),2));
        right_mdl_lr_dff{k}{j} = fitlm(location_onelarvae{k}{j}(value_right_hb,1),mean(trace_dff_min{k}{j}(value_right_hb,:),2));
    end
end

% figure();
% [sort_80p,idx_sort] = sort(All_Glocation_forscatter{k}(ind_80p,1)');
% sort_80p_int = All_GSUMF_forscatter{k}(idx_sort);
% plot(sort_80p,sort_80p_int,'green'); hold on;
% 
%     [sort_20p,idx_sort] = sort(All_Glocation_forscatter{k}(ind_20p,1)');
%     sort_20p_int = sumFreq_forscatter{k}(idx_sort);
%     plot(sort_20p,sort_20p_int,'blue'); hold on;
   
