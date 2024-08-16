clear;clc;
rng(42);

% stats for subnuclei difference (v-d)
% stats for v - d for all groups
% generate graphs for subnuclei and v-d all groups 

cd '/Users/kowteckfong/Desktop/MP InVivoCaImaging/pERK data/';
directoryPath = "/Users/kowteckfong/Desktop/MP InVivoCaImaging/pERK data/";
groups = {'DMSO';'TAK';'JNJ';'MKTAK';'MKJNJ';'MK'};
for l = 1:size(groups,1)
    group_path = directoryPath + groups{l} + '/';
    group_folder = dir(group_path);

visibility = 'off';
save_mode = 0;

%% Pool all neurons for a group AND by plane
num_samples = {group_folder.name};
sample_size = size(num_samples,2) - 3;
plane1 = [];
plane2 = [];
plane3 = [];
plane4 = [];
plane5 = [];

for i = 1:sample_size
    terk_matrix = []; % temporary store for a trial
    perk_matrix = [];
    aplane_perk_terk = {}; % temporary store
    current_group_path = group_path + group_folder(i + 3).name;
    group_1 = dir(fullfile(current_group_path, '*.csv')); % list all trials for the group
    k = 1; % to collect neurons separated by each plane
    xyz_locs_terk = [NaN,NaN,NaN];
    xyz_locs_perk = [NaN,NaN,NaN];
    exp_date = [];
    for j = 1:2:5*2 % collecting for all plane
        current_trial = current_group_path + "/" + group_1(j).name; % load for each plane, terk
        terk = readmatrix(current_trial);
        nan_row = any(isnan(terk),2);
        terk_without_nan = terk(~nan_row,:);
        terk_matrix = [terk_matrix; terk_without_nan]; % remove nan
        mean_aplane_terk = mean(terk_without_nan(:,3),'all');

        xyz_locs_terk = [xyz_locs_terk;terk_without_nan(:,5),terk_without_nan(:,6),ones(size(terk_without_nan,1),1)*k]; % xyz locations

        current_trial = current_group_path + "/" + group_1(j+1).name; % load for each plane, perk
        perk = readmatrix(current_trial);
        nan_row = any(isnan(perk),2);
        perk_without_nan = perk(~nan_row,:);
        perk_matrix = [perk_matrix; perk_without_nan]; % remove nan
        aplane_perk_terk{k} = perk_without_nan(:,3) ./ mean_aplane_terk; % pERK/tERK ratio with mean tERK of respective plane

%         xyz_locs_perk = [xyz_locs_perk;perk_without_nan(:,5),perk_without_nan(:,6),ones(size(perk_without_nan,1),1)*k]; % xyz locations
        xyz_locs_perk = [xyz_locs_perk;normalize(perk_without_nan(:,5),'range',[1 7]),normalize(perk_without_nan(:,6),'range',[1 7]),ones(size(perk_without_nan,1),1)*k];

        exp_prop = split(group_1(j).name,'_');
        exp_date = [exp_date;exp_prop{1}];
        k = k + 1;
    end
    
    % assign data to respective group
    org_group = split(group_path,'/');
    if org_group(end-1) == 'DMSO'
        mean_terk = mean(terk_matrix(:,3),'all'); % mean of all terk for each group
        DMSO_perk_terk{i} = perk_matrix(:,3) ./ mean_terk; % pERK/tERK ratio with mean terk from all planes
    
        DMSO_byplane_perk_terk{i} = [aplane_perk_terk{1};aplane_perk_terk{2};aplane_perk_terk{3};
            aplane_perk_terk{4};aplane_perk_terk{5}];

        DMSO_samples_planes{i,1} = aplane_perk_terk{1}; DMSO_samples_planes{i,2} = aplane_perk_terk{2}; DMSO_samples_planes{i,3} = aplane_perk_terk{3}; 
        DMSO_samples_planes{i,4} = aplane_perk_terk{4}; DMSO_samples_planes{i,5} = aplane_perk_terk{5}; 
    
        plane1 = [plane1;aplane_perk_terk{1}]; 
        plane2 = [plane2;aplane_perk_terk{2}];
        plane3 = [plane3;aplane_perk_terk{3}];
        plane4 = [plane4;aplane_perk_terk{4}];
        plane5 = [plane5;aplane_perk_terk{5}];
        DMSO_planes{1} = plane1; DMSO_planes{2} = plane2; DMSO_planes{3} = plane3; %need to redo this, all wrong
        DMSO_planes{4} = plane4; DMSO_planes{5} = plane5;

        DMSO_xyz{i,1} = xyz_locs_terk(2:end,:);
        DMSO_xyz{i,2} = xyz_locs_perk(2:end,:);

        date_DMSO = exp_date;
    elseif org_group(end-1) == 'TAK';
        mean_terk = mean(terk_matrix(:,3),'all'); % mean of all terk for each group
        TAK_perk_terk{i} = perk_matrix(:,3) ./ mean_terk; % pERK/tERK ratio with mean terk from all planes
    
        TAK_byplane_perk_terk{i} = [aplane_perk_terk{1};aplane_perk_terk{2};aplane_perk_terk{3};
            aplane_perk_terk{4};aplane_perk_terk{5}];

        TAK_samples_planes{i,1} = aplane_perk_terk{1}; TAK_samples_planes{i,2} = aplane_perk_terk{2}; TAK_samples_planes{i,3} = aplane_perk_terk{3}; 
        TAK_samples_planes{i,4} = aplane_perk_terk{4}; TAK_samples_planes{i,5} = aplane_perk_terk{5}; 
    
        plane1 = [plane1;aplane_perk_terk{1}];
        plane2 = [plane2;aplane_perk_terk{2}];
        plane3 = [plane3;aplane_perk_terk{3}];
        plane4 = [plane4;aplane_perk_terk{4}];
        plane5 = [plane5;aplane_perk_terk{5}];
        TAK_planes{1} = plane1; TAK_planes{2} = plane2; TAK_planes{3} = plane3; 
        TAK_planes{4} = plane4; TAK_planes{5} = plane5;

        TAK_xyz{i,1} = xyz_locs_terk(2:end,:);
        TAK_xyz{i,2} = xyz_locs_perk(2:end,:);

        date_TAK = exp_date;
    elseif org_group(end-1) == 'JNJ';
        mean_terk = mean(terk_matrix(:,3),'all'); % mean of all terk for each group
        JNJ_perk_terk{i} = perk_matrix(:,3) ./ mean_terk; % USE: collecting for all plane
    
        JNJ_byplane_perk_terk{i} = [aplane_perk_terk{1};aplane_perk_terk{2};aplane_perk_terk{3};
            aplane_perk_terk{4};aplane_perk_terk{5}];

        JNJ_samples_planes{i,1} = aplane_perk_terk{1}; JNJ_samples_planes{i,2} = aplane_perk_terk{2}; JNJ_samples_planes{i,3} = aplane_perk_terk{3}; 
        JNJ_samples_planes{i,4} = aplane_perk_terk{4}; JNJ_samples_planes{i,5} = aplane_perk_terk{5}; 
    
        plane1 = [plane1;aplane_perk_terk{1}]; % USE
        plane2 = [plane2;aplane_perk_terk{2}];
        plane3 = [plane3;aplane_perk_terk{3}];
        plane4 = [plane4;aplane_perk_terk{4}];
        plane5 = [plane5;aplane_perk_terk{5}];
        JNJ_planes{1} = plane1; JNJ_planes{2} = plane2; JNJ_planes{3} = plane3; 
        JNJ_planes{4} = plane4; JNJ_planes{5} = plane5;

        JNJ_xyz{i,1} = xyz_locs_terk(2:end,:);
        JNJ_xyz{i,2} = xyz_locs_perk(2:end,:);

        date_JNJ = exp_date;
    elseif org_group(end-1) == 'MKTAK';
        mean_terk = mean(terk_matrix(:,3),'all'); % mean of all terk for each group
        MKTAK_perk_terk{i} = perk_matrix(:,3) ./ mean_terk; % USE: collecting for all plane
    
        MKTAK_byplane_perk_terk{i} = [aplane_perk_terk{1};aplane_perk_terk{2};aplane_perk_terk{3};
            aplane_perk_terk{4};aplane_perk_terk{5}];

        MKTAK_samples_planes{i,1} = aplane_perk_terk{1}; MKTAK_samples_planes{i,2} = aplane_perk_terk{2}; MKTAK_samples_planes{i,3} = aplane_perk_terk{3}; 
        MKTAK_samples_planes{i,4} = aplane_perk_terk{4}; MKTAK_samples_planes{i,5} = aplane_perk_terk{5}; 
    
        plane1 = [plane1;aplane_perk_terk{1}]; % USE
        plane2 = [plane2;aplane_perk_terk{2}];
        plane3 = [plane3;aplane_perk_terk{3}];
        plane4 = [plane4;aplane_perk_terk{4}];
        plane5 = [plane5;aplane_perk_terk{5}];
        MKTAK_planes{1} = plane1; MKTAK_planes{2} = plane2; MKTAK_planes{3} = plane3; 
        MKTAK_planes{4} = plane4; MKTAK_planes{5} = plane5;

        MKTAK_xyz{i,1} = xyz_locs_terk(2:end,:);
        MKTAK_xyz{i,2} = xyz_locs_perk(2:end,:);

        date_MKTAK = exp_date;
    elseif org_group(end-1) == 'MKJNJ';
        mean_terk = mean(terk_matrix(:,3),'all'); % mean of all terk for each group
        MKJNJ_perk_terk{i} = perk_matrix(:,3) ./ mean_terk; % USE: collecting for all plane
    
        MKJNJ_byplane_perk_terk{i} = [aplane_perk_terk{1};aplane_perk_terk{2};aplane_perk_terk{3};
            aplane_perk_terk{4};aplane_perk_terk{5}];

        MKJNJ_samples_planes{i,1} = aplane_perk_terk{1}; MKJNJ_samples_planes{i,2} = aplane_perk_terk{2}; MKJNJ_samples_planes{i,3} = aplane_perk_terk{3}; 
        MKJNJ_samples_planes{i,4} = aplane_perk_terk{4}; MKJNJ_samples_planes{i,5} = aplane_perk_terk{5}; 
    
        plane1 = [plane1;aplane_perk_terk{1}]; % USE
        plane2 = [plane2;aplane_perk_terk{2}];
        plane3 = [plane3;aplane_perk_terk{3}];
        plane4 = [plane4;aplane_perk_terk{4}];
        plane5 = [plane5;aplane_perk_terk{5}];
        MKJNJ_planes{1} = plane1; MKJNJ_planes{2} = plane2; MKJNJ_planes{3} = plane3; 
        MKJNJ_planes{4} = plane4; MKJNJ_planes{5} = plane5;

        MKJNJ_xyz{i,1} = xyz_locs_terk(2:end,:);
        MKJNJ_xyz{i,2} = xyz_locs_perk(2:end,:);

        date_MKJNJ = exp_date;
    elseif org_group(end-1) == 'MK';
        mean_terk = mean(terk_matrix(:,3),'all'); % mean of all terk for each group
        MK_perk_terk{i} = perk_matrix(:,3) ./ mean_terk; % USE: collecting for all plane
    
        MK_byplane_perk_terk{i} = [aplane_perk_terk{1};aplane_perk_terk{2};aplane_perk_terk{3};
            aplane_perk_terk{4};aplane_perk_terk{5}];

        MK_samples_planes{i,1} = aplane_perk_terk{1}; MK_samples_planes{i,2} = aplane_perk_terk{2}; MK_samples_planes{i,3} = aplane_perk_terk{3}; 
        MK_samples_planes{i,4} = aplane_perk_terk{4}; MK_samples_planes{i,5} = aplane_perk_terk{5}; 
    
        plane1 = [plane1;aplane_perk_terk{1}]; % USE
        plane2 = [plane2;aplane_perk_terk{2}];
        plane3 = [plane3;aplane_perk_terk{3}];
        plane4 = [plane4;aplane_perk_terk{4}];
        plane5 = [plane5;aplane_perk_terk{5}];
        MK_planes{1} = plane1; MK_planes{2} = plane2; MK_planes{3} = plane3; 
        MK_planes{4} = plane4; MK_planes{5} = plane5;

        MK_xyz{i,1} = xyz_locs_terk(2:end,:);
        MK_xyz{i,2} = xyz_locs_perk(2:end,:);

        date_MK = exp_date;
    end
end
end

%% Identify ventral and dorsal
[int_ven_DMSO,int_dor_DMSO,iv_d_all,id_d_all] = sortvd(DMSO_xyz,DMSO_samples_planes);
[int_ven_JNJ,int_dor_JNJ,iv_j_all,id_j_all] = sortvd(JNJ_xyz,JNJ_samples_planes);
[int_ven_TAK,int_dor_TAK,iv_t_all,id_t_all] = sortvd(TAK_xyz,TAK_samples_planes);
[int_ven_MK,int_dor_MK,iv_m_all,id_m_all] = sortvd(MK_xyz,MK_samples_planes);
[int_ven_MKJNJ,int_dor_MKJNJ,iv_mj_all,id_mj_all] = sortvd(MKJNJ_xyz,MKJNJ_samples_planes);
[int_ven_MKTAK,int_dor_MKTAK,iv_mt_all,id_mt_all] = sortvd(MKTAK_xyz,MKTAK_samples_planes);

% LMM - v vs d
% [lme_vd_DMSO,md_vd_DMSO] = fixmm(int_ven_DMSO,int_dor_DMSO);
% p_value = lme_vd_DMSO.Coefficients.pValue(3);
% es_sign = sign(lme_vd_DMSO.Coefficients.Estimate(3));
% p_log_vd_DMSO = -1 * log10(p_value) * es_sign;
% 
% [lme_vd_JNJ,md_vd_JNJ] = fixmm(int_ven_JNJ,int_dor_JNJ);
% p_value = lme_vd_JNJ.Coefficients.pValue(3);
% es_sign = sign(lme_vd_JNJ.Coefficients.Estimate(3));
% p_log_vd_JNJ = -1 * log10(p_value) * es_sign;
% 
% [lme_vd_TAK,md_vd_TAK] = fixmm(int_ven_TAK,int_dor_TAK);
% p_value = lme_vd_TAK.Coefficients.pValue(3);
% es_sign = sign(lme_vd_TAK.Coefficients.Estimate(3));
% p_log_vd_TAK = -1 * log10(p_value) * es_sign;
% 
% [lme_vd_MK,md_vd_MK] = fixmm(int_ven_MK,int_dor_MK);
% p_value = lme_vd_MK.Coefficients.pValue(3);
% es_sign = sign(lme_vd_MK.Coefficients.Estimate(3));
% p_log_vd_MK = -1 * log10(p_value) * es_sign;
% 
% [lme_vd_MKJNJ,md_vd_MKJNJ] = fixmm(int_ven_MKJNJ,int_dor_MKJNJ);
% p_value = lme_vd_MKJNJ.Coefficients.pValue(3);
% es_sign = sign(lme_vd_MKJNJ.Coefficients.Estimate(3));
% p_log_vd_MKJNJ = -1 * log10(p_value) * es_sign;
% 
% [lme_vd_MKTAK,md_vd_MKTAK] = fixmm(int_ven_MKTAK,int_dor_MKTAK);
% p_value = lme_vd_MKTAK.Coefficients.pValue(3);
% es_sign = sign(lme_vd_MKTAK.Coefficients.Estimate(3));
% p_log_vd_MKTAK = -1 * log10(p_value) * es_sign;
% 
% figure('Name','Ventral vs dorsal difference at stimulated condition', 'visible','on'); clf;
% plot([p_log_vd_DMSO,p_log_vd_JNJ,p_log_vd_TAK,p_log_vd_MK,p_log_vd_MKJNJ,p_log_vd_MKTAK],'LineWidth',1.5,'Color','blue','Marker','.','MarkerSize',15); hold on;
% xticks([1,2,3,4,5,6]);
% xticklabels({'DMSO','JNJ','TAK','MK','MKJNJ','MKTAK'});
% xlim([0.5 6.5]);
% ylim([0 24]);
% ylabel('-log_{10}(p)','Interpreter', 'tex')
% ax = gca;
% ax.XAxis.FontSize = 13;
% ax.YAxis.FontSize = 13;
% ax.XAxis.LineWidth = 1.5;
% ax.YAxis.LineWidth = 1.5;
% ax.Box = "off";

%% ven vs dorsal within group
[h_d,p_d,~,stats_d] = ttest2(iv_d_all,id_d_all); 
[h_j,p_j,~,stats_j] = ttest2(iv_j_all,id_d_all);
[h_t,p_t,~,stats_t] = ttest2(iv_t_all,id_t_all);
[h_m,p_m,~,stats_m] = ttest2(iv_m_all,id_m_all);
[h_mj,p_mj,~,stats_mj] = ttest2(iv_mj_all,id_mj_all);
[h_mt,p_mt,~,stats_mt] = ttest2(iv_mt_all,id_mt_all);

% [f,bc] = boxplot2(iv_d_all,id_d_all, 'Ventral','Dorsal','emptyname'); hold on;
% ylim([0 4])
% plot([1.1,1.9],[3.5 3.5],'LineWidth',1.5,'Color',[0 0 0]);
% text(1.35,3.6,'***','FontSize',24); hold off;
% [f,bc] = boxplot2(iv_j_all,id_d_all, 'Ventral','Dorsal','emptyname'); hold on;
% ylim([0 4])
% plot([1.1,1.9],[3.5 3.5],'LineWidth',1.5,'Color',[0 0 0]);
% text(1.35,3.6,'***','FontSize',24); hold off;
% [f,bc] = boxplot2(iv_t_all,id_t_all, 'Ventral','Dorsal','emptyname'); hold on;
% ylim([0 4])
% plot([1.1,1.9],[3.5 3.5],'LineWidth',1.5,'Color',[0 0 0]);
% text(1.35,3.6,'***','FontSize',24); hold off;
% [f,bc] = boxplot2(iv_m_all,id_m_all, 'Ventral','Dorsal','emptyname'); hold on;
% ylim([0 4])
% plot([1.1,1.9],[3.5 3.5],'LineWidth',1.5,'Color',[0 0 0]);
% text(1.35,3.6,'***','FontSize',24); hold off;
% [f,bc] = boxplot2(iv_mj_all,id_mj_all, 'Ventral','Dorsal','emptyname'); hold on;
% ylim([0 4])
% plot([1.1,1.9],[3.5 3.5],'LineWidth',1.5,'Color',[0 0 0]);
% text(1.4,3.7,'ns','FontSize',17); hold off;
% [f,bc] = boxplot2(iv_mt_all,id_mt_all, 'Ventral','Dorsal','emptyname'); hold on;
% ylim([0 4])
% plot([1.1,1.9],[3.5 3.5],'LineWidth',1.5,'Color',[0 0 0]);
% text(1.35,3.6,'***','FontSize',24); hold off;

%% subnuclei difference
[dmso_mean,dmso_sem,dmso_meandiff] = subnuclei_diff(iv_d_all,id_d_all);
[jnj_mean,jnj_sem,jnj_meandiff] = subnuclei_diff(iv_j_all,id_j_all);
[tak_mean,tak_sem,tak_meandiff] = subnuclei_diff(iv_t_all,id_t_all);
[mk_mean,mk_sem,mk_meandiff] = subnuclei_diff(iv_m_all,id_m_all);
[mkjnj_mean,mkjnj_sem,mkjnj_meandiff] = subnuclei_diff(iv_mj_all,id_mj_all);
[mktak_mean,mktak_sem,mktak_meandiff] = subnuclei_diff(iv_mt_all,id_mt_all);

% anova
% DMSO v JNJ v TAK
test_data = [dmso_meandiff;jnj_meandiff;tak_meandiff];
test_group = [ones(1,length(dmso_meandiff))*1, ones(1,length(jnj_meandiff))*2, ones(1,length(tak_meandiff))*3];
[p_dmso,tbl_dmso,stats_dmso] = anova1(test_data',test_group,'off');
p_multcomp_dmso = multcompare(stats_dmso,"CriticalValueType","bonferroni");

% MK v MKJNJ v MKTAK
test_data = [mk_meandiff;mkjnj_meandiff;mktak_meandiff];
test_group = [ones(1,length(mk_meandiff))*1, ones(1,length(mkjnj_meandiff))*2, ones(1,length(mktak_meandiff))*3];
[p_mk,tbl_mk,stats_mk] = anova1(test_data',test_group,'off');
p_multcomp_mk = multcompare(stats_mk,"CriticalValueType","bonferroni");

% DMSO v MK
[hyp_vd,pro_vd,ci_vd,stats_vd] = ttest2(dmso_meandiff,mk_meandiff,'Vartype','unequal');

% graph for subnuclei difference and v-d within group
f_d = figure('Name','Ventral and Dorsal Intensities Difference by Group DMSOJNJTAK','Visible','off');clf; hold on;
y = [dmso_mean,jnj_mean,tak_mean];
err_bar = [dmso_sem,jnj_sem,tak_sem];
errorbar(1:3,y,err_bar,'LineWidth',1.5,'Color','blue');
ylabel('Subnuclei pERK/tERK Difference (vHb - dHb)','FontSize', 20)
xlim([0.5 3.5])
ylim([0 0.6])
xticks([1,2,3])
xticklabels({'DMSO','JNJ','TAK'})
ax = gca;
ax.FontName = 'Arial';
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
ax.XAxis.LineWidth = 1.5;
ax.YAxis.LineWidth = 1.5;
ax.Box = "off";
plot([1.05,1.95],[0.51 0.51],'LineWidth',1.5,'Color',[0 0 0]);
text(1.4,0.53,'ns','FontSize',17); 
plot([1.05,2.95],[0.55 0.55],'LineWidth',1.5,'Color',[0 0 0]);
text(1.9,0.57,'ns','FontSize',17);
plot([2.05,2.95],[0.51 0.51],'LineWidth',1.5,'Color',[0 0 0]);
text(2.4,0.52,'*','FontSize',24); hold off;

f_m = figure('Name','Ventral and Dorsal Intensities Difference by Group MKJNJTAK','Visible','off');clf; hold on;
y = [mk_mean,mkjnj_mean,mktak_mean];
err_bar = [mk_sem,mkjnj_sem,mktak_sem];
errorbar(1:3,y,err_bar,'LineWidth',1.5,'Color','blue');
ylabel('Subnuclei pERK/tERK Difference (vHb - dHb)','FontSize', 20)
xlim([0.5 3.5])
ylim([0 0.8])
xticks([1,2,3])
xticklabels({'MK','MKJNJ','MKTAK'})
ax = gca;
ax.FontName = 'Arial';
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
ax.XAxis.LineWidth = 1.5;
ax.YAxis.LineWidth = 1.5;
ax.Box = "off";
plot([1.05,1.95],[0.7 0.7],'LineWidth',1.5,'Color',[0 0 0]);
text(1.35,0.71,'***','FontSize',24); 
plot([1.05,2.95],[0.75 0.75],'LineWidth',1.5,'Color',[0 0 0]);
text(1.9,0.77,'ns','FontSize',17);
plot([2.05,2.95],[0.7 0.7],'LineWidth',1.5,'Color',[0 0 0]);
text(2.40,0.72,'ns','FontSize',17); hold off;

%% stats for between group and graphs
% ven vs ven
[h_vj,p_vj,~,stats_vj] = ttest2(iv_j_all,iv_d_all);
[h_vt,p_vt,~,stats_vt] = ttest2(iv_t_all,iv_d_all);
[h_vm,p_vm,~,stats_vm] = ttest2(iv_m_all,iv_d_all);
[h_vmj,p_vmj,~,stats_vmj] = ttest2(iv_mj_all,iv_m_all);
[h_vmt,p_vmt,~,stats_vmt] = ttest2(iv_mt_all,iv_m_all);

% anova ven v ven
test_data = [iv_d_all;iv_j_all;iv_t_all];
test_group = [ones(1,length(iv_d_all))*1, ones(1,length(iv_j_all))*2, ones(1,length(iv_t_all))*3];
[p_djt,tbl_djt,stats_djt] = anova1(test_data,test_group','off');
p_djt_bon = multcompare(stats_djt,"CriticalValueType","bonferroni");

test_data = [iv_m_all;iv_mj_all;iv_mt_all];
test_group = [ones(1,length(iv_m_all))*1, ones(1,length(iv_mj_all))*2, ones(1,length(iv_mt_all))*3];
[p_mjt,tbl_mjt,stats_mjt] = anova1(test_data,test_group','off');
p_mjt_bon = multcompare(stats_mjt,"CriticalValueType","bonferroni");

% dor v dor
test_data = [id_d_all;id_j_all;id_t_all];
test_group = [ones(1,length(id_d_all))*1, ones(1,length(id_j_all))*2, ones(1,length(id_t_all))*3];
[p_dor_djt,tbl_dor_djt,stats_dor_djt] = anova1(test_data,test_group','off');
p_dor_djt_bon = multcompare(stats_dor_djt,"CriticalValueType","bonferroni");

test_data = [id_m_all;id_mj_all;id_mt_all];
test_group = [ones(1,length(id_m_all))*1, ones(1,length(id_mj_all))*2, ones(1,length(id_mt_all))*3];
[p_dor_mjt,tbl_dor_mjt,stats_dor_mjt] = anova1(test_data,test_group','off');
p_dor_mjt_bon = multcompare(stats_dor_mjt,"CriticalValueType","bonferroni");

return

% two-way anova (CANNOT because of perfect multicolinearity)
% % DMSO JNJ TAK MK
% % ventral
% test_data_DMSO = [iv_d_all;iv_j_all;iv_t_all;iv_m_all];
% test_group_DMSO = [ones(1,length(iv_d_all))*1, ones(1,length(iv_j_all))*2, ones(1,length(iv_t_all))*3, ones(1,length(iv_m_all))*4];
% obj1_leng = sum([double(length(iv_d_all)),double(length(iv_j_all)),double(length(iv_t_all))]);
% test_factor_DMSO(1:length(iv_d_all)) = "normal";
% test_factor_DMSO(length(iv_d_all)+1:obj1_leng) = "enhance";
% test_factor_DMSO(obj1_leng+1:obj1_leng+length(iv_m_all)) = "impair";
% [~,~,stats_DMSO] = anovan(test_data_DMSO,{test_group_DMSO test_factor_DMSO},"Model",2, ...
%     "Varnames",["test_group_DMSO", "test_factor_DMSO"]);
% [results,~,~,gnames] = multcompare(stats_DMSO,"Dimension",[1 2]);


% ven
% [f_v_d,bc_v_d] = boxplot3(iv_d_all,iv_j_all,iv_t_all,'DMSO','JNJ','TAK','pERK/tERK ven DMSO'); hold on;
% ylim([0 4])
% plot([1.1,1.9],[3.5 3.5],'LineWidth',1.5,'Color',[0 0 0]);
% text(1.3,3.55,'***','FontSize',24,'FontWeight','bold');
% plot([1.1,2.9],[3.8 3.8],'LineWidth',1.5,'Color',[0 0 0]);
% text(1.8,3.85,'***','FontSize',24,'FontWeight','bold');
% plot([2.1,2.9],[3.5 3.5],'LineWidth',1.5,'Color',[0 0 0]);
% text(2.35,3.55,'**','FontSize',24,'FontWeight','bold');
% 
% [f_v_m,bc_v_m] = boxplot3(iv_m_all,iv_mj_all,iv_mt_all,'MK','MKJNJ','MKTAK','pERK/tERK ven MK'); hold on;
% ylim([0 4])
% plot([1.1,1.9],[3.5 3.5],'LineWidth',1.5,'Color',[0 0 0]);
% text(1.3,3.55,'***','FontSize',24,'FontWeight','bold');
% plot([1.1,2.9],[3.8 3.8],'LineWidth',1.5,'Color',[0 0 0]);
% text(1.8,3.85,'***','FontSize',24,'FontWeight','bold');
% plot([2.1,2.9],[3.5 3.5],'LineWidth',1.5,'Color',[0 0 0]);
% text(2.35,3.6,'ns','FontSize',17,'FontWeight','bold');

% [f_v_control,bc_v_control] = boxplot2(iv_d_all,iv_m_all,'DMSO','MK','pERK/tERK dor MK'); hold on;
% ylim([0 4])
% plot([1.1,1.9],[3.5 3.5],'LineWidth',1.5,'Color',[0 0 0]);
% text(1.3,3.55,'***','FontSize',24,'FontWeight','bold');

% dor vs dor
[h_dj,p_dj,~,stats_dj] = ttest2(id_j_all,id_d_all);
[h_dt,p_dt,~,stats_dt] = ttest2(id_t_all,id_d_all);
[h_dm,p_dm,~,stats_dm] = ttest2(id_m_all,id_d_all);
[h_dmj,p_dmj,~,stats_dmj] = ttest2(id_mj_all,id_m_all);
[h_dmt,p_dmt,~,stats_dmt] = ttest2(id_mt_all,id_m_all);

% [f_d_d,bc_d_d] = boxplot3(id_d_all,id_j_all,id_t_all,'DMSO','JNJ','TAK','pERK/tERK dor DMSO'); hold on;
% ylim([0 4])
% plot([1.1,1.9],[3.5 3.5],'LineWidth',1.5,'Color',[0 0 0]);
% text(1.3,3.55,'***','FontSize',24,'FontWeight','bold');
% plot([1.1,2.9],[3.8 3.8],'LineWidth',1.5,'Color',[0 0 0]);
% text(1.8,3.85,'***','FontSize',24,'FontWeight','bold');
% plot([2.1,2.9],[3.5 3.5],'LineWidth',1.5,'Color',[0 0 0]);
% text(2.35,3.6,'ns','FontSize',17,'FontWeight','bold');
% 
% [f_d_m,bc_d_m] = boxplot3(id_m_all,id_mj_all,id_mt_all,'MK','MKJNJ','MKTAK','pERK/tERK dor MK'); hold on;
% ylim([0 4])
% plot([1.1,1.9],[3.5 3.5],'LineWidth',1.5,'Color',[0 0 0]);
% text(1.3,3.55,'***','FontSize',24,'FontWeight','bold');
% plot([1.1,2.9],[3.8 3.8],'LineWidth',1.5,'Color',[0 0 0]);
% text(1.8,3.85,'***','FontSize',24,'FontWeight','bold');
% plot([2.1,2.9],[3.5 3.5],'LineWidth',1.5,'Color',[0 0 0]);
% text(2.35,3.55,'**','FontSize',24,'FontWeight','bold');

% [f_d_control,bc_d_control] = boxplot2(id_d_all,id_m_all,'DMSO','MK','pERK/tERK dor MK'); hold on;
% ylim([0 4])
% plot([1.1,1.9],[3.5 3.5],'LineWidth',1.5,'Color',[0 0 0]);
% text(1.3,3.55,'***','FontSize',24,'FontWeight','bold');



% %% Linear mixed model - ven vs dor
% % ventral
% [lme_d_JNJ,md_v_JNJ] = fixmm(int_ven_JNJ,int_ven_DMSO);
% p_value = lme_d_JNJ.Coefficients.pValue(3);
% es_sign = sign(lme_d_JNJ.Coefficients.Estimate(3));
% p_log_v_JNJ = -1 * log10(p_value) * es_sign;
% 
% [lme_d_TAK,md_v_TAK] = fixmm(int_ven_TAK,int_ven_DMSO);
% p_value = lme_d_TAK.Coefficients.pValue(3);
% es_sign = sign(lme_d_TAK.Coefficients.Estimate(3));
% p_log_v_TAK = -1 * log10(p_value) * es_sign;
% 
% [lme_vd_MK,md_v_MK] = fixmm(int_ven_MK,int_ven_DMSO);
% p_value = lme_vd_MK.Coefficients.pValue(3);
% es_sign = sign(lme_vd_MK.Coefficients.Estimate(3));
% p_log_v_MK = -1 * log10(p_value) * es_sign;
% 
% [lme_vd_MKJNJ,md_v_MKJNJ] = fixmm(int_ven_MKJNJ,int_ven_MK);
% p_value = lme_vd_MKJNJ.Coefficients.pValue(3);
% es_sign = sign(lme_vd_MKJNJ.Coefficients.Estimate(3));
% p_log_v_MKJNJ = -1 * log10(p_value) * es_sign;
% 
% [lme_vd_MKTAK,md_v_MKTAK] = fixmm(int_ven_MKTAK,int_ven_MK);
% p_value = lme_vd_MKTAK.Coefficients.pValue(3);
% es_sign = sign(lme_vd_MKTAK.Coefficients.Estimate(3));
% p_log_v_MKTAK = -1 * log10(p_value) * es_sign;
% 
% % dorsal
% [lme_d_JNJ,md_d_JNJ] = fixmm(int_dor_JNJ,int_dor_DMSO);
% p_value = lme_d_JNJ.Coefficients.pValue(3);
% es_sign = sign(lme_d_JNJ.Coefficients.Estimate(3));
% p_log_d_JNJ = -1 * log10(p_value) * es_sign;
% 
% [lme_d_TAK,md_d_TAK] = fixmm(int_dor_TAK,int_dor_DMSO);
% p_value = lme_d_TAK.Coefficients.pValue(3);
% es_sign = sign(lme_d_TAK.Coefficients.Estimate(3));
% p_log_d_TAK = -1 * log10(p_value) * es_sign;
% 
% [lme_vd_MK,md_d_MK] = fixmm(int_dor_MK,int_dor_DMSO);
% p_value = lme_vd_MK.Coefficients.pValue(3);
% es_sign = sign(lme_vd_MK.Coefficients.Estimate(3));
% p_log_d_MK = -1 * log10(p_value) * es_sign;
% 
% [lme_vd_MKJNJ,md_d_MKJNJ] = fixmm(int_dor_MKJNJ,int_dor_MK);
% p_value = lme_vd_MKJNJ.Coefficients.pValue(3);
% es_sign = sign(lme_vd_MKJNJ.Coefficients.Estimate(3));
% p_log_d_MKJNJ = -1 * log10(p_value) * es_sign;
% 
% [lme_d_MKTAK,md_d_MKTAK] = fixmm(int_dor_MKTAK,int_dor_MK);
% p_value = lme_d_MKTAK.Coefficients.pValue(3);
% es_sign = sign(lme_d_MKTAK.Coefficients.Estimate(3));
% p_log_d_MKTAK = -1 * log10(p_value) * es_sign;
% 
% cg_mat = [...
%     p_log_v_JNJ,p_log_d_JNJ;...
%     p_log_v_TAK,p_log_d_TAK;...
%     p_log_v_MK,p_log_d_MK;...
%     p_log_v_MKJNJ,p_log_d_MKJNJ;...
%     p_log_v_MKTAK,p_log_d_MKTAK;];
% 
% cg = clustergram(cg_mat);
% cg.ColumnLabels = {'ventral','dorsal'};
% cg.RowLabels = {'JNJ','TAK','MK','MKJNJ','MKTAK'};
% cg.DisplayRange = 3;
% cg.ShowDendrogram = 'off';
% 
% %% Linear mixed model -- habenula as a whole
% % DMSO_pt_all = [DMSO_perk_terk{1};DMSO_perk_terk{2};DMSO_perk_terk{3}];
% % MK_pt_all = [MK_perk_terk{1};MK_perk_terk{2};MK_perk_terk{3}];
% % MKJNJ_pt_all = [MKJNJ_perk_terk{1};MKJNJ_perk_terk{2};MKJNJ_perk_terk{3}];
% % MKTAK_pt_all = [MKTAK_perk_terk{1};MKTAK_perk_terk{2};MKTAK_perk_terk{3}];
% % TAK_pt_all = [TAK_perk_terk{1};TAK_perk_terk{2};TAK_perk_terk{3}];
% % JNJ_pt_all = [JNJ_perk_terk{1};JNJ_perk_terk{2};JNJ_perk_terk{3}];
% % 
% % beta_MK = mean(MK_pt_all,'all') - mean(DMSO_pt_all,'all');
% % beta_MKTAK = mean(MKJNJ_pt_all,'all') - mean(DMSO_pt_all,'all');
% % beta_MKJNJ = mean(MKTAK_pt_all,'all') - mean(DMSO_pt_all,'all');
% % beta_TAK = mean(JNJ_pt_all,'all') - mean(DMSO_pt_all,'all');
% % beta_JNJ = mean(TAK_pt_all,'all') - mean(DMSO_pt_all,'all');
% % 
% % neu_int = [DMSO_pt_all;MK_pt_all;MKJNJ_pt_all;MKTAK_pt_all;JNJ_pt_all;TAK_pt_all];
% % grp_varia = [ones(1,length(DMSO_pt_all))*1,ones(1,length(MK_pt_all))*2,ones(1,length(MKJNJ_pt_all))*3,...
% %     ones(1,length(MKTAK_pt_all))*4,ones(1,length(TAK_pt_all))*5,ones(1,length(JNJ_pt_all))*6];
% % 
% % grp_beta = [ones(1,length(DMSO_pt_all))*1,ones(1,length(MK_pt_all))*beta_MK,ones(1,length(MKJNJ_pt_all))*beta_MKJNJ,...
% %     ones(1,length(MKTAK_pt_all))*beta_MKTAK,ones(1,length(TAK_pt_all))*beta_TAK,ones(1,length(JNJ_pt_all))*beta_JNJ];
% % 
% % grp_date = [ones(1,length(DMSO_pt_all))*1,ones(1,length(MK_pt_all))*2,ones(1,length(MKJNJ_pt_all))*3,...
% %     ones(1,length(MKTAK_pt_all))*3,ones(1,length(TAK_pt_all))*4,ones(1,length(JNJ_pt_all))*1];
% 
% [lme_JNJ,md_JNJ] = fixmm(JNJ_samples_planes,DMSO_samples_planes);
% p_value = lme_JNJ.Coefficients.pValue(3);
% es_sign = sign(lme_JNJ.Coefficients.Estimate(3));
% p_log_JNJ = -1 * log10(p_value) * es_sign;
% 
% [lme_TAK,md_TAK] = fixmm(TAK_samples_planes,DMSO_samples_planes);
% p_value = lme_TAK.Coefficients.pValue(3);
% es_sign = sign(lme_TAK.Coefficients.Estimate(3));
% p_log_TAK = -1 * log10(p_value) * es_sign;
% 
% [lme_MK,md_MK] = fixmm(MK_samples_planes,DMSO_samples_planes);
% p_value = lme_MK.Coefficients.pValue(3);
% es_sign = sign(lme_MK.Coefficients.Estimate(3));
% p_log_MK = -1 * log10(p_value) * es_sign;
% 
% [lme_MKJNJ,md_MKJNJ] = fixmm(MKJNJ_samples_planes,MK_samples_planes);
% p_value = lme_MKJNJ.Coefficients.pValue(3);
% es_sign = sign(lme_MKJNJ.Coefficients.Estimate(3));
% p_log_MKJNJ = -1 * log10(p_value) * es_sign;
% 
% [lme_MKTAK,md_MKTAK] = fixmm(MKTAK_samples_planes,MK_samples_planes);
% p_value = lme_MKTAK.Coefficients.pValue(3);
% es_sign = sign(lme_MKTAK.Coefficients.Estimate(3));
% p_log_MKTAK = -1 * log10(p_value) * es_sign;
% 
% pheatmap_matrix = [p_log_JNJ,p_log_TAK,p_log_MK,p_log_MKJNJ,p_log_MKTAK];
% % clustergram(pheatmap_matrix_DMSO)

%% all funcs

function [int_ven,int_dor,int_ven_vect,int_dor_vect] = sortvd(locs,int)
    sample_size = size(locs,2);
    int_p2_ven = []; int_p3_ven = []; int_p4_ven = [];
    int_p2_dor = []; int_p3_dor = []; int_p4_dor = [];
    int_ven_vect = []; int_dor_vect = [];
    for k = 1:sample_size
        perk_locs = locs{k,2};
        for h = 2:4
        grp_locs = find(perk_locs(:,3) == h); % plane 2
        x_locs = perk_locs(grp_locs,1); % only x_locs
        ven_x = find(x_locs < 1 | x_locs > 6);
        dor_x = find((x_locs > 2 & x_locs < 3) | (x_locs > 4 & x_locs < 5));
            if h == 2
                int_plane = int{k,h};
                int_p2_ven = [int_p2_ven; int_plane(ven_x)]; % separate plane, all samples
                int_p2_dor = [int_p2_dor; int_plane(dor_x)];
            elseif h == 3
                int_plane = int{k,h};
                int_p3_ven = [int_p3_ven; int_plane(ven_x)];
                int_p3_dor = [int_p3_dor; int_plane(dor_x)];
            elseif h == 4
                int_plane = int{k,h};
                int_p4_ven = [int_p4_ven; int_plane(ven_x)];
                int_p4_dor = [int_p4_dor; int_plane(dor_x)];
            end
        end
        int_ven{k,2} = int_p2_ven; int_ven{k,3} = int_p3_ven; int_ven{k,4} = int_p4_ven;
        int_dor{k,2} = int_p2_dor; int_dor{k,3} = int_p3_dor; int_dor{k,4} = int_p4_dor;
        int_ven_vect = [int_ven_vect; int_p2_ven;int_p3_ven;int_p4_ven];
        int_dor_vect = [int_dor_vect; int_p2_dor;int_p3_dor;int_p4_dor];
    end
end

function [lme,mean_diff,m_t,m_c] = fixmm(test_group,control_group)
    sample_size_test = size(test_group,1);
    test_planes = []; vol_test = []; sub_test = [];
    for i = 1:sample_size_test
        test_planes = [test_planes; test_group{i,2};test_group{i,3};test_group{i,4}];
        vol_test = [vol_test, ones(1,length(test_group{i,2}))*1,ones(1,length(test_group{i,3}))*2,ones(1,length(test_group{i,4}))*3];
        sub_test = [sub_test, ones(1,length(test_group{i,2}))*i,ones(1,length(test_group{i,3}))*i,ones(1,length(test_group{i,4}))*i];
    end
    sample_size_control = size(control_group,1);
    control_planes = []; vol_control = []; sub_control = [];
    for j = 1:sample_size_control
        control_planes = [control_planes; control_group{j,2};control_group{j,3};control_group{j,4}];
        vol_control = [vol_control, ones(1,length(control_group{j,2}))*1,ones(1,length(control_group{j,3}))*2,ones(1,length(control_group{j,4}))*3];
        sub_control = [sub_control, ones(1,length(control_group{j,2}))*j,ones(1,length(control_group{j,3}))*j,ones(1,length(control_group{j,4}))*j];
    end

    rept_sub = [sub_control,sub_test];
    date = [ones(1, length(control_planes))*1,ones(1, length(test_planes))*2];
    grp_varia = date;
    beta = ones(1,length(date)) * (mean(control_planes,'all') - mean(test_planes,'all'));
    mean_diff = mean(control_planes,'all') - mean(test_planes,'all');
    m_c = mean(control_planes,'all');
    m_t = mean(test_planes,'all');

    neu_int = [control_planes;test_planes];
    vol_tc = [vol_control,vol_test];

    neu_int_table = table(neu_int,date',beta',vol_tc',grp_varia',rept_sub','VariableNames',{'Intensity', 'GroupDate','GroupBeta','Volume','GroupName','RepeatedSub'});
    lme = fitlme(neu_int_table, 'Intensity ~ GroupName + Volume + (1|RepeatedSub)');
end

% four input
function [f,bc] = boxplot2(input1,input2, input1_name,input2_name,fig_name)
    uniformColor = 'blue';
    empty_identifiers = [ones(length(input1),1)*1; ones(length(input2),1)*2];
    identifiers_name = {input1_name,input2_name}; % CHANGE name
    identifiers = categorical(identifiers_name(empty_identifiers));
    amp_matrix_MK = [input1;input2];
    f = figure('Name',fig_name,'Visible','on');clf
    bc = boxchart(identifiers,amp_matrix_MK,'Notch','off','JitterOutliers','on','MarkerStyle','.'); hold on;
    bc.BoxFaceColor = uniformColor;
    bc.MarkerColor = uniformColor;
    ylabel('pERK/tERK value', 'FontSize', 20);
    ax = gca;
    ax.XAxis.Categories = categorical(identifiers_name);
    ax.FontName = 'Arial';
    ax.XAxis.FontSize = 15;
    ax.YAxis.FontSize = 16;
    ax.XAxis.LineWidth = 1.5;
    ax.YAxis.LineWidth = 1.5;
    ax.Box = "off";
    ax.TickLength = [0 0];
    pbaspect([0.5 0.7 0.7]);
end

%three inputs
function [f,bc] = boxplot3(input1,input2,input3, input1_name,input2_name,input3_name,fig_name)
    uniformColor = 'blue';
    empty_identifiers = [ones(length(input1),1)*1; ones(length(input2),1)*2;ones(length(input3),1)*3];
    identifiers_name = {input1_name,input2_name,input3_name}; % CHANGE name
    identifiers = categorical(identifiers_name(empty_identifiers));
    amp_matrix_MK = [input1;input2;input3];
    f = figure('Name',fig_name,'Visible','on');clf
    bc = boxchart(identifiers,amp_matrix_MK,'Notch','off','JitterOutliers','on','MarkerStyle','.'); hold on;
    bc.BoxFaceColor = uniformColor;
    bc.MarkerColor = uniformColor;
    ylabel('pERK/tERK value', 'FontSize', 20);
    ax = gca;
    ax.XAxis.Categories = categorical(identifiers_name);
    ax.FontName = 'Arial';
    ax.XAxis.FontSize = 15;
    ax.YAxis.FontSize = 16;
    ax.XAxis.LineWidth = 1.5;
    ax.YAxis.LineWidth = 1.5;
    ax.Box = "off";
    ax.TickLength = [0 0];
    pbaspect([0.5 0.7 0.7]);
end

function [mean_vect,sem_vect] = mean_sem(meandiff)
mean_vect = mean(meandiff,'all'); 
sd_vect = std(meandiff);
n = length(meandiff);
sem_vect = sd_vect / sqrt(n);
end

function [mean_vect,sem_vect,dmso_meandiff] = subnuclei_diff(iv_d_all,id_d_all)
v_indices = randperm(length(iv_d_all), 90);
d_indices = randperm(length(id_d_all), 90);
dmso_meandiff = iv_d_all(v_indices) - id_d_all(d_indices);
[mean_vect,sem_vect] = mean_sem(dmso_meandiff);
end