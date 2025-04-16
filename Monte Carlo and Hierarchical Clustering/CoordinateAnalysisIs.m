lfpr_cell = [];
hfpr_cell = [];
All_Is_Data_No_Empty = [];
lengths_lfprcell = [];
for nnr = 1:10
    if nnr == 6
    nnr = nnr + 1;
    else
        All_Is_Data_No_Empty{nnr} = All_Is_Data{nnr};
    end
end

for ll = 1:9
    lfpr_cell{ll} = All_Is_Data_No_Empty{ll}(:,3);
    hfpr_cell{ll} = All_Is_Data_No_Empty{ll}(:,4)
    lengths_lfprcell{ll} = length(lfpr_cell{ll});
end

total_az = cell2mat(lengths_lfprcell);
sum_total_az = sum(total_az);
for nnnn = 1:numel(lfpr_cell)
    if nnnn == 1
    total_lfpr = cell2mat(lfpr_cell(nnnn));
    else
    templfpr = cell2mat(lfpr_cell(nnnn))
    total_lfpr = vertcat(total_lfpr,templfpr)
    end
end
hist_Exp = expfit(total_lfpr);
%%%%%Make a
for nnnn = 1:numel(lfpr_cell)
    if nnnn == 1
    total_lfpr = cell2mat(lfpr_cell(nnnn));
    else
    templfpr = cell2mat(lfpr_cell(nnnn))
    total_lfpr = vertcat(total_lfpr,templfpr)
    end
end
low_lf_thresh_upper = quantile(total_lfpr, 0.33);
med_lf_thresh_upper = quantile(total_lfpr, 0.66);
%%%For me 2/9/2024 Make a low med and high pr lf_pr cell array with
%%%coordinates for each NMJ. Run g function. Test for CSR.
data_coord_pr = [];
for ii = 1:9
    data_coord_pr{ii} = [All_Is_Data_No_Empty{ii}(:,1) All_Is_Data_No_Empty{ii}(:,2) All_Is_Data_No_Empty{ii}(:,3) All_Is_Data_No_Empty{ii}(:,4)] ;
end
hf_low_pr_coord = [];
hf_med_pr_coord = [];
hf_high_pr_coord = [];

for ii = 1:9
  low_ind = find(All_Is_Data_No_Empty{ii}(:,4)<=low_lf_thresh_upper);
  med_ind = find(All_Is_Data_No_Empty{ii}(:,4)>low_lf_thresh_upper & All_Is_Data_No_Empty{ii}(:,4)<=med_lf_thresh_upper);
  high_ind = find(All_Is_Data_No_Empty{ii}(:,4)>med_lf_thresh_upper);
  hf_low_pr_coord{ii} = [All_Is_Data_No_Empty{ii}(low_ind,1) All_Is_Data_No_Empty{ii}(low_ind,2) All_Is_Data_No_Empty{ii}(low_ind,3) All_Is_Data_No_Empty{ii}(low_ind,4)];
  hf_med_pr_coord{ii} = [All_Is_Data_No_Empty{ii}(med_ind,1) All_Is_Data_No_Empty{ii}(med_ind,2) All_Is_Data_No_Empty{ii}(med_ind,3) All_Is_Data_No_Empty{ii}(med_ind,4)];
  hf_high_pr_coord{ii} = [All_Is_Data_No_Empty{ii}(high_ind,1) All_Is_Data_No_Empty{ii}(high_ind,2) All_Is_Data_No_Empty{ii}(high_ind,3) All_Is_Data_No_Empty{ii}(high_ind,4)];
end
%%Calculate G function for each group for each NMJ
hGfunctionaAll = []
for oooo = 1:9
    LowPoint = [hf_low_pr_coord{oooo}(:,1) hf_low_pr_coord{oooo}(:,2)];
    MedPoint = [hf_med_pr_coord{oooo}(:,1) hf_med_pr_coord{oooo}(:,2)];
    HighPoint = [hf_high_pr_coord{oooo}(:,1) hf_high_pr_coord{oooo}(:,2)];
    maxD_L = computeMaxDistance(LowPoint);
    maxD_M = computeMaxDistance(MedPoint);
    maxD_H = computeMaxDistance(HighPoint);
    hfrRange_L = 0:1:maxD_L;
    hfrRange_M = 0:1:maxD_M;
    hfrRange_H = 0:1:maxD_H;
    hfgFunct_L = computeGFunction(LowPoint, hfrRange_L);
    hfgFunct_M = computeGFunction(MedPoint, hfrRange_M);
    hfgFunct_H = computeGFunction(HighPoint, hfrRange_H);
    hGfunctionaAll{oooo}.GF_L = hfgFunct_L;
    hGfunctionaAll{oooo}.GF_M = hfgFunct_M;
    hGfunctionaAll{oooo}.GF_H = hfgFunct_H;
    hGfunctionaAll{oooo}.Range_L = hfrRange_L;
    hGfunctionaAll{oooo}.Range_M = hfrRange_M;
    hGfunctionaAll{oooo}.Range_H = hfrRange_H;
end
%Plot Low Pr for all NMJ
% figure;
% for lll = 1:10
% 
% plot(GfunctionaAll{lll}.Range_L, GfunctionaAll{lll}.GF_L, 'LineWidth', 2);
% xlabel('Distance (r)');
% ylabel('G(r)');
% title('Nearest Neighbors G Function');
% grid on;
% hold on
% 
% end
hflength_mat = double.empty;
hflength_matM = double.empty;
hflength_matH = double.empty;
for jjj = 1:9
length_mat_temp1 = length(hGfunctionaAll{jjj}.GF_L);
hflength_mat = [hflength_mat length_mat_temp1];
end
for jjj = 1:9
length_mat_temp2 = length(hGfunctionaAll{jjj}.GF_M);
hflength_matM = [hflength_matM length_mat_temp2];
end
for jjj = 1:9
length_mat_temp3 = length(hGfunctionaAll{jjj}.GF_H);
hflength_matH = [hflength_matH length_mat_temp3];
end
hfrRange_L_Min = min(hflength_mat);
hfrRange_M_Min = min(hflength_matM);
hfrRange_H_Min = min(hflength_matH);

hfGfunctionLAll_Minimized = zeros(9,hfrRange_L_Min);
hfGfunctionMAll_Minimized = zeros(9,hfrRange_M_Min);
hfGfunctionHAll_Minimized = zeros(9,hfrRange_H_Min);


for rrr = 1:9
hfGfunctionLAll_Minimized(rrr,1:hfrRange_L_Min) = hGfunctionaAll{rrr}.GF_L(1,1:hfrRange_L_Min);
hfGfunctionMAll_Minimized(rrr,1:hfrRange_M_Min) = hGfunctionaAll{rrr}.GF_M(1,1:hfrRange_M_Min);
hfGfunctionHAll_Minimized(rrr,1:hfrRange_H_Min) = hGfunctionaAll{rrr}.GF_H(1,1:hfrRange_H_Min);
end

figure(1)
stdshade(hfGfunctionLAll_Minimized,0.3);
xlabel('Distance (r)');
ylabel('G(r)');
title('Nearest Neighbors G Function for Low Pr Sites');
grid on;

figure(2)
stdshade(hfGfunctionMAll_Minimized,0.3);
xlabel('Distance (r)');
ylabel('G(r)');
title('Nearest Neighbors G Function for Medium Pr Sites');
grid on;

figure(3)
stdshade(hfGfunctionHAll_Minimized,0.3);
xlabel('Distance (r)');
ylabel('G(r)');
title('Nearest Neighbors G Function for High Pr Sites');
grid on;
% Generate random points for demonstration

% Generate a random number from an exponential distribution
randomNumber = exprnd(hist_Exp);
SimuPr = zeros(1,400);
for iii = 1:400
randomNumber = exprnd(hist_Exp);
SimuPr(iii) = randomNumber;
end

hfsim_data_coord_pr = [];
for ii = 1:9
    sim_x_coord = zeros(length(All_Is_Data_No_Empty{ii}(:,1)),1);
    sim_y_coord = zeros(length(All_Is_Data_No_Empty{ii}(:,1)),1);
    sim_hf_pr = zeros(length(All_Is_Data_No_Empty{ii}(:,1)),1);
    for nnnn = 1:length(All_Is_Data_No_Empty{ii}(:,1))
        sim_x_coord(nnnn,1) = All_Is_Data_No_Empty{ii}(nnnn,1);
        sim_y_coord(nnnn,1) = All_Is_Data_No_Empty{ii}(nnnn,2);
        sim_hf_pr(nnnn,1) = exprnd(hist_Exp);
    end
    hfsim_data_coord_pr{ii} = [sim_x_coord sim_y_coord sim_hf_pr] ;
end

hfsim_lf_low_pr_coord = [];
hfsim_lf_med_pr_coord = [];
hfsim_lf_high_pr_coord = [];

for ii = 1:9
  low_ind = find(hfsim_data_coord_pr{ii}(:,3)<=low_lf_thresh_upper);
  med_ind = find(hfsim_data_coord_pr{ii}(:,3)>low_lf_thresh_upper & hfsim_data_coord_pr{ii}(:,3)<=med_lf_thresh_upper);
  high_ind = find(hfsim_data_coord_pr{ii}(:,3)>med_lf_thresh_upper);
  hfsim_lf_low_pr_coord{ii} = [hfsim_data_coord_pr{ii}(low_ind,1) hfsim_data_coord_pr{ii}(low_ind,2) hfsim_data_coord_pr{ii}(low_ind,3)];
  hfsim_lf_med_pr_coord{ii} = [hfsim_data_coord_pr{ii}(med_ind,1) hfsim_data_coord_pr{ii}(med_ind,2) hfsim_data_coord_pr{ii}(med_ind,3)];
  hfsim_lf_high_pr_coord{ii} = [hfsim_data_coord_pr{ii}(high_ind,1) hfsim_data_coord_pr{ii}(high_ind,2) hfsim_data_coord_pr{ii}(high_ind,3)];
end

%%Calculate G function for each group for each NMJ
hfSimGfunctionAll = []
for oooo = 1:9
    hfSimLowPoint = [hfsim_lf_low_pr_coord{oooo}(:,1) hfsim_lf_low_pr_coord{oooo}(:,2)];
    hfSimMedPoint = [hfsim_lf_med_pr_coord{oooo}(:,1) hfsim_lf_med_pr_coord{oooo}(:,2)];
    hfSimHighPoint = [hfsim_lf_high_pr_coord{oooo}(:,1) hfsim_lf_high_pr_coord{oooo}(:,2)]
    hfsimmaxD_L = computeMaxDistance(hfSimLowPoint);
    hfsimmaxD_M = computeMaxDistance(hfSimMedPoint);
    hfsimmaxD_H = computeMaxDistance(hfSimHighPoint);
    hfsimrRange_L = 0:1:hfsimmaxD_L;
    hfsimrRange_M = 0:1:hfsimmaxD_M;
    hfsimrRange_H = 0:1:hfsimmaxD_H;
    hfsimgFunct_L = computeGFunction(hfSimLowPoint, hfsimrRange_L);
    hfsimgFunct_M = computeGFunction(hfSimMedPoint, hfsimrRange_M);
    hfsimgFunct_H = computeGFunction(hfSimHighPoint, hfsimrRange_H);
    hfSimGfunctionAll{oooo}.Sim_GF_L = hfsimgFunct_L;
    hfSimGfunctionAll{oooo}.Sim_GF_M = hfsimgFunct_M;
    hfSimGfunctionAll{oooo}.Sim_GF_H = hfsimgFunct_H;
    hfSimGfunctionAll{oooo}.Sim_Range_L = hfsimgFunct_L;
    hfSimGfunctionAll{oooo}.Sim_Range_M = hfsimgFunct_M;
    hfSimGfunctionAll{oooo}.Sim_Range_H = hfsimgFunct_H;
end
%LowPrSim
hfSim_length_mat = double.empty;
for jjj = 1:9
hfsim_length_mat_temp = length(hfSimGfunctionAll{jjj}.Sim_GF_L);
hfSim_length_mat = [hfSim_length_mat hfsim_length_mat_temp];
end
hfsim_rRange_L_Min = min(hfSim_length_mat);

%%MedPrSim
hfSim_length_mat2 = double.empty;
for jjj = 1:9
hfsim_length_mat_temp2 = length(hfSimGfunctionAll{jjj}.Sim_GF_M);
hfSim_length_mat2 = [hfSim_length_mat2 hfsim_length_mat_temp2];
end
hfsim_rRange_M_Min = min(hfSim_length_mat2);

%HighPrSim
hfSim_length_mat3 = double.empty;
for jjj = 1:9
hfsim_length_mat_temp3 = length(hfSimGfunctionAll{jjj}.Sim_GF_H);
hfSim_length_mat3 = [hfSim_length_mat3 hfsim_length_mat_temp3];
end
hfsim_rRange_H_Min = min(hfSim_length_mat3);

hfSim_GfunctionLAll_Minimized = zeros(9,hfsim_rRange_L_Min);
hfSim_GfunctionMAll_Minimized = zeros(9,hfsim_rRange_M_Min);
hfSim_GfunctionHAll_Minimized = zeros(9,hfsim_rRange_H_Min);

for rrr = 1:9
hfSim_GfunctionLAll_Minimized(rrr,1:hfsim_rRange_L_Min) = hfSimGfunctionAll{rrr}.Sim_GF_L(1,1:hfsim_rRange_L_Min);
hfSim_GfunctionMAll_Minimized(rrr,1:hfsim_rRange_M_Min) = hfSimGfunctionAll{rrr}.Sim_GF_M(1,1:hfsim_rRange_M_Min);
hfSim_GfunctionHAll_Minimized(rrr,1:hfsim_rRange_H_Min) = hfSimGfunctionAll{rrr}.Sim_GF_H(1,1:hfsim_rRange_H_Min);

end

stdshade(hfSim_GfunctionLAll_Minimized,0.3);
hold on
stdshade(hfSim_GfunctionMAll_Minimized,0.3);
hold on
stdshade(hfSim_GfunctionHAll_Minimized,0.3); 

xlabel('Distance (r)');
ylabel('G(r)');
title('Simulated Nearest Neighbors G Function for Low Pr Sites');
grid on;

%%Plotall
%%%Plot both overlaid

figure(1)
stdshade(hfGfunctionLAll_Minimized,0.3,'m');
xlabel('Distance (pixels)');
ylabel('G(r)');
title('Nearest Neighbors G Function for Low Pr Sites');
grid off;
hold on
stdshade(hfSim_GfunctionLAll_Minimized,0.3,[.7 .7 .7]);
legend('Observed data Std Dev','Observed data mean','Simulated data Std Dev','Simulated data mean','Location', 'northwest' )
xlim([0 130])

figure(2)
stdshade(hfGfunctionMAll_Minimized,0.3,'m');
xlabel('Distance (pixels)');
ylabel('G(r)');
title('Nearest Neighbors G Function for Medium Pr Sites');
grid off;
hold on
stdshade(hfSim_GfunctionMAll_Minimized,0.3,[.7 .7 .7]);
legend('Observed data Std Dev','Observed data mean','Simulated data Std Dev','Simulated data mean','Location', 'northwest' )
xlim([0 130])

figure(3)
stdshade(hfGfunctionHAll_Minimized,0.3,'m');
xlabel('Distance (pixels)');
ylabel('G(r)');
title('Nearest Neighbors G Function for High Pr Sites');
grid off;
hold on
stdshade(hfSim_GfunctionHAll_Minimized,0.3,[.7 .7 .7]);
legend('Observed data Std Dev','Observed data mean','Simulated data Std Dev','Simulated data mean','Location', 'northwest' )
xlim([0 1200])

%% 3a Cluster data analysis (spatial pattern)

clear all
close all
clc

MC_interations=100;

data_all=cell(1);
S=struct;
data_all =All_Is_Data_No_Empty; 
%Load data
for i=1:9 %number of shtees
S(i).table=All_Is_Data_No_Empty{i};
end

th_dist=650;
for i=1:length(data_all) %for each NMJ
data=data_all{i};
    x=data(:,1);
    y=data(:,2);
    coord=[x y];
    plast_abs=data(:,4);
    plast_bin=data(:,5);
    plast02=data(:,6);
    plast5=data(:,7);

    plast=plast5;

Dist=pdist(coord);
Sq=squareform(Dist);

Z = linkage(coord,'ward');
T = cluster(Z,'cutoff',th_dist,'Criterion','distance'); %Index of cluster for each position
n_clusters=length(unique(T)); %Number of cluster

% Force a pattern (debugging, check the script works)
% T=sort(T);
% plast02=sort(plast02);
% plast5=sort(plast5);

% Get centroid of all clusters
centroid=[];
    for j=1:n_clusters
    idx=[];
    idx=find(T==j);

    cent_x=mean(x(idx)); cent_y=mean(y(idx));
    centroid(j,1)=cent_x; centroid(j,2)=cent_y;
    end

%Area and boundary of each cluster
k_c=[];area_c=[];x_events=[];y_events=[];
for j=1:n_clusters
    idx=[]; c=[];
    idx=find(T==j);
 
    x_tmp=x(idx,:);y_tmp=y(idx,:);
    x_events{j}=x_tmp; y_events{j}=y_tmp;
    
    if size(x_tmp,1)<=2%if there are only 2 events or less, use boundary instead of hull
    [k_c{j},area_c{j}] = boundary(x_tmp,y_tmp);
    else
    [k_c{j},area_c{j}] = convhull(x_tmp,y_tmp);  
    end
end

% Shuffled position
coord_rand=[];
rand_table=[];
for k=1:MC_interations %define number of iterations
r_idx=randperm(size(coord,1));
r_idx=r_idx';
%tmp=coord;
%coord_rand{k}=tmp([r_idx],:);
rand_table(:,k)=r_idx; %This is the index of the shuffled position, use it to extract Pr
end

% Extract values of each cluster, like mean Pr, Std Pr, Difference of 
% cluster Pr to global pr, difference between Pr 0.2 and 5hz

mean_Pr02=[];
std_Pr02=[];
mean_Pr5=[];
std_Pr5=[];
mean_Pr_abs=[];
std_Pr_abs=[];
fac_frac=[];

for n=1:n_clusters %for the number of cluster for each cell
Pr02_mean=plast02(find(T==n));
mean_Pr02(n)=mean(Pr02_mean);
std_Pr02(n)=std(Pr02_mean);

Pr5=plast5(find(T==n));
mean_Pr5(n)=mean(Pr5);
std_Pr5(n)=std(Pr5);

Pr_abs=plast_abs(find(T==n));
mean_Pr_abs(n)=mean(Pr_abs);
std_Pr_abs(n)=std(Pr_abs);
fac_frac(n)=length(find(Pr_abs>0))/length(Pr_abs);% facilitating fraction;
end


mean_Pr02_rand_all=[];
std_Pr02_rand_all=[];
mean_Pr5_rand_all=[];
std_Pr5_rand_all=[];
mean_Pr_abs_rand_all=[];
std_Pr_abs_rand_all=[];
fac_frac_rand_all=[];
%Extract values for suffled versions
for n=1:n_clusters 
    for k=1:size(rand_table,2)% number of MC iterations
    idx=rand_table(:,k);
        
    plast02_rand=plast02(idx);% this is where the index gets shuffled
    Pr02_mean_rand=plast02_rand(find(T==n));
    mean_Pr02_rand_all(n,k)=mean(Pr02_mean_rand);
    std_Pr02_rand_all(n,k)=std(Pr02_mean_rand);

    plast5_rand=plast5(idx);% this is where the index gets shuffled
    Pr5_rand=plast5_rand(find(T==n));
    mean_Pr5_rand_all(n,k)=mean(Pr5_rand);
    std_Pr5_rand_all(n,k)=std(Pr5_rand);

    plast_abs_rand=plast_abs(idx);% this is where the index gets shuffled
    Pr_abs_rand=plast_abs_rand(find(T==n));
    mean_Pr_abs_rand_all(n,k)=mean(Pr_abs_rand);
    std_Pr_abs_rand_all(n,k)=std(Pr_abs_rand);
    fac_frac_rand_all(n,k)=length(find(Pr_abs_rand>0))/length(Pr_abs_rand);% facilitating fraction;
    end
end
mean_Pr02_rand=[];
std_Pr02_rand=[];
mean_Pr5_rand=[];
std_Pr5_rand=[];
mean_Pr_abs_rand=[];
std_Pr_abs_rand=[];

mean_Pr02_rand=mean(mean_Pr02_rand_all,2);
std_Pr02_rand=mean(std_Pr02_rand_all,2);

mean_Pr5_rand=mean(mean_Pr5_rand_all,2);
std_Pr5_rand=mean(std_Pr5_rand_all,2);

mean_Pr_abs_rand=mean(mean_Pr_abs_rand_all,2);
std_Pr_abs_rand=mean(std_Pr_abs_rand_all,2);
fac_frac_rand=mean(fac_frac_rand_all,2);

%Pass values to struct
S(i).n_clusters=n_clusters;
S(i).i_clusters=T;
S(i).coord_rand=rand_table;

S(i).mean_Pr02=mean_Pr02;
S(i).std_Pr02=std_Pr02;
S(i).mean_Pr5=mean_Pr5;
S(i).std_Pr5=std_Pr5;
S(i).mean_Pr_abs=mean_Pr_abs;
S(i).std_Pr_abs=std_Pr_abs;
S(i).fac_frac=fac_frac;

S(i).mean_Pr02_rand_all=mean_Pr02_rand_all;
S(i).std_Pr02_rand_all=std_Pr02_rand_all;
S(i).mean_Pr5_rand_all=mean_Pr5_rand_all;
S(i).std_Pr5_rand_all=std_Pr5_rand_all;
S(i).mean_Pr_abs_rand_all=mean_Pr_abs_rand_all;
S(i).std_Pr_abs_rand_all=std_Pr_abs_rand_all;
S(i).fac_frac_rand_all=fac_frac_rand_all;

S(i).mean_Pr02_rand=mean_Pr02_rand;
S(i).std_Pr02_rand=std_Pr02_rand;
S(i).mean_Pr5_rand=mean_Pr5_rand;
S(i).std_Pr5_rand=std_Pr5_rand;
S(i).mean_Pr_abs_rand=mean_Pr_abs_rand;
S(i).std_Pr_abs_rand=std_Pr_abs_rand;
S(i).fac_frac_rand=fac_frac_rand;
end

%% 3b Plot data

close all
clearvars -except S
figure(1)
sgtitle('Is spatial plasticity analysis')

for nnnn = 1:9
    for c=nnnn%1:10 %choose cell to plot
%c=3; %cell

Pr02_mean=S(c).mean_Pr02;
Pr02_mean_rand=S(c).mean_Pr02_rand;

Pr5_mean=S(c).mean_Pr5;
Pr5_mean_rand=S(c).mean_Pr5_rand;

Pr_abs_mean=S(c).mean_Pr_abs;
Pr_abs_mean_rand=S(c).mean_Pr_abs_rand;

Pr02_std=S(c).std_Pr02;
Pr02_std_rand=S(c).std_Pr02_rand;

Pr5_std=S(c).std_Pr5;
Pr5_std_rand=S(c).std_Pr5_rand;

Pr_abs_std=S(c).std_Pr_abs;
Pr_abs_std_rand=S(c).std_Pr_abs_rand;

fac_frac=S(c).fac_frac;
fac_frac_rand=S(c).fac_frac_rand;

f1=figure(1);%------------------------------------------------------

subplot(3,3,1)
for i=1:size(Pr02_mean,2) %for each cluster
plot([1 2], [Pr02_mean(i) Pr02_mean_rand(i)],'b--o'); hold on
end
title('Mean Pr at 0.2 Hz')
xticks([1 2])
xticklabels({'Observed data','Simulated data'})

subplot(3,3,2)
for i=1:size(Pr5_mean,2) %for each cluster
plot([1 2], [Pr5_mean(i) Pr5_mean_rand(i)],'b--o'); hold on
end
title('Mean Pr at 5 Hz')
xticks([1 2])
xticklabels({'Observed data','Simulated data'})

subplot(3,3,3)
for i=1:size(Pr_abs_mean,2) %for each cluster
plot([1 2], [Pr_abs_mean(i) Pr_abs_mean_rand(i)],'b--o'); hold on
end
title('Mean Pr difference')
xticks([1 2])
xticklabels({'Observed data','Simulated data'})

subplot(3,3,4)
for i=1:size(Pr02_std,2) %for each cluster
plot([1 2], [Pr02_std(i) Pr02_std_rand(i)],'r--o'); hold on
end
title('Std Pr at 0.2 Hz')
xticks([1 2])
xticklabels({'Observed data','Simulated data'})

subplot(3,3,5)
for i=1:size(Pr5_std,2) %for each cluster
plot([1 2], [Pr5_std(i) Pr5_std_rand(i)],'r--o'); hold on

end
xticks([1 2])
xticklabels({'Observed data','Simulated data'})
title('Std Pr at 5 Hz')

subplot(3,3,6)
for i=1:size(Pr_abs_std,2) %for each cluster
plot([1 2], [Pr_abs_std(i) Pr_abs_std_rand(i)],'r--o'); hold on
end
title('Std Pr difference')
xticks([1 2])
xticklabels({'Observed data','Simulated data'})
subplot(3,3,7)
for i=1:size(fac_frac,2) %for each cluster
plot([1 2], [fac_frac(i) fac_frac_rand(i)],'g--o'); hold on
end
title('Fraction of Facilitating Synapses')
xticklabels({'Original data','MC version'})


figure(2)
for i=1:size(Pr02_mean,2) %for each cluster
plot([1 2], [Pr02_mean(i) Pr5_mean(i)],'b--o'); hold on
end
xticks([1 2])
xticklabels({'Observed data','Simulated data'})
    end
    hold on
end
title('Is spatial plasticity analysis')
figure(3)
scatter(coord(:,1),coord(:,2));

Pr02_mean_real_conc = [];
for kkkk = 1:9
ex = S(kkkk).mean_Pr02;
Pr02_mean_real_conc = [Pr02_mean_real_conc ex];
end
Pr02_mean_real_conc = Pr02_mean_real_conc.';

Pr5_mean_real_conc = [];
for kkkk = 1:9
ex = S(kkkk).mean_Pr5;
Pr5_mean_real_conc = [Pr5_mean_real_conc ex];
end
Pr5_mean_real_conc = Pr5_mean_real_conc.';

Pr_abs_mean_real_conc = [];
for kkkk = 1:9
ex = S(kkkk).mean_Pr_abs;
Pr_abs_mean_real_conc = [Pr_abs_mean_real_conc ex];
end
Pr_abs_mean_real_conc = Pr_abs_mean_real_conc.';

Pr02_std_mean_real_conc = [];
for kkkk = 1:9
ex = S(kkkk).std_Pr02;
Pr02_std_mean_real_conc = [Pr02_std_mean_real_conc ex];
end
Pr02_std_mean_real_conc = Pr02_std_mean_real_conc.';

Pr5_std_mean_real_conc = [];
for kkkk = 1:9
ex = S(kkkk).std_Pr5;
Pr5_std_mean_real_conc = [Pr5_std_mean_real_conc ex];
end
Pr5_std_mean_real_conc = Pr5_std_mean_real_conc.';

Pr_abs_std_real_conc = [];
for kkkk = 1:9
ex = S(kkkk).std_Pr_abs;
Pr_abs_std_real_conc = [Pr_abs_std_real_conc ex];
end
Pr_abs_std_real_conc = Pr_abs_std_real_conc.';

%Simulated data

Pr02_mean_rand_conc = [];
for kkkk = 1:9
ex = S(kkkk).mean_Pr02_rand;
Pr02_mean_rand_conc = [Pr02_mean_rand_conc ex.'];
end
Pr02_mean_rand_conc = Pr02_mean_rand_conc.';

Pr5_mean_rand_conc = [];
for kkkk = 1:9
ex = S(kkkk).mean_Pr5_rand;
Pr5_mean_rand_conc = [Pr5_mean_rand_conc ex.'];
end
Pr5_mean_rand_conc = Pr5_mean_rand_conc.';

Pr_abs_mean_rand_conc = [];
for kkkk = 1:9
ex = S(kkkk).mean_Pr_abs_rand;
Pr_abs_mean_rand_conc = [Pr_abs_mean_rand_conc ex.'];
end
Pr_abs_mean_rand_conc = Pr_abs_mean_rand_conc.';

Pr02_std_mean_rand_conc = [];
for kkkk = 1:9
ex = S(kkkk).std_Pr02_rand;
Pr02_std_mean_rand_conc = [Pr02_std_mean_rand_conc ex.'];
end
Pr02_std_mean_rand_conc = Pr02_std_mean_rand_conc.';

Pr5_std_mean_rand_conc = [];
for kkkk = 1:9
ex = S(kkkk).std_Pr5_rand;
Pr5_std_mean_rand_conc = [Pr5_std_mean_rand_conc ex.'];
end
Pr5_std_mean_rand_conc = Pr5_std_mean_rand_conc.';

Pr_abs_std_rand_conc = [];
for kkkk = 1:9
ex = S(kkkk).std_Pr_abs_rand;
Pr_abs_std_rand_conc = [Pr_abs_std_rand_conc ex.'];
end
Pr_abs_std_rand_conc = Pr_abs_std_rand_conc.';

for nnnn = 1:9
    for c=nnnn%1:10 %choose cell to plot
%c=3; %cell

Pr02_mean=S(c).mean_Pr02;
Pr02_mean_rand=S(c).mean_Pr02_rand;

Pr5_mean=S(c).mean_Pr5;
Pr5_mean_rand=S(c).mean_Pr5_rand;

Pr_abs_mean=S(c).mean_Pr_abs;
Pr_abs_mean_rand=S(c).mean_Pr_abs_rand;

Pr02_std=S(c).std_Pr02;
Pr02_std_rand=S(c).std_Pr02_rand;

Pr5_std=S(c).std_Pr5;
Pr5_std_rand=S(c).std_Pr5_rand;

Pr_abs_std=S(c).std_Pr_abs;
Pr_abs_std_rand=S(c).std_Pr_abs_rand;

fac_frac=S(c).fac_frac;
fac_frac_rand=S(c).fac_frac_rand;

for i=1:size(Pr_abs_mean,2) %for each cluster
plot([1 2], [Pr_abs_mean(i) Pr_abs_mean_rand(i)],'b--o'); hold on
end

Pr02_mean_real_conc = [Pr02_mean_real_conc;
Pr5_mean_real_conc = [];
Pr_abs_mean_real_conc = [];
Pr02_std_mean_real_conc = [];
Pr5_std_mean_real_conc = [];
Pr02_abs_std_real_conc = [];

Pr02_mean_rand_conc = [];
Pr5_mean_rand_conc = [];
Pr_abs_mean_rand_conc = [];
Pr02_std_mean_rand_conc = [];
Pr5_std_mean_rand_conc = [];
Pr02_abs_std_rand_conc = [];
    end
end

writetable(struct2table(S), 'Is_Plasticity_Spatial.csv')
