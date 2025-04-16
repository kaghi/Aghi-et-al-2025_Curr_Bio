lfpr_cell = [];
hfpr_cell = []
lengths_lfprcell = []
for ll = 1:10
    lfpr_cell{ll} = data_all{ll}(:,6);
    hfpr_cell{ll} = data_all{ll}(:,7)
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
for ii = 1:10
    data_coord_pr{ii} = [data_all{ii}(:,1) data_all{ii}(:,2) data_all{ii}(:,6) data_all{ii}(:,7)] ;
end
hf_low_pr_coord = [];
hf_med_pr_coord = [];
hf_high_pr_coord = [];

for ii = 1:10
  low_ind = find(data_all{ii}(:,7)<=low_lf_thresh_upper);
  med_ind = find(data_all{ii}(:,7)>low_lf_thresh_upper & data_all{ii}(:,6)<=med_lf_thresh_upper);
  high_ind = find(data_all{ii}(:,7)>med_lf_thresh_upper);
  hf_low_pr_coord{ii} = [data_all{ii}(low_ind,1) data_all{ii}(low_ind,2) data_all{ii}(low_ind,6) data_all{ii}(low_ind,7)];
  hf_med_pr_coord{ii} = [data_all{ii}(med_ind,1) data_all{ii}(med_ind,2) data_all{ii}(med_ind,6) data_all{ii}(med_ind,7)];
  hf_high_pr_coord{ii} = [data_all{ii}(high_ind,1) data_all{ii}(high_ind,2) data_all{ii}(high_ind,6) data_all{ii}(high_ind,7)];
end
%%Calculate G function for each group for each NMJ
hGfunctionaAll = []
for oooo = 1:10
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
for jjj = 1:10
length_mat_temp1 = length(hGfunctionaAll{jjj}.GF_L);
hflength_mat = [hflength_mat length_mat_temp1];
end
for jjj = 1:10
length_mat_temp2 = length(hGfunctionaAll{jjj}.GF_M);
hflength_matM = [hflength_matM length_mat_temp2];
end
for jjj = 1:10
length_mat_temp3 = length(hGfunctionaAll{jjj}.GF_H);
hflength_matH = [hflength_matH length_mat_temp3];
end
hfrRange_L_Min = min(hflength_mat);
hfrRange_M_Min = min(hflength_matM);
hfrRange_H_Min = min(hflength_matH);

hfGfunctionLAll_Minimized = zeros(10,hfrRange_L_Min);
hfGfunctionMAll_Minimized = zeros(10,hfrRange_M_Min);
hfGfunctionHAll_Minimized = zeros(10,hfrRange_H_Min);


for rrr = 1:10
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
for ii = 1:10
    sim_x_coord = zeros(length(data_all{ii}(:,1)),1);
    sim_y_coord = zeros(length(data_all{ii}(:,1)),1);
    sim_hf_pr = zeros(length(data_all{ii}(:,1)),1);
    for nnnn = 1:length(data_all{ii}(:,1))
        sim_x_coord(nnnn,1) = data_all{ii}(nnnn,1);
        sim_y_coord(nnnn,1) = data_all{ii}(nnnn,2);
        sim_hf_pr(nnnn,1) = exprnd(hist_Exp);
    end
    hfsim_data_coord_pr{ii} = [sim_x_coord sim_y_coord sim_hf_pr] ;
end

hfsim_lf_low_pr_coord = [];
hfsim_lf_med_pr_coord = [];
hfsim_lf_high_pr_coord = [];

for ii = 1:10
  low_ind = find(hfsim_data_coord_pr{ii}(:,3)<=low_lf_thresh_upper);
  med_ind = find(hfsim_data_coord_pr{ii}(:,3)>low_lf_thresh_upper & hfsim_data_coord_pr{ii}(:,3)<=med_lf_thresh_upper);
  high_ind = find(hfsim_data_coord_pr{ii}(:,3)>med_lf_thresh_upper);
  hfsim_lf_low_pr_coord{ii} = [hfsim_data_coord_pr{ii}(low_ind,1) hfsim_data_coord_pr{ii}(low_ind,2) hfsim_data_coord_pr{ii}(low_ind,3)];
  hfsim_lf_med_pr_coord{ii} = [hfsim_data_coord_pr{ii}(med_ind,1) hfsim_data_coord_pr{ii}(med_ind,2) hfsim_data_coord_pr{ii}(med_ind,3)];
  hfsim_lf_high_pr_coord{ii} = [hfsim_data_coord_pr{ii}(high_ind,1) hfsim_data_coord_pr{ii}(high_ind,2) hfsim_data_coord_pr{ii}(high_ind,3)];
end

%%Calculate G function for each group for each NMJ
hfSimGfunctionAll = []
for oooo = 1:10
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

hfSim_length_mat = double.empty;
for jjj = 1:10
hfsim_length_mat_temp = length(hfSimGfunctionAll{jjj}.Sim_GF_L);
hfSim_length_mat = [hfSim_length_mat hfsim_length_mat_temp];
end
hfsim_rRange_L_Min = min(hfSim_length_mat);

%%MedPrSim
hfSim_length_mat2 = double.empty;
for jjj = 1:10
hfsim_length_mat_temp2 = length(hfSimGfunctionAll{jjj}.Sim_GF_M);
hfSim_length_mat2 = [hfSim_length_mat2 hfsim_length_mat_temp2];
end
hfsim_rRange_M_Min = min(hfSim_length_mat2);

%HighPrSim
hfSim_length_mat3 = double.empty;
for jjj = 1:10
hfsim_length_mat_temp3 = length(hfSimGfunctionAll{jjj}.Sim_GF_H);
hfSim_length_mat3 = [hfSim_length_mat3 hfsim_length_mat_temp3];
end
hfsim_rRange_H_Min = min(hfSim_length_mat3);

hfSim_GfunctionLAll_Minimized = zeros(10,hfsim_rRange_L_Min);
hfSim_GfunctionMAll_Minimized = zeros(10,hfsim_rRange_M_Min);
hfSim_GfunctionHAll_Minimized = zeros(10,hfsim_rRange_H_Min);

for rrr = 1:10
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
xlim([0 1200])
figure(2)
stdshade(hfGfunctionMAll_Minimized,0.3,'m');
xlabel('Distance (pixels)');
ylabel('G(r)');
title('Nearest Neighbors G Function for Medium Pr Sites');
grid off;
hold on
stdshade(hfSim_GfunctionMAll_Minimized,0.3,[.7 .7 .7]);
legend('Observed data Std Dev','Observed data mean','Simulated data Std Dev','Simulated data mean','Location', 'northwest' )
xlim([0 1200])

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

