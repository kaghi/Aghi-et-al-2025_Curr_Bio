DmGluRA_Ib = table2array(MergedDmGluRARNAiIb);
DmGluRA_Is = table2array(MergedDmGluRARNAiIs);
Attp2_Ib = table2array(MergedAttp2Ib);
Attp2_Is = table2array(MergedAttp2Is);

%%DmGluRA first
%%First divide into quartiles based on pre agonist pr 
DB_Q1 = 0;
DB_Q2 = 0.02;
DB_Q3 = 0.05;
DS_Q1 = 0.007500;
DS_Q2 = 0.03000;
DS_Q3 = 0.09250;

DB_Q1_ind = find(DmGluRA_Ib(:,1)==DB_Q1);
DB_Q2_ind = find(DmGluRA_Ib(:,1)>DB_Q1 & DmGluRA_Ib(:,1)<=DB_Q2);
DB_Q3_ind =find(DmGluRA_Ib(:,1)>DB_Q2 & DmGluRA_Ib(:,1)<=DB_Q3);
DB_Q4_ind =find(DmGluRA_Ib(:,1)>DB_Q3);

DS_Q1_ind = find(DmGluRA_Is(:,1)<=DS_Q1);
DS_Q2_ind = find(DmGluRA_Is(:,1)>DS_Q1 & DmGluRA_Is(:,1)<=DS_Q2);
DS_Q3_ind =find(DmGluRA_Is(:,1)>DS_Q2 & DmGluRA_Is(:,1)<=DS_Q3);
DS_Q4_ind =find(DmGluRA_Is(:,1)>DS_Q3);

Pr_DB_Q1 = [DmGluRA_Ib(DB_Q1_ind,1) DmGluRA_Ib(DB_Q1_ind,2)] ;
Pr_DB_Q2 = [DmGluRA_Ib(DB_Q2_ind,1) DmGluRA_Ib(DB_Q2_ind,2)] ;
Pr_DB_Q3 = [DmGluRA_Ib(DB_Q3_ind,1) DmGluRA_Ib(DB_Q3_ind,2)] ;
Pr_DB_Q4 = [DmGluRA_Ib(DB_Q4_ind,1) DmGluRA_Ib(DB_Q4_ind,2)] ;

Pr_DS_Q1 = [DmGluRA_Is(DS_Q1_ind,1) DmGluRA_Is(DS_Q1_ind,2)];
Pr_DS_Q2 = [DmGluRA_Is(DS_Q2_ind,1) DmGluRA_Is(DS_Q2_ind,2)] ;
Pr_DS_Q3 = [DmGluRA_Is(DS_Q3_ind,1) DmGluRA_Is(DS_Q3_ind,2)] ;;
Pr_DS_Q4 = [DmGluRA_Is(DS_Q4_ind,1) DmGluRA_Is(DS_Q4_ind,2)] ;

Diff_DB_Q1 = Pr_DB_Q1(:,2) - Pr_DB_Q1(:,1);
Diff_DB_Q2 = Pr_DB_Q2(:,2) - Pr_DB_Q2(:,1);
Diff_DB_Q3 = Pr_DB_Q3(:,2) - Pr_DB_Q3(:,1);
Diff_DB_Q4 = Pr_DB_Q4(:,2) - Pr_DB_Q4(:,1);

Diff_DS_Q1 = Pr_DS_Q1(:,2) - Pr_DS_Q1(:,1);
Diff_DS_Q2 = Pr_DS_Q2(:,2) - Pr_DS_Q2(:,1);
Diff_DS_Q3 = Pr_DS_Q3(:,2) - Pr_DS_Q3(:,1);
Diff_DS_Q4 = Pr_DS_Q4(:,2) - Pr_DS_Q4(:,1);

%%Attp2 next
%%First divide into quartiles based on pre agonist pr 
AB_Q1 = 0.01;
AB_Q2 = 0.03;
AB_Q3 = 0.06;
AS_Q1 = 0.01;
AS_Q2 = 0.05000;
AS_Q3 = 0.1175;

AB_Q1_ind = find(Attp2_Ib(:,1)<=AB_Q1);
AB_Q2_ind = find(Attp2_Ib(:,1)>AB_Q1 & Attp2_Ib(:,1)<=AB_Q2);
AB_Q3_ind =find(Attp2_Ib(:,1)>AB_Q2 & Attp2_Ib(:,1)<=AB_Q3);
AB_Q4_ind =find(Attp2_Ib(:,1)>AB_Q3);

AS_Q1_ind = find(Attp2_Is(:,1)<=AS_Q1);
AS_Q2_ind = find(Attp2_Is(:,1)>AS_Q1 & Attp2_Is(:,1)<=AS_Q2);
AS_Q3_ind =find(Attp2_Is(:,1)>AS_Q2 & Attp2_Is(:,1)<=AS_Q3);
AS_Q4_ind =find(Attp2_Is(:,1)>AS_Q3);

Pr_AB_Q1 = [Attp2_Ib(AB_Q1_ind,1) Attp2_Ib(AB_Q1_ind,2)] ;
Pr_AB_Q2 = [Attp2_Ib(AB_Q2_ind,1) Attp2_Ib(AB_Q2_ind,2)] ;
Pr_AB_Q3 = [Attp2_Ib(AB_Q3_ind,1) Attp2_Ib(AB_Q3_ind,2)] ;
Pr_AB_Q4 = [Attp2_Ib(AB_Q4_ind,1) Attp2_Ib(AB_Q4_ind,2)] ;

Pr_AS_Q1 = [Attp2_Is(AS_Q1_ind,1) Attp2_Is(AS_Q1_ind,2)];
Pr_AS_Q2 = [Attp2_Is(AS_Q2_ind,1) Attp2_Is(AS_Q2_ind,2)] ;
Pr_AS_Q3 = [Attp2_Is(AS_Q3_ind,1) Attp2_Is(AS_Q3_ind,2)] ;;
Pr_AS_Q4 = [Attp2_Is(AS_Q4_ind,1) Attp2_Is(AS_Q4_ind,2)] ;

Diff_AB_Q1 = Pr_AB_Q1(:,2) - Pr_AB_Q1(:,1);
Diff_AB_Q2 = Pr_AB_Q2(:,2) - Pr_AB_Q2(:,1);
Diff_AB_Q3 = Pr_AB_Q3(:,2) - Pr_AB_Q3(:,1);
Diff_AB_Q4 = Pr_AB_Q4(:,2) - Pr_AB_Q4(:,1);

Diff_AS_Q1 = Pr_AS_Q1(:,2) - Pr_AS_Q1(:,1);
Diff_AS_Q2 = Pr_AS_Q2(:,2) - Pr_AS_Q2(:,1);
Diff_AS_Q3 = Pr_AS_Q3(:,2) - Pr_AS_Q3(:,1);
Diff_AS_Q4 = Pr_AS_Q4(:,2) - Pr_AS_Q4(:,1);

writematrix(Diff_DB_Q1,'DmGluRA_Quartiles_Differences_Q1_Ib.csv');
writematrix(Diff_DB_Q2,'DmGluRA_Quartiles_Differences.csv','WriteMode','append');




% Pre_D_Ib_Rank = sortrows(DmGluRA_Ib,1);
% counter = 1;
% for ii = 1:length(DmGluRA_Ib)-1;
%     if Pre_D_Ib_Rank(ii,1) == Pre_D_Ib_Rank(ii+1,1);
%         Pre_D_Ib_Rank(ii,3) = counter;
%     else
%         Pre_D_Ib_Rank(ii,3) = counter;
%         counter = counter + 1;
%     end
% end
% Pre_D_Ib_Rank(length(DmGluRA_Ib),3) = counter;
% 
% Post_D_Ib_Rank = sortrows(DmGluRA_Ib,2);
% counter = 1;
% for ii = 1:length(DmGluRA_Ib)-1;
%     if Post_D_Ib_Rank(ii,2) == Post_D_Ib_Rank(ii+1,2);
%         Post_D_Ib_Rank(ii,3) = counter;
%     else
%         Post_D_Ib_Rank(ii,3) = counter;
%         counter = counter + 1;
%     end
% end
% Post_D_Ib_Rank(length(DmGluRA_Ib),3) = counter;
% 
% RankDiff_D_Ib = [];
% for ii = 1:length(DmGluRA_Ib)
%     RankDiff_D_Ib(ii,1) = Post_D_Ib_Rank(ii,3) - Pre_D_Ib_Rank(ii,3); 
% end
% avg_RankDiff_D_Ib = mean(RankDiff_D_Ib);
% std_RankDiff_D_Ib = std(RankDiff_D_Ib);
% Pr_change_D_Ib = DmGluRA_Ib(:,2) - DmGluRA_Ib(:,1); 
% avg_Pr_change_D_Ib = mean(Pr_change_D_Ib);
% std_Pr_change_D_Ib = std(Pr_change_D_Ib); 
%  
% Pr_change_D_Is = DmGluRA_Is(:,2) - DmGluRA_Is(:,1); 
% avg_Pr_change_D_Is = mean(Pr_change_D_Is);
% std_Pr_change_D_Is = std(Pr_change_D_Is); 
% 
% Pr_change_A_Ib = Attp2_Ib(:,2) - Attp2_Ib(:,1); 
% avg_Pr_change_A_Ib = mean(Pr_change_A_Ib);
% std_Pr_change_A_Ib = std(Pr_change_A_Ib); 
% 
% Pr_change_A_Is = Attp2_Is(:,2) - Attp2_Is(:,1); 
% avg_Pr_change_A_Is = mean(Pr_change_A_Is);
% std_Pr_change_A_Is = std(Pr_change_A_Is);

ecdf(DmGluRA_Ib(:,1))

%% %bin within quintile and compute rank order 
