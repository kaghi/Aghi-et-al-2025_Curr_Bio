%Remove zero rows 
Non_Zero_Ind = find(Pre_D_Ib_Rank(:,1)>0);
Pre_D_Ranked_Nonzero(:,1) =Pre_D_Ib_Rank(Non_Zero_Ind,1);
Pre_D_Ranked_Nonzero(:,2) =Pre_D_Ib_Rank(Non_Zero_Ind,2);
raw_change = (Pre_D_Ranked_Nonzero(:,2)-Pre_D_Ranked_Nonzero(:,1));

for nn = 1:length(raw_change);
per_change(nn,1) = raw_change(nn)/(Pre_D_Ranked_Nonzero(nn,1));
end

per_change = per_change*100;
pre_pr_and_pr_change = [Pre_D_Ranked_Nonzero(:,1) per_change];


%Attp2
%Remove zero rows 

mattat = Attp2_Ib;
Ranked_Attp2 = sortrows(mattat,1);
Non_Zero_Ind_A = find(Ranked_Attp2(:,1)>0);
Pre_A_Ranked_Nonzero(:,1) =Ranked_Attp2(Non_Zero_Ind_A,1);
Pre_A_Ranked_Nonzero(:,2) =Ranked_Attp2(Non_Zero_Ind_A,2);
raw_change_A = (Pre_A_Ranked_Nonzero(:,2)-Pre_A_Ranked_Nonzero(:,1));

for nn = 1:length(raw_change_A);
per_changeA(nn,1) = raw_change_A(nn)/(Pre_A_Ranked_Nonzero(nn,1));
end

per_changeA = per_changeA*100;
pre_pr_and_pr_changeA = [Pre_A_Ranked_Nonzero(:,1) per_changeA];


scatter(pre_pr_and_pr_changeA(:,1),pre_pr_and_pr_changeA(:,2),5,'Filled')
hold on 
scatter(pre_pr_and_pr_change(:,1),pre_pr_and_pr_change(:,2),5,'Filled')

pre_A_x = pre_pr_and_pr_changeA(:,1);
pre_A_y = pre_pr_and_pr_changeA(:,2);

X = pre_pr_and_pr_changeA;
A = X(:,1);
MC = arrayfun(@(i) mean(X(A == i, 2:end)), unique(A), 'UniformOutput', false);
M = reshape(cell2mat(MC), size(X,2)-1, [])'
C = unique(A)

[ids, ~, rows] = unique(X(:, 1));
stddevs1 = accumarray(rows, X(:, 2), [], @std);

for nn = 1:length(M);
envelopesstdv1(nn,1) = [M(nn,1)-stddevs1(nn,1)]
envelopesstdv1(nn,2) = [M(nn,1)+stddevs1(nn,1)]
end

x3 = [C;flip(C)];
inBetween2 = [envelopesstdv1(:,1); flip(envelopesstdv1(:,2))];



X2 = pre_pr_and_pr_change;
A2 = X2(:,1);
MC2 = arrayfun(@(i) mean(X2(A2 == i, 2:end)), unique(A2), 'UniformOutput', false);
M2 = reshape(cell2mat(MC2), size(X2,2)-1, [])'
C2 = unique(A2)
plot(C2,M2)
hold on 
plot(envelopestdv)
envelopestdv


[ids, ~, rows] = unique(X2(:, 1));
stddevs = accumarray(rows, X2(:, 2), [], @std);

for nn = 1:length(M2);
envelopestdv(nn,1) = [M2(nn,1)-stddevs(nn,1)]
envelopestdv(nn,2) = [M2(nn,1)+stddevs(nn,1)]
end

x2 = [C2;flip(C2)];
inBetween = [envelopestdv(:,1); flip(envelopestdv(:,2))];

%plot it all
fill(x3, inBetween2, [0.4660 0.6740 0.1880],'LineStyle','none','FaceAlpha',0.2);
hold on;
plot(C,M,'color',[0.4660 0.6740 0.1880], 'LineWidth', 2);

hold on 

fill(x2, inBetween, [0.4940 0.1840 0.5560],'LineStyle','none','FaceAlpha',0.2);
hold on;
plot(C2,M2,'color',[0.4940 0.1840 0.5560], 'LineWidth', 2);

legend('attP2 stdev','attP2','DmGluRA-RNAi stdev','DmGluRA-RNAi')
ylabel('Percent change')
xlabel('Pre-LY354740 Pr')


%Running average 

M4 = movmean(M,4)
plot(C,M4);
M5 = movmean(M2,4)
plot(C2,M5);
hold on 

env1_sm = movmean(envelopestdv(:,1),3);
env2_sm = movmean(envelopestdv(:,2),3);

InBetween_smoothed_1 = [env1_sm; flip(env2_sm)];

fill(x2, InBetween_smoothed_1, [0.4660 0.6740 0.1880],'LineStyle','none','FaceAlpha',0.2);
hold on;
plot(C2,M4,'color',[0.4660 0.6740 0.1880], 'LineWidth', 2);

