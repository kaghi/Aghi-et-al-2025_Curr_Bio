allpr = table2array(All5HzTrainsZach);
allpr = [total_lfpr allpr];
sortedpr = sortrows(allpr,1);
labelcolumn = ones(length(sortedpr),6);
for ii = 1:length(labelcolumn)
    labelcolumn(ii,2) = 2;
    labelcolumn(ii,3) = 3;
    labelcolumn(ii,4) = 4;
    labelcolumn(ii,5) = 5;
    labelcolumn(ii,6) = 6;
end
sortedlfpr =  sortedpr(:,1);
sortedhfpr1 = sortedpr(:,2);
sortedhfpr2 = sortedpr(:,3);
sortedhfpr3 = sortedpr(:,4);
sortedhfpr4 = sortedpr(:,5);
sortedhfpr5 = sortedpr(:,6);

figure
hold on
scatter(labelcolumn(:,1),sortedlfpr,'filled')
scatter(labelcolumn(:,2),sortedhfpr1,'filled');
scatter(labelcolumn(:,3),sortedhfpr2,'filled');
scatter(labelcolumn(:,4),sortedhfpr3,'filled')
scatter(labelcolumn(:,5),sortedhfpr4,'filled');
scatter(labelcolumn(:,6),sortedhfpr5,'filled');
x1 = labelcolumn(:,1);
x2 = labelcolumn(:,2);
x3 = labelcolumn(:,3);
x4 = labelcolumn(:,4);
x5 = labelcolumn(:,5);
x6 = labelcolumn(:,6);
y1 = sortedlfpr;
y2 = sortedhfpr1;
y3 = sortedhfpr2;
y4 = sortedhfpr3;
y5 = sortedhfpr4;
y6 = sortedhfpr5;
plot([x1(:)';x2(:)';x3(:)';x4(:)';x5(:)';x6(:)'], [y1(:)';y2(:)';y3(:)';y4(:)';y5(:)';y6(:)'], 'k-')

%%%Divide based on every x values
%%Want every 218 values
indie = [];
for ii = 1:6
    temp = 218;
    val = temp*ii;
    indie(ii) = val;
end
indie = [1 indie];
colormat = ["#1F2041", "#4B3F72","#FFC857","#119DA4","#19647E","#6D72C3"];

for mmm = 1:6
ind = indie(mmm);
ind2 = indie(mmm+1);
if ind2<length(labelcolumn);
figure(mmm)
hold on
scatter(labelcolumn(ind:ind2,1),sortedlfpr(ind:ind2,1),2,'filled',"MarkerEdgeColor",colormat(mmm),"MarkerFaceColor",colormat(mmm))
scatter(labelcolumn(ind:ind2,2),sortedhfpr1(ind:ind2,1),2,'filled',"MarkerEdgeColor",colormat(mmm),"MarkerFaceColor",colormat(mmm));
scatter(labelcolumn(ind:ind2,3),sortedhfpr2(ind:ind2,1),2,'filled',"MarkerEdgeColor",colormat(mmm),"MarkerFaceColor",colormat(mmm));
scatter(labelcolumn(ind:ind2,4),sortedhfpr3(ind:ind2,1),2,'filled',"MarkerEdgeColor",colormat(mmm),"MarkerFaceColor",colormat(mmm))
scatter(labelcolumn(ind:ind2,5),sortedhfpr4(ind:ind2,1),2,'filled',"MarkerEdgeColor",colormat(mmm),"MarkerFaceColor",colormat(mmm));
scatter(labelcolumn(ind:ind2,6),sortedhfpr5(ind:ind2,1),2,'filled',"MarkerEdgeColor",colormat(mmm),"MarkerFaceColor",colormat(mmm));
set(gca,'xticklabel',[])
set(gcf,'paperunits','inches')
set(gcf,'position',[0 0 200 400])
ax = gca;
ax.FontSize = 15; 
ylim([0 0.6])
x1 = labelcolumn(ind:ind2,1);
x2 = labelcolumn(ind:ind2,2);
x3 = labelcolumn(ind:ind2,3);
x4 = labelcolumn(ind:ind2,4);
x5 = labelcolumn(ind:ind2,5);
x6 = labelcolumn(ind:ind2,6);
y1 = sortedlfpr(ind:ind2,1);
y2 = sortedhfpr1(ind:ind2,1);
y3 = sortedhfpr2(ind:ind2,1);
y4 = sortedhfpr3(ind:ind2,1);
y5 = sortedhfpr4(ind:ind2,1);
y6 = sortedhfpr5(ind:ind2,1);
plot([x1(:)';x2(:)';x3(:)';x4(:)';x5(:)';x6(:)'], [y1(:)';y2(:)';y3(:)';y4(:)';y5(:)';y6(:)'], 'k-','Color',colormat(mmm))
else
ind2 = length(labelcolumn);
figure(mmm)
hold on
scatter(labelcolumn(ind:ind2,1),sortedlfpr(ind:ind2,1),2,'filled',"MarkerEdgeColor",colormat(mmm),"MarkerFaceColor",colormat(mmm));
scatter(labelcolumn(ind:ind2,2),sortedhfpr1(ind:ind2,1),2,'filled',"MarkerEdgeColor",colormat(mmm),"MarkerFaceColor",colormat(mmm));
scatter(labelcolumn(ind:ind2,3),sortedhfpr2(ind:ind2,1),2,'filled',"MarkerEdgeColor",colormat(mmm),"MarkerFaceColor",colormat(mmm));
scatter(labelcolumn(ind:ind2,4),sortedhfpr3(ind:ind2,1),2,'filled',"MarkerEdgeColor",colormat(mmm),"MarkerFaceColor",colormat(mmm));
scatter(labelcolumn(ind:ind2,5),sortedhfpr4(ind:ind2,1),2,'filled',"MarkerEdgeColor",colormat(mmm),"MarkerFaceColor",colormat(mmm));
scatter(labelcolumn(ind:ind2,6),sortedhfpr5(ind:ind2,1),2,'filled',"MarkerEdgeColor",colormat(mmm),"MarkerFaceColor",colormat(mmm));
set(gca,'xticklabel',[])
set(gcf,'paperunits','inches')
set(gcf,'position',[0 0 200 400])
ax = gca;
ax.FontSize = 15; 
ylim([0 0.6])
x1 = labelcolumn(ind:ind2,1);
x2 = labelcolumn(ind:ind2,2);
x3 = labelcolumn(ind:ind2,3);
x4 = labelcolumn(ind:ind2,4);
x5 = labelcolumn(ind:ind2,5);
x6 = labelcolumn(ind:ind2,6);
y1 = sortedlfpr(ind:ind2,1);
y2 = sortedhfpr1(ind:ind2,1);
y3 = sortedhfpr2(ind:ind2,1);
y4 = sortedhfpr3(ind:ind2,1);
y5 = sortedhfpr4(ind:ind2,1);
y6 = sortedhfpr5(ind:ind2,1);
plot([x1(:)';x2(:)';x3(:)';x4(:)';x5(:)';x6(:)'], [y1(:)';y2(:)';y3(:)';y4(:)';y5(:)';y6(:)'], 'k-','Color',colormat(mmm))
end
end

%%%%Look at average behavior per group
delta_pr = sortedhfpr1 - sortedlfpr;
delta_pr_groups = cell(1,6);
for nn = 1:6
    ind = indie(nn);
    ind2 = indie(nn+1);
    delta_pr_groups{nn} = delta_pr(ind:ind2,1);
end

range_vals = sortedlfpr(indie,1);
