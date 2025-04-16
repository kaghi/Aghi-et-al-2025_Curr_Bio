
Attp2_Ib_Sorted = sortrows(Attp2_Ib,1);
Attp2_Is_Sorted = sortrows(Attp2_Is,1);
DmGluRA_Ib_Sorted = sortrows(DmGluRA_Ib,1);
DmGluRA_Is_Sorted = sortrows(DmGluRA_Is,1);

%Attp2 Ib
ABSortInd = 1:(length(Attp2_Ib_Sorted)/5):(length(Attp2_Ib_Sorted)+1);
ABSortFifths =[];
for ii = 1:5
    start_ind = ABSortInd(ii);
    end_ind = ABSortInd(ii+1)-1;
    ABSortFifths{ii} = Attp2_Ib_Sorted(start_ind:end_ind,:)
end
labelcolumn1 = ones(length(ABSortFifths{1}),2);
labelcolumn1(:,2) = labelcolumn1(:,1)+labelcolumn1(:,2);

absize = size(ABSortFifths);
colormat = ["#1F2041", "#4B3F72","#FFC857","#119DA4","#19647E"]
tiledlayout(1,5);
for iiii = 1:5
nexttile
hold on
x1 = labelcolumn1(:,1);
y1 = ABSortFifths{iiii}(:,1);
x2 = labelcolumn1(:,2);
y2 = ABSortFifths{iiii}(:,2);
scatter(labelcolumn1(:,1),ABSortFifths{iiii}(:,1),2,'filled',"MarkerEdgeColor",colormat(iiii),"MarkerFaceColor",colormat(iiii));
scatter(labelcolumn1(:,2),ABSortFifths{iiii}(:,2),2,'filled',"MarkerEdgeColor",colormat(iiii),"MarkerFaceColor",colormat(iiii));
plot([x1(:)';x2(:)'], [y1(:)';y2(:)'], 'k-','Color',colormat(iiii),'LineWidth',0.1 )
xticks([1 2])
xticklabels({'Pre-LY354740','Post-LY354740'})
ylim([0 0.6])
end

%Attp2 Is
AsSortInd = 1:(length(Attp2_Is_Sorted)/5):(length(Attp2_Is_Sorted)+1);
AsSortInd = floor(AsSortInd);
AsSortFifths =[];
for ii = 1:5
    start_ind = AsSortInd(ii);
    end_ind = AsSortInd(ii+1)-1;
    AsSortFifths{ii} = Attp2_Is_Sorted(start_ind:end_ind,:)
end
maxlengthAS = []
for nnnn = 1:5
maxlengthAS(nnnn) = length(AsSortFifths{nnnn});
end
maxmaxind = max(maxlengthAS);
labelcolumn1 = ones(maxmaxind,2);
labelcolumn1(:,2) = labelcolumn1(:,1)+labelcolumn1(:,2);

absize = size(AsSortFifths);
colormat = ["#1F2041", "#4B3F72","#FFC857","#119DA4","#19647E"]
tiledlayout(1,5);
for iiii = 1:5
nexttile
hold on
xind = length(AsSortFifths{iiii}(:,1));
x1 = labelcolumn1(1:xind,1);
y1 = AsSortFifths{iiii}(:,1);
x2 = labelcolumn1(1:xind,2);
y2 = AsSortFifths{iiii}(:,2);
scatter(labelcolumn1(1:xind,1),AsSortFifths{iiii}(:,1),2,'filled',"MarkerEdgeColor",colormat(iiii),"MarkerFaceColor",colormat(iiii));
scatter(labelcolumn1(1:xind,2),AsSortFifths{iiii}(:,2),2,'filled',"MarkerEdgeColor",colormat(iiii),"MarkerFaceColor",colormat(iiii));
plot([x1(:)';x2(:)'], [y1(:)';y2(:)'], 'k-','Color',colormat(iiii),'LineWidth',0.25 )
xticks([1 2])
xticklabels({'Pre-LY354740','Post-LY354740'})
ylim([0 0.6])
end

%DmGluRA-RNAi Ib
DBSortInd = 1:(length(DmGluRA_Ib_Sorted)/5):(length(DmGluRA_Ib_Sorted)+1);
DBSortInd = floor(DBSortInd);
DBSortFifths =[];
for ii = 1:5
    start_ind = DBSortInd(ii);
    end_ind = DBSortInd(ii+1)-1;
    DBSortFifths{ii} = DmGluRA_Ib_Sorted(start_ind:end_ind,:)
end
maxlengthAB = []
for nnnn = 1:5
maxlengthAB(nnnn) = length(DBSortFifths{nnnn});
end
maxmaxind = max(maxlengthAB);
labelcolumn1 = ones(maxmaxind,2);
labelcolumn1(:,2) = labelcolumn1(:,1)+labelcolumn1(:,2);

dbsize = size(DBSortFifths);
colormat = ["#1F2041", "#4B3F72","#FFC857","#119DA4","#19647E"]
tiledlayout(1,5);
for iiii = 1:5
nexttile
hold on
xind = length(DBSortFifths{iiii}(:,1));
x1 = labelcolumn1(1:xind,1);
y1 = DBSortFifths{iiii}(:,1);
x2 = labelcolumn1(1:xind,2);
y2 = DBSortFifths{iiii}(:,2);
scatter(labelcolumn1(1:xind,1),DBSortFifths{iiii}(:,1),2,'filled',"MarkerEdgeColor",colormat(iiii),"MarkerFaceColor",colormat(iiii));
scatter(labelcolumn1(1:xind,2),DBSortFifths{iiii}(:,2),2,'filled',"MarkerEdgeColor",colormat(iiii),"MarkerFaceColor",colormat(iiii));
plot([x1(:)';x2(:)'], [y1(:)';y2(:)'], 'k-','Color',colormat(iiii),'LineWidth',0.25 )
xticks([1 2])
xticklabels({'Pre-LY354740','Post-LY354740'})
ylim([0 0.6])
end

%DmGluRA-RNAi Is
DsSortInd = 1:(length(DmGluRA_Is_Sorted)/5):(length(DmGluRA_Is_Sorted)+1);
DsSortInd = floor(DsSortInd);
DSSortFifths =[];
for ii = 1:5
    start_ind = DsSortInd(ii);
    end_ind = DsSortInd(ii+1)-1;
    DSSortFifths{ii} = DmGluRA_Is_Sorted(start_ind:end_ind,:)
end
maxlengthDS = []
for nnnn = 1:5
maxlengthDS(nnnn) = length(DSSortFifths{nnnn});
end
maxmaxind = max(maxlengthDS);
labelcolumn1 = ones(maxmaxind,2);
labelcolumn1(:,2) = labelcolumn1(:,1)+labelcolumn1(:,2);

dssize = size(DSSortFifths);
colormat = ["#1F2041", "#4B3F72","#FFC857","#119DA4","#19647E"]
tiledlayout(1,5);
for iiii = 1:5
nexttile
hold on
xind = length(DSSortFifths{iiii}(:,1));
x1 = labelcolumn1(1:xind,1);
y1 = DSSortFifths{iiii}(:,1);
x2 = labelcolumn1(1:xind,2);
y2 = DSSortFifths{iiii}(:,2);
scatter(labelcolumn1(1:xind,1),DSSortFifths{iiii}(:,1),2,'filled',"MarkerEdgeColor",colormat(iiii),"MarkerFaceColor",colormat(iiii));
scatter(labelcolumn1(1:xind,2),DSSortFifths{iiii}(:,2),2,'filled',"MarkerEdgeColor",colormat(iiii),"MarkerFaceColor",colormat(iiii));
plot([x1(:)';x2(:)'], [y1(:)';y2(:)'], 'k-','Color',colormat(iiii),'LineWidth',0.25 )
xticks([1 2])
xticklabels({'Pre-LY354740','Post-LY354740'})
ylim([0 0.6])
end

