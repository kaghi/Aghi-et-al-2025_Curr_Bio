AnimalNo = 7;
coordim = size(ExperimentSet_Reduced(3).All_Pair_Structures(AnimalNo).QuaSOR_STORM_Pair_Structure_Sorted_Verified);
StormAZCoords2 = cell(coordim(2),2);
for nn = 1:(coordim(2));
StormAZCoords2{nn,1} = ExperimentSet_Reduced(3).All_Pair_Structures(AnimalNo).QuaSOR_STORM_Pair_Structure_Sorted_Verified(nn).STORM_Coord_Orig;
StormAZCoords2{nn,2} = ExperimentSet_Reduced(3).All_Pair_Structures(AnimalNo).QuaSOR_STORM_Pair_Structure_Sorted_Verified(nn).STORM_Matched_AZ_BoutonType;
end
AZCoordX = zeros(coordim(2),1);
AZCoordY = zeros(coordim(2),1);
BoutType = zeros(coordim(2),1);
for ggg = 1:(coordim(2));
    AZCoordX(ggg,1) = StormAZCoords2{ggg,1}(1);
    AZCoordY(ggg,1) = StormAZCoords2{ggg,1}(2);
    TempBout = StormAZCoords2{ggg,2};
    BoutType(ggg,1) = TempBout;
end

BoutIbInd = find(BoutType(:,1)==1);
AZCoordXIb = AZCoordX(1:length(BoutIbInd),1);
AZCoordYIb = AZCoordY(1:length(BoutIbInd),1);

Fig2SampleCoordsX = AZCoordXIb;
Fig2SampleCoordsY = AZCoordYIb;
Fig2SampleCoordsLFPr = ExperimentSet_Reduced(3).Grouped_Data_Reduced(AnimalNo).Verified_Quantifications.BoutonSorted(1).All_QuaSOR_Data(4).Recording(1).Evoked_Pr;
Fig2SampleCoordsLFPr =Fig2SampleCoordsLFPr.';
% x0=10;
% y0=10;
% width=269*4 ;
% height=696; 
c = Fig2SampleCoordsLFPr;
scatter(Fig2SampleCoordsY,Fig2SampleCoordsX,30,c,'filled');
% set(gcf,'position',[x0,y0,width,height])
% make data labels:
%text(Fig2SampleCoordsY,Fig2SampleCoordsX,sprintfc(' %d',1:numel(Fig2SampleCoordsX)))

ax = gca;
ax.TickDir = 'out';
ax.TickLength = [0 0];
ax.XTickLabel = [];
ax.YTickLabel = [];
box on
grid off
c = colorbar;
w = c.LineWidth;
c.LineWidth = 1.5;
c.FontSize = 14;
colormap hot
y1 = ylabel(c,'Probability of release (Pr)','FontSize',16,'Rotation',270);
yl.Position(1) = min(xlim(c)) + 1;
yl.VerticalAlignment = 'bottom';
ax.PositionConstraint = 'position';
c.Position(1) = c.Position(1) + 0.05;
c.Position(1) = -0.05;
