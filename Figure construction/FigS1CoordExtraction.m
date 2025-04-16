coordim = size(ExperimentSet_Reduced(3).All_Pair_Structures(8).QuaSOR_STORM_Pair_Structure_Sorted_Verified);
StormAZCoords2 = cell(coordim(2),2);
for nn = 1:(coordim(2));
StormAZCoords2{nn,1} = ExperimentSet_Reduced(3).All_Pair_Structures(8).QuaSOR_STORM_Pair_Structure_Sorted_Verified(nn).STORM_Coord_Orig;
StormAZCoords2{nn,2} = ExperimentSet_Reduced(3).All_Pair_Structures(8).QuaSOR_STORM_Pair_Structure_Sorted_Verified(nn).STORM_Matched_AZ_BoutonType;
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

%%Full NMJ
scatter(AZCoordXIb,AZCoordYIb,[],'filled');

axis equal
xlim([0 155])
ylim([0 260])

%Pick 30 AZs
scatter(AZCoordXIb(1:30),AZCoordYIb(1:30),[],'filled');
axis equal
xlim([60 155])
ylim([0 60])
dx = 0.2;
dy = 0.2;
text(AZCoordXIb(1:30)+dx,AZCoordYIb(1:30)+dy,sprintfc(' %d',1:30),'FontSize', 14)

%%Extract spikes for specific 30 AZ
spiketimes = Quasor_Refined_AZ_Stim_Matched_Events.LF;
spiketimes = Quasor_Refined_AZ_Stim_Matched_Events.HF1;
spiketimes = Quasor_Refined_AZ_Stim_Matched_Events.HF2;

ibonly = cell(164,1);
trials = [0:5:100];

for nn = 1:164;
  if isempty(spiketimes{nn})
  ibonly{nn} = 0;
  else
  ibonly{nn} = spiketimes{nn};
  end
  ibonly{nn} = ibonly{nn}.'
end

[xPoints, yPoints] = plotSpikeRaster(ibonly);
for jj = 1:length(xPoints)
   if xPoints(jj,1) == 0
   xPoints(jj,1) = NaN;
   else
   xPoints(jj,1)  = xPoints(jj,1);
   end
end

scatter(xPoints,yPoints,'.','black');
xlim([0 100])
ylim([0 165])
xlabel('Stimulus Number')
ylabel('Active Zone ID')
title('0.2 Hz')
ax = gca;
xticks(trials);

    
