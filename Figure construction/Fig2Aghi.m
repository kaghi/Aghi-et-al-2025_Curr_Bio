

%Substitute in the modality you want as a raster
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

scatter(xPoints,yPoints,'.','black');
xlim([1 100])
ylim([0 165])
xlabel('Stimulus Number')
ylabel('Active Zone ID')
title('5 Hz Train 1')
ax = gca;
xticks(trials);