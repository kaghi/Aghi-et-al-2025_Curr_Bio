%Requires plotSpikeRaster.m by:
%% AUTHOR    : Jeffrey Chiou
%% $DATE     : 07-Feb-2014 12:15:47 $
%% $Revision : 1.2 $
%% DEVELOPED : 8.1.0.604 (R2013a)
%% FILENAME  : plotSpikeRaster.m


close
[baseName, folder] = uigetfile();
fullFileName = fullfile(folder, baseName)
T = xlsread(baseName);
R = size(T);
endl = R(2);
%for Spont (comment out if not using modality)
% Tarray = T(:,7:endl);
% Tarray_t = Tarray.';
% rastdata_t = logical(Tarray_t);
% %for 0.2 Hz
Tarray = T(:,6:100);
Tarray_t = Tarray.';
rastdata_t = logical(Tarray_t);
rastdata = logical(Tarray);
rastdata_t = transpose(rastdata);
LineFormat.LineWidth = 0.5;
[xPoints yPoints]= plotSpikeRaster(rastdata_t,'PlotType','vertline');
% xlabel('Frame Number')
% ylabel('AZ')
% title('Spontaneous Release Events')
xlabel('Stimulus Number')
ylabel('AZ')
title('0.2 Hz Release Events')
newStr = erase(baseName,".xlsx")
saveas(gcf, fullfile(folder, newStr),'emf');
