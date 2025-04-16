


%Train 1
AnimalNo = 7

DimNoLF = size(ExperimentSet_Reduced(3).Grouped_Data_Reduced(AnimalNo).Verified_Quantifications.BoutonSorted(1).All_QuaSOR_Data(1).Recording(1).AllEvents)
RasterValsKLF = cell(DimNoLF(2),1);

DimNoHF = size(ExperimentSet_Reduced(3).Grouped_Data_Reduced(AnimalNo).Verified_Quantifications.BoutonSorted(1).All_QuaSOR_Data(4).Recording(1).AllEvents)
RasterValsKHF = cell(DimNoHF(2),1);

%LF Train
for i = 1:DimNoLF(2)
cellz = ExperimentSet_Reduced(3).Grouped_Data_Reduced(AnimalNo).Verified_Quantifications.BoutonSorted(1).All_QuaSOR_Data(1).Recording(1).AllEvents(i).QuaSOR_All_Location_Coords_Stim_Matched;
if isempty(cellz)
RasterValsKLF{i,1} = 0
else
convcell = ExperimentSet_Reduced(3).Grouped_Data_Reduced(AnimalNo).Verified_Quantifications.BoutonSorted(1).All_QuaSOR_Data(1).Recording(1).AllEvents(i).QuaSOR_All_Location_Coords_Stim_Matched(:,3);    
RasterValsKLF{i,1} = convcell.' ;
end
end
 plotSpikeRaster(RasterValsKLF,'PlotType','vertline')
 xlabel('Stimulus Number')
ylabel('Active Zone ID')
title('0.2 Hz')

 
 
%HF Train 1
for i = 1:DimNoHF(2)
cellz = ExperimentSet_Reduced(3).Grouped_Data_Reduced(AnimalNo).Verified_Quantifications.BoutonSorted(1).All_QuaSOR_Data(4).Recording(1).AllEvents(i).QuaSOR_All_Location_Coords_Stim_Matched;
if isempty(cellz)
RasterValsKHF{i,1} = 0
else
convcell = ExperimentSet_Reduced(3).Grouped_Data_Reduced(AnimalNo).Verified_Quantifications.BoutonSorted(1).All_QuaSOR_Data(4).Recording(1).AllEvents(i).QuaSOR_All_Location_Coords_Stim_Matched(:,3);    
RasterValsKHF{i,1} = convcell.' ;
end
end
 plotSpikeRaster(RasterValsKHF,'PlotType','vertline')
 xlabel('Stimulus Number')
ylabel('Active Zone ID')
title('5 Hz (Train 1)')