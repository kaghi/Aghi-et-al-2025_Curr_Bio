xcoord = ExperimentSet_Reduced(3).Grouped_Data_Reduced(1).Verified_Quantifications.All_QuaSOR_Data(1).Recording.AllEvents(1).QuaSOR_All_Location_Coords_Stim_Matched(:,1);
ycoord = ExperimentSet_Reduced(3).Grouped_Data_Reduced(1).Verified_Quantifications.All_QuaSOR_Data(1).Recording.AllEvents(1).QuaSOR_All_Location_Coords_Stim_Matched(:,2);
xcoordav = mean(xcoord);
ycoordav = mean(ycoord);
scatter(xcoord,ycoord)
hold on 
scatter(xcoordav, ycoordav, 'filled')


%Train 1
AnimalNo = 3
DimNo = size(Grouped_Data_Reduced(AnimalNo).Verified_Quantifications.BoutonSorted(2).All_QuaSOR_Data(4).Recording(1).AllEvents)
RasterVals = cell(DimNo(2),1);

DimNoLF = size(Grouped_Data_Reduced(AnimalNo).Verified_Quantifications.BoutonSorted(2).All_QuaSOR_Data(1).Recording(1).AllEvents)
RasterValsKLF = cell(DimNoLF(2),1);


DimNoHF = size(Grouped_Data_Reduced(AnimalNo).Verified_Quantifications.BoutonSorted(2).All_QuaSOR_Data(4).Recording(1).AllEvents)
RasterValsKHF = cell(DimNo(2),1);

%LF Train
for i = 1:DimNo(2)
cellz = Grouped_Data_Reduced(AnimalNo).Verified_Quantifications.BoutonSorted(2).All_QuaSOR_Data(1).Recording(1).AllEvents(i).QuaSOR_All_Location_Coords_Matched;
if isempty(cellz)
RasterValsKLF{i,1} = 0
else
convcell = Grouped_Data_Reduced(AnimalNo).Verified_Quantifications.BoutonSorted(2).All_QuaSOR_Data(1).Recording(1).AllEvents(i).QuaSOR_All_Location_Coords_Matched(:,3);    
RasterValsKLF{i,1} = convcell.' ;
end
end
 plotSpikeRaster(RasterValsKLF,'PlotType','vertline')
 
 
 LowerB = 168
 UpperB = 198
%HF Train
 for i = 1:DimNoHF(2)
cellz = Grouped_Data_Reduced(AnimalNo).Verified_Quantifications.All_QuaSOR_Data(4).Recording(1).AllEvents(i).QuaSOR_All_Location_Coords_Matched;
if isempty(cellz)
RasterValsKHF{i,1} = 0
else
convcell = Grouped_Data_Reduced(AnimalNo).Verified_Quantifications.All_QuaSOR_Data(4).Recording(1).AllEvents(i).QuaSOR_All_Location_Coords_Matched(:,3);    
RasterValsKHF{i,1} = convcell.' ;
end
 end

 %Compute Pr for all sites in LF vs HF stimulation%
PrLF = zeros(length(RasterValsKLF),1);
for i = 1:DimNoLF(2)
    if RasterValsKLF{i} == 0
        PrLF(i) = 0;
    else
    PrLF(i) = length(RasterValsKLF{i});
    PrLF(i) = PrLF(i)/200;
    end
end
PrHF = zeros(length(RasterValsKHF),1);
for i = 1:DimNoLF(2)
    if RasterValsKHF{i} == 0
        PrHF(i) = 0;
    else
    PrHF(i) = length(RasterValsKHF{i});
    PrHF(i) = PrHF(i)/100;
    end
end

histogram(PrLF, 'FaceColor', [0 0.4470 0.7410], 'EdgeColor', [0  0 0]);
hold on
histogram(PrHF, 'FaceColor',[0.8500 0.3250 0.0980], 'EdgeColor', [0  0 0] );
legend('0.2Hz Stimulation', '5 Hz Stimulation Train 1')
xlabel('Pr')
ylabel('count')

RasterValsMostDistalLF = RasterValsKLF(LowerB:UpperB,1);
RasterValsMostDistalHF = RasterValsKHF(LowerB:UpperB,1);

 plotSpikeRaster(RasterValsMostDistalLF,'PlotType','vertline')
  plotSpikeRaster(RasterValsMostDistalHF,'PlotType','vertline','LineFormat', LineFormat)

%Compute Prs to Compare LF and HF stimulation of the same sites%
MostDistalPrLF = zeros(length(RasterValsMostDistalLF),1);
for i = LowerB: UpperB
    if RasterValsMostDistalLF{i} == 0
        MostDistalPrLF(i) = 0;
    else
    MostDistalPrLF(i) = length(RasterValsMostDistalLF{i});
    MostDistalPrLF(i) = MostDistalPrLF(i)/100;
    end
end
MostDistalPrHF = zeros(length(RasterValsMostDistalHF),1);
for i = LowerB: UpperB
    if RasterValsMostDistalHF{i} == 0
        MostDistalPrHF(i) = 0;
    else
    MostDistalPrHF(i) = length(RasterValsMostDistalHF{i});
    MostDistalPrHF(i) = MostDistalPrHF(i)/100;
    end
end


histogram(MostDistalPrLF,'FaceColor', [0 0.4470 0.7410], 'EdgeColor', [0  0 0]);
hold on
histogram(MostDistalPrHF,'FaceColor',[0.8500 0.3250 0.0980], 'EdgeColor', [0  0 0]);
legend('0.2Hz Stimulation', '5 Hz Stimulation Train 1')
xlabel('Pr')
ylabel('count')


MostDistalPrLF()
length(RasterValsMostDistalLF{})
%Extract coordinates of most distal bouton
 RasterCoordsMostDistalLF = cell(UpperB,1);
 for i = LowerB:UpperB
cellz2 = Grouped_Data_Reduced(AnimalNo).Verified_Quantifications.All_QuaSOR_Data(1).AllEvents(i).QuaSOR_All_Location_Coords_Matched;
if isempty(cellz2)
RasterCoordsMostDistalLF{i,1} = 0
else
avgx = mean(cellz2(:,1));
avgy = mean(cellz2(:,2));
cellcoord = zeros(1,2);
cellcoord(1) = avgx;
cellcoord(2) = avgy;
RasterCoordsMostDistalLF{i,1} = cellcoord ;
end
 end
 
 %Extract coordinates of most distal bouton for HF 
  RasterCoordsMostDistalHF = cell(UpperB,1);
  for i = LowerB:UpperB
cellz2 = Grouped_Data_Reduced(AnimalNo).Verified_Quantifications.All_QuaSOR_Data(4).AllEvents(i).QuaSOR_All_Location_Coords_Matched;
if isempty(cellz2)
RasterCoordsMostDistalHF{i,1} = 0
else
avgx = mean(cellz2(:,1));
avgy = mean(cellz2(:,2));
cellcoord = zeros(1,2);
cellcoord(1) = avgx;
cellcoord(2) = avgy;
RasterCoordsMostDistalHF{i,1} = cellcoord ;
end
 end
 
%HF values for the same sites in Train 1 
for i = 1:DimNoLF(2)
cellz = Grouped_Data_Reduced(AnimalNo).Verified_Quantifications.All_QuaSOR_Data(1).Recording(1).AllEvents(i).QuaSOR_All_Location_Coords_Matched;
if isempty(cellz)
RasterValsKLF{i,1} = 0
else
convcell = Grouped_Data_Reduced(AnimalNo).Verified_Quantifications.All_QuaSOR_Data(1).Recording(1).AllEvents(i).QuaSOR_All_Location_Coords_Matched(:,3);    
RasterValsKLF{i,1} = convcell.' ;
end
end

 
 
Ind = (UpperB - LowerB)+ 1 
%LF Coord Extraction
xvalsLF = zeros(Ind,1);
yvalsLF = zeros(Ind,1);
for nn = LowerB:UpperB
if RasterCoordsMostDistalLF{nn,1} == 0 
    xvalsLF(nn,1)= 0;
    yvalsLF(nn,1)= 0;
else
xvalsLF(nn,1) = RasterCoordsMostDistalLF{nn,1}(1,1);
yvalsLF(nn,1) = RasterCoordsMostDistalLF{nn,1}(1,2);
end

end

%HF Coord Extraction
xvalsHF = zeros(Ind,1);
yvalsHF = zeros(Ind,1);
for nn = LowerB:UpperB
if RasterCoordsMostDistalHF{nn,1} == 0 
    xvalsHF(nn,1)= 0;
    yvalsHF(nn,1)= 0;
else
xvalsHF(nn,1) = RasterCoordsMostDistalHF{nn,1}(1,1);
yvalsHF(nn,1) = RasterCoordsMostDistalHF{nn,1}(1,2);
end

end


truexvalsLF = xvalsLF(LowerB:UpperB,1);
trueyvalsLF = yvalsLF(LowerB:UpperB,1);

truexvalsHF = xvalsHF(LowerB:UpperB,1);
trueyvalsHF = yvalsHF(LowerB:UpperB,1);

a = [1:Ind]';
b = num2str(a); 
c = cellstr(b);
dx = 0.5; 
dy = 0.5; % displacement so the text does not overlay the data points
scatter(truexvalsLF,trueyvalsLF,'filled')
text(truexvalsLF+dx,trueyvalsLF+dy,c);
hold on
scatter(truexvalsHF,trueyvalsHF,'filled')
text(truexvalsHF+dx,trueyvalsHF+dy,c);

 %Train 1
for i = 1:DimNo(2)
cellz = Grouped_Data_Reduced(AnimalNo).Verified_Quantifications.BoutonSorted(1).All_QuaSOR_Data(4).Recording(1).AllEvents(i).QuaSOR_All_Location_Coords_Stim_Matched;
if isempty(cellz)
RasterVals{i,1} = 0
else
convcell = Grouped_Data_Reduced(AnimalNo).Verified_Quantifications.BoutonSorted(1).All_QuaSOR_Data(4).Recording(1).AllEvents(i).QuaSOR_All_Location_Coords_Stim_Matched(:,3);    
RasterVals{i,1} = convcell.' ;
end
end

%Train 2  
for i = 1:DimNo(2)
cellz = Grouped_Data_Reduced(AnimalNo).Verified_Quantifications.All_QuaSOR_Data(4).Recording(2).AllEvents(i).QuaSOR_All_Location_Coords_Stim_Matched;
if isempty(cellz)
RasterVals{i,1} = 0
else
convcell = Grouped_Data_Reduced(AnimalNo).Verified_Quantifications.All_QuaSOR_Data(4).Recording(2).AllEvents(i).QuaSOR_All_Location_Coords_Stim_Matched(:,3);    
RasterVals{i,1} = convcell.' ;
end
end

%Train 3
for i = 1:DimNo(2)
cellz = Grouped_Data_Reduced(AnimalNo).Verified_Quantifications.All_QuaSOR_Data(4).Recording(3).AllEvents(i).QuaSOR_All_Location_Coords_Stim_Matched;
if isempty(cellz)
RasterVals{i,1} = 0
else
convcell = Grouped_Data_Reduced(AnimalNo).Verified_Quantifications.All_QuaSOR_Data(4).Recording(3).AllEvents(i).QuaSOR_All_Location_Coords_Stim_Matched(:,3);    
RasterVals{i,1} = convcell.' ;
end
end

%Train 4
for i = 1:DimNo(2)
cellz = Grouped_Data_Reduced(AnimalNo).Verified_Quantifications.All_QuaSOR_Data(4).Recording(4).AllEvents(i).QuaSOR_All_Location_Coords_Stim_Matched;
if isempty(cellz)
RasterVals{i,1} = 0
else
convcell = Grouped_Data_Reduced(AnimalNo).Verified_Quantifications.All_QuaSOR_Data(4).Recording(4).AllEvents(i).QuaSOR_All_Location_Coords_Stim_Matched(:,3);    
RasterVals{i,1} = convcell.' ;
end
end

%Train 5 
for i = 1:DimNo(2)
cellz = Grouped_Data_Reduced(AnimalNo).Verified_Quantifications.All_QuaSOR_Data(4).Recording(5).AllEvents(i).QuaSOR_All_Location_Coords_Stim_Matched;
if isempty(cellz)
RasterVals{i,1} = 0
else
convcell = Grouped_Data_Reduced(AnimalNo).Verified_Quantifications.All_QuaSOR_Data(4).Recording(5).AllEvents(i).QuaSOR_All_Location_Coords_Stim_Matched(:,3);    
RasterVals{i,1} = convcell.' ;
end
end

% for ii = 1:151
% times = RasterVals{ii,1};
% rasterplot(times,151,100);
% hold on
% end
C = 1:DimNo(2);
LineFormat.LineWidth = 2
LineFormat.Color = [0.3 0.3 0.3]
LineFormat.LineStyle = "-"
 plotSpikeRaster(RasterVals,'PlotType','vertline','LineFormat', LineFormat)
 
 RasterCoords = cell(DimNo(2),1);
 for i = 1:DimNo(2)
cellz2 = Grouped_Data_Reduced(AnimalNo).Verified_Quantifications.All_QuaSOR_Data(4).AllEvents(i).QuaSOR_All_Location_Coords_Matched;
if isempty(cellz2)
RasterCoords{i,1} = 0
else
avgx = mean(cellz2(:,1));
avgy = mean(cellz2(:,2));
cellcoord = zeros(1,2);
cellcoord(1) = avgx;
cellcoord(2) = avgy;
RasterCoords{i,1} = cellcoord ;
end
 end

xvals = zeros(DimNo(2),1);
yvals = zeros(DimNo(2),1);
for nn = 1:DimNo(2)
if RasterCoords{nn,1} == 0 
    xvals(nn,1)= 0;
    yvals(nn,1)= 0;
else
xvals(nn,1) = RasterCoords{nn,1}(1,1);
yvals(nn,1) = RasterCoords{nn,1}(1,2);
end

end

C = 1:DimNo(2);
scatter(xvals,yvals,40,C,'filled')
hAx = gca;
cbar = colorbar;
colormap(winter) %your original map...

%Let us calculate the time between successes for all identified sites
Time_diff = cell(length(RasterVals),1);
for ii = 1:length(RasterVals)
    
    