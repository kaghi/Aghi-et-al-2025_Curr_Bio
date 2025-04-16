%Load from Quasor_Data file
AreaBouton = nnz(AllBoutonsRegion);
NumberEvents = QuaSOR_Data.Total_Coords;

%Load from Analysis_Setup file
PixelConversion = ScaleBar.ScaleFactor;

ActualArea = AreaBouton*PixelConversion;

%QD
QuantalDensity = NumberEvents/ActualArea;
QuantalDensityPerStim = QuantalDensity/100;

%QD over time (for episodic data)
FrameNo = QuaSOR_Data.All_Location_Coords_byEpisodeNum(:,3);
IndCell = [];
Counter = 0;
ParsedInd = max(FrameNo);
for kk = 1:ParsedInd
    StartInd = kk;
    Frames = find(FrameNo == StartInd); 
    IndCell{kk} = Frames;
end

%QD over time (for 5 Hz  data)
FrameNo = QuaSOR_Data.All_Location_Coords_byOverallStim(:,3);
IndCell = [];
Counter = 0;
ParsedInd = max(FrameNo);
for kk = 1:ParsedInd
    StartInd = kk;
    Frames = find(FrameNo == StartInd); 
    IndCell{kk} = Frames;
end
% %QD over time (for streaming data), every 10 frames
% FrameNo = QuaSOR_Data.All_Location_Coords_byEpisodeFrame(:,3);
% IndCell = [];
% Counter = 0;
% ParsedInd = length(FrameNo)/10;
% for kk = 1:ParsedInd
%     StartInd = Counter;
%     EndInd = Counter+10;
%     Frames = find(FrameNo > StartInd & FrameNo<= EndInd); 
%     IndCell{kk} = Frames;
%     Counter = Counter + 10;
% end
% %QD over time (for streaming data), every 5 frames
% FrameNo = QuaSOR_Data.All_Location_Coords_byEpisodeFrame(:,3);
% IndCell = [];
% Counter = 0;
% ParsedInd = length(FrameNo)/5;
% for kk = 1:ParsedInd
%     StartInd = Counter;
%     EndInd = Counter+10;
%     Frames = find(FrameNo > StartInd & FrameNo<= EndInd); 
%     IndCell{kk} = Frames;
%     Counter = Counter + 10;
% end

%plot results
QD_Cell = [];
True_QD_Cell = [];
for lll = 1:ParsedInd
    QD_Cell{lll} = length(IndCell{lll});
    True_QD_Cell{lll} = QD_Cell{lll}/ActualArea;
end

for rrrr = 1:ParsedInd
bar(rrrr,True_QD_Cell{rrrr},'FaceColor',[0 .5 .5]);
hold on
end
ylabel('Quantal Density')
%for episodic
xlabel('Episode Number')
%for 5 Hz 
xlabel('Trial number')

%for streaming
xlabel('Time (0.5s)')

QD_Mat = cell2mat(QD_Cell);
QD_Mat = QD_Mat.';
True_QD_Mat = cell2mat(True_QD_Cell);
True_QD_Mat = True_QD_Mat.';

clear all
close all
