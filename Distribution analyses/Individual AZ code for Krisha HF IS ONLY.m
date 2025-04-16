        
%Collecting trace data for individual AZs 
%Some notes: 
%Modality_Merge.BoutonArray(1).QuaSOR_Auto_AZ.QuaSOR_AutoQuant.AZ_Struct(1).DeltaFF0_MeanTrace--> trace for all the frames (across all episodes) 
%Modality_Merge.BoutonArray(1).QuaSOR_Auto_AZ.QuaSOR_AutoQuant.Episode(1).AZ_Struct--> trace for every single episode across all the "Auto AZs" not sure what the
%difference between DeltaFF0 mean and max is supposed to be [Coord,Coord_Round, PixelList, Coord_Orig, Mask (empty), DeltaFF0_MeanTrace, DeltaFF0_MaxTrace] 

%% START HERE: Part 1: 
%GOAL: make an AZ structure similar to the Auto_AZ one [Coord, Coord_Round, PixelList] 
%QuaSOR_Airy_Pair_Structure(1).QuaSOR_Prelim_Quant_ROI_RegionProps has the centroid data plus the pixel list 
%(note that QuaSOR_Airy_Pair Structure_Refined does not have the centroid with the ROI region) 
 
%TO LOAD: (1) QuaSOR_Airy_Pair_Structure, (2) QuaSOR_Parameters. (1)
%Located in Modality_Merge, (2) in _Quasor_Data.mat. 

AZ_Temp_Structure = ([]) %creates an empty structure 
for i = 1:size(QuaSOR_Airy_Pair_Structure,2)
    AZ_Temp_Structure(i).Coord = QuaSOR_Airy_Pair_Structure(i).QuaSOR_Prelim_Quant_ROI_RegionProps.Centroid
    AZ_Temp_Structure(i).Coord_Round = round(AZ_Temp_Structure(i).Coord)
    AZ_Temp_Structure(i).PixelList = QuaSOR_Airy_Pair_Structure(i).QuaSOR_Prelim_Quant_ROI_RegionProps.PixelList
    AZ_Temp_Structure(i).Coord_Orig = round(AZ_Temp_Structure(i).Coord_Round/QuaSOR_Parameters.UpScaling.QuaSOR_UpScaleFactor)
    AZ_Temp_Structure(i).Bouton = QuaSOR_Airy_Pair_Structure(i).Bouton
end
%make a directory for each AZ
selpath = uigetdir;
addpath(genpath('K:\QuaSOR Airy Data'));
cd(selpath)
mkdir('AZs')
%% Part 2: Adding remapping coordinates to the AZ_Temp_Structure

%TO LOAD: BoutonMerge.BoutonArray (make sure that it has the
%XCoordinate_ReMapping field. Located _QuaSOR_Data.mat. 
Ib_X_Coord_ReMap = BoutonMerge.BoutonArray(1).XCoordinate_ReMapping
Ib_Y_Coord_ReMap = BoutonMerge.BoutonArray(1).YCoordinate_ReMapping
Is_X_Coord_ReMap = BoutonMerge.BoutonArray(2).XCoordinate_ReMapping
Is_Y_Coord_ReMap = BoutonMerge.BoutonArray(2).YCoordinate_ReMapping

for i = 1:size(QuaSOR_Airy_Pair_Structure,2)
    
    if AZ_Temp_Structure(i).Bouton == 1 
        AZ_Temp_Structure(i).Coord_Orig_Remapped(1) = ...
            AZ_Temp_Structure(i).Coord_Orig(1) - Ib_X_Coord_ReMap
        AZ_Temp_Structure(i).Coord_Orig_Remapped(2) = ...
            AZ_Temp_Structure(i).Coord_Orig(2) - Ib_Y_Coord_ReMap
    else
        AZ_Temp_Structure(i).Coord_Orig_Remapped(1) = ...
            AZ_Temp_Structure(i).Coord_Orig(1) - Is_X_Coord_ReMap
        AZ_Temp_Structure(i).Coord_Orig_Remapped(2) = ...
            AZ_Temp_Structure(i).Coord_Orig(2) - Is_Y_Coord_ReMap
    end
end
%% Part 3A: Adding traces for Ib (loading variables for Ib ONLY right now). Skip if no Ib
%TO LOAD: ImageAnalysisSetup.mat file for Ib and EpisodeStructure file for Ib 


Image_Width = Image_Width %depends on how you crop it. comes from image analysis setup 
Image_Height = Image_Height %depends on how you crop it 
[DeltaFF0_Col DeltaFF0_Row] = meshgrid(1:Image_Width,...
    1:Image_Height);
ZerosImage=zeros(Image_Height,Image_Width,'logical');
main_bouton = 1 %1 = Ib 

for EpisodeNumber_Load=1:ImagingInfo.NumEpisodes
    fprintf(['Extracting Episode ',num2str(EpisodeNumber_Load),' Traces...'])

    if SplitEpisodeFiles %for HF data 
        FileSuffix=['_DeltaFData_Ep_',num2str(EpisodeNumber_Load),'.mat'];
        fprintf(['Loading: ',FileSuffix,'...'])
        load([CurrentScratchDir,StackSaveName,FileSuffix],'EpisodeStruct')
        EpisodeNumber=1;
    else
        EpisodeNumber=EpisodeNumber_Load;
   end
    TempDeltaFF0=EpisodeStruct(EpisodeNumber).ImageArrayReg_Episode_DeltaFF0; 
    %ImageArrayReg_EpisodeDFF0 (contains all the frames in the cropped image for every single episode).
    %For example, this is imagesc(EpisodeStruct(1).ImageArrayReg_Episode_DeltaFF0(:,:,4).
    %First Episode, 4th frame 
    %%%%%%%%%%%%%%%%%%%%%%%%

    %QuaSOR AZ Centered Traces
    fprintf('by AZ Pooled Mod...')
    AZ_Struct=AZ_Temp_Structure;
    TempMask=zeros(Image_Height,Image_Width,'uint16');
    QuaSOR_Auto_AZ.QuaSOR_AutoQuant.AZTrace.RegionRadius_px = 2
    QuaSOR_Auto_AZ.QuaSOR_AutoQuant.AZTrace.RegionRadius_nm = 350
    %Filling out the AZ Struct with regions (temp masks centered on
    %the structure verified AZs (not the Auto AZs) 
    for AZ=1:length(AZ_Struct)
        if AZ_Struct(AZ).Bouton == main_bouton %1 is Ib 2 is Is. Make sure to change for every single one 

        %AZ_Struct(AZ).Coord_Orig=AZ_Struct(AZ).Coord;
        %AZ_Struct(AZ).Coord_Orig=round(AZ_Struct(AZ).Coord_Orig/QuaSOR_Parameters.UpScaling.QuaSOR_UpScaleFactor);
            AZ_Struct(AZ).Mask =    (DeltaFF0_Row - AZ_Struct(AZ).Coord_Orig_Remapped(2)).^2 + ...
                (DeltaFF0_Col - AZ_Struct(AZ).Coord_Orig_Remapped(1)).^2 <= ...
                QuaSOR_Auto_AZ.QuaSOR_AutoQuant.AZTrace.RegionRadius_px.^2; %region radius is 2 pixels^2 = 4 pixels 
            AZ_Struct(AZ).Mask=logical(AZ_Struct(AZ).Mask);
            TempMask=TempMask+uint16(AZ_Struct(AZ).Mask);
            AZ_Struct(AZ).DeltaFF0_MeanTrace=[]; 
            AZ_Struct(AZ).DeltaFF0_MaxTrace=[];
        end

    end


      
    for AZ=1:length(AZ_Struct) %used to be a parfor loop 
        TempMaxVals=[];
        TempMeanTrace=[];
        if AZ_Struct(AZ).Bouton == main_bouton
            for f1=1:size(TempDeltaFF0,3)
                TempImage=TempDeltaFF0(:,:,f1);
                TempVals=TempImage(AZ_Struct(AZ).Mask);
                TempMaxVals(f1)=max(TempVals(:));
                TempMeanTrace(f1)=nanmean(TempVals(:));
            end
            AZ_Struct(AZ).DeltaFF0_MeanTrace=TempMeanTrace;
            AZ_Struct(AZ).DeltaFF0_MaxTrace=TempMaxVals;
        end

    end

   QuaSOR_Matched_AZ_Structure.Ib.Episode(EpisodeNumber_Load).AZ_Struct = AZ_Struct %creates QuaSOR_Matched_AZ_Structure 


end
%% Part 3B: Same as above except now we're adding the Is traces 
% TO LOAD: Only reload the ImageAnalysis file and the EpisodeStruc for Is.
% Don't clear anything else!

Image_Width = Image_Width %depends on how you crop it 
Image_Height = Image_Height %depends on how you crop it 
[DeltaFF0_Col DeltaFF0_Row] = meshgrid(1:Image_Width,...
    1:Image_Height);
ZerosImage=zeros(Image_Height,Image_Width,'logical');
main_bouton = 2 %2 when investigating Is 
other_bouton = 1 %1 when investigating Is

for EpisodeNumber_Load=1:ImagingInfo.NumEpisodes
    fprintf(['Extracting Episode ',num2str(EpisodeNumber_Load),' Traces...'])

    if SplitEpisodeFiles %for HF data 
        FileSuffix=['_DeltaFData_Ep_',num2str(EpisodeNumber_Load),'.mat'];
        fprintf(['Loading: ',FileSuffix,'...'])
        load([CurrentScratchDir,StackSaveName,FileSuffix],'EpisodeStruct')
        EpisodeNumber=1;
    else
        EpisodeNumber=EpisodeNumber_Load;
   end
    TempDeltaFF0=EpisodeStruct(EpisodeNumber).ImageArrayReg_Episode_DeltaFF0; 

    %QuaSOR AZ Centered Traces
    fprintf('by AZ Pooled Mod...')
    AZ_Struct=AZ_Temp_Structure;
    TempMask=zeros(Image_Height,Image_Width,'uint16');
    QuaSOR_Auto_AZ.QuaSOR_AutoQuant.AZTrace.RegionRadius_px = 2
    QuaSOR_Auto_AZ.QuaSOR_AutoQuant.AZTrace.RegionRadius_nm = 350

    for AZ=1:length(AZ_Struct)
        if AZ_Struct(AZ).Bouton == main_bouton %1 is Ib 2 is Is. Make sure to change for every single one 
            AZ_Struct(AZ).Mask =    (DeltaFF0_Row - AZ_Struct(AZ).Coord_Orig_Remapped(2)).^2 + ...
                (DeltaFF0_Col - AZ_Struct(AZ).Coord_Orig_Remapped(1)).^2 <= ...
                QuaSOR_Auto_AZ.QuaSOR_AutoQuant.AZTrace.RegionRadius_px.^2; %region radius is 2 pixels^2 = 4 pixels 
            AZ_Struct(AZ).Mask=logical(AZ_Struct(AZ).Mask);
            TempMask=TempMask+uint16(AZ_Struct(AZ).Mask);
            AZ_Struct(AZ).DeltaFF0_MeanTrace=[]; 
            AZ_Struct(AZ).DeltaFF0_MaxTrace=[];
        end

    end

 
    for AZ=1:length(AZ_Struct) %used to be a parfor loop 
        TempMaxVals=[];
        TempMeanTrace=[];
        if AZ_Struct(AZ).Bouton == main_bouton
            for f1=1:size(TempDeltaFF0,3)
                TempImage=TempDeltaFF0(:,:,f1);
                TempVals=TempImage(AZ_Struct(AZ).Mask);
                TempMaxVals(f1)=max(TempVals(:));
                TempMeanTrace(f1)=nanmean(TempVals(:));
            end
            AZ_Struct(AZ).DeltaFF0_MeanTrace=TempMeanTrace;
            AZ_Struct(AZ).DeltaFF0_MaxTrace=TempMaxVals;
        else
            continue
        end

    end

   QuaSOR_Matched_AZ_Structure.Is.Episode(EpisodeNumber_Load).AZ_Struct = AZ_Struct %adding Is data to the QuaSOR_Matched_AZ_Structure

end
%% Part 4: Merging traces from Ib and Is within QuaSOR_Matched_AZ_Structure
for AZ=1:length(AZ_Struct)
    QuaSOR_Matched_AZ_Structure.Ib.AZ_Struct(AZ).DeltaFF0_MeanTrace=[];
    for EpisodeNumber_Load=1:ImagingInfo.NumEpisodes
        QuaSOR_Matched_AZ_Structure.Ib.AZ_Struct(AZ).DeltaFF0_MeanTrace=...
            horzcat(QuaSOR_Matched_AZ_Structure.Ib.AZ_Struct(AZ).DeltaFF0_MeanTrace,...
            QuaSOR_Matched_AZ_Structure.Ib.Episode(EpisodeNumber_Load).AZ_Struct(AZ).DeltaFF0_MeanTrace);
    end
end
for AZ=1:length(AZ_Struct)
    QuaSOR_Matched_AZ_Structure.Is.AZ_Struct(AZ).DeltaFF0_MeanTrace=[];
    for EpisodeNumber_Load=1:ImagingInfo.NumEpisodes
        QuaSOR_Matched_AZ_Structure.Is.AZ_Struct(AZ).DeltaFF0_MeanTrace=...
            horzcat(QuaSOR_Matched_AZ_Structure.Is.AZ_Struct(AZ).DeltaFF0_MeanTrace,...
            QuaSOR_Matched_AZ_Structure.Is.Episode(EpisodeNumber_Load).AZ_Struct(AZ).DeltaFF0_MeanTrace);
    end
end


QuaSOR_Matched_AZ_Structure.All_Boutons = []
QuaSOR_Matched_AZ_Structure.All_Boutons = [QuaSOR_Matched_AZ_Structure.Ib.AZ_Struct(1:212), QuaSOR_Matched_AZ_Structure.Is.AZ_Struct(213:end)] %change for every animal
        %%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1: length(QuaSOR_Matched_AZ_Structure.All_Boutons)
    QuaSOR_Matched_AZ_Structure.All_Boutons(i).AZ_ID = i
    QuaSOR_Matched_AZ_Structure.All_Boutons(i).Bouton = QuaSOR_Matched_AZ_Structure.Ib.Episode(1).AZ_Struct(i).Bouton
end

%% Part 5: Making the Airy_Alignment_Image_Color to plot the things in Part 6
% TO LOAD: Image_Data_Matched.mat file (usually found within the structural
% folder. Found in QuaSOR Airy Data/_AiryData.mat
%Just put PixelValueScalar = 1
Airy_Alignment_Contrast = 1.5 % usually 0.5
Airy_Alignment_Image=ImageData_Matched{AZ_Detection_Channel}.MaxProj;
Airy_Alignment_Image_Cont=Airy_Alignment_Image;
MaxValue=ImageData_Matched{AZ_Detection_Channel}.MaxValue;
MaxValue_Cont=MaxValue*Airy_Alignment_Contrast*PixelValueScalar;
Airy_Alignment_Image_Cont(Airy_Alignment_Image_Cont>MaxValue_Cont)=MaxValue_Cont;
Airy_Alignment_Image_Color=...
    single(grs2rgb(Airy_Alignment_Image_Cont,makeColorMap([0 0 0],...
    ColorDefinitions(ImageData_Matched{AZ_Detection_Channel}.Color),...
    ImageData_Matched{AZ_Detection_Channel}.MaxValue*Airy_Alignment_Contrast*PixelValueScalar)));
Airy_Alignment_Image_Color=single(ColorMasking(Airy_Alignment_Image_Color,~ImageData_MaxProj_StitchingMask_Matched_Crop,StitchingMask_Color));
figure, imshow(double(Airy_Alignment_Image_Color),[],'border','tight')


%% Part 6: Plotting. Will create Brp and traces on the same plot for each AZ of interest and automatically save them to a folder 
%Many of the subplots are commented out. Feel free to play around with this
%LOAD Airy_Alignment_Image_Color which is a carryover from Part 5 
%NOTE: This is a crude way of indexing manually but for AZ_nums it's 1: whatever your last index is for Ib AZs
%for Is it's first index of Is:length(QuaSOR_Matched_AZ_Structure.All_Boutons
for AZ_num = 243:243 %whatever AZ or index you want (see note above for Ib vs Is ranges)
    pre_traces = []
    NumEp = length(QuaSOR_Matched_AZ_Structure.Is.Episode);
    for episode_num = 1:NumEp %number of indices. variable should be in ImageAnalysis Setup which we loaded in earlier
        x = 1:10; %number of frames  
        %comment out depending if it's Ib or Is (ignore the "pre" prefix. that's for when I was splitting up the data 
       %pre_traces = [pre_traces; QuaSOR_Matched_AZ_Structure.Ib.Episode(episode_num).AZ_Struct(AZ_num).DeltaFF0_MeanTrace];
        pre_traces = [pre_traces QuaSOR_Matched_AZ_Structure.Is.Episode(episode_num).AZ_Struct(AZ_num).DeltaFF0_MeanTrace];
    end
    mean_pre_trace = mean(pre_traces);


    %Airy part
   
    x_center_airy = QuaSOR_Airy_Pair_Structure(AZ_num).Airy_Coord(1)
    y_center_airy = QuaSOR_Airy_Pair_Structure(AZ_num).Airy_Coord(2)
    constant = 20
    x_low_airy = x_center_airy - constant 
    x_high_airy = x_center_airy + constant
    y_low_airy = y_center_airy - constant 
    y_high_airy = y_center_airy + constant

    Airy_Alignment_Image_Color = insertMarker(Airy_Alignment_Image_Color,[QuaSOR_Airy_Pair_Structure(AZ_num).Airy_Coord_Match(1)...
        QuaSOR_Airy_Pair_Structure(AZ_num).Airy_Coord_Match(2)], 'size',1);

    AZ_Airy_Alignment_Image_Color_Brp = Airy_Alignment_Image_Color(y_low_airy:y_high_airy,x_low_airy:x_high_airy,:);

    figure, 
    subplot(3,4,1)
    imshow(AZ_Airy_Alignment_Image_Color_Brp)
    title('Airy Brp')
    txt = ['AZ #', num2str(AZ_num)]
    ylabel(txt)
    box on 

    
    subplot(3,4,[5 8])
    %subplot(2,2,[3,4])
    %for 0.2 Hz
        %x = [1:NumIndices*10]; 
%axis([1 NumIndices*10 -0.2 1.2])
%For 5 Hz
xlimmy = length(QuaSOR_Matched_AZ_Structure.Is.Episode(episode_num).AZ_Struct(AZ_num).DeltaFF0_MeanTrace);
axis([1 xlimmy -0.2 1.2])
%     xline(MarkerSetInfo.Markers.MarkerStart*10,'--', 'LineWidth', 1, 'Color', [0.8500 0.3250 0.0980])
    hold on
%     txt = ['Pre Pr is ',num2str(QuaSOR_Matched_AZ_Structure.All_Boutons(AZ_num).Pre_Pr)];
%     txt1 = ['Post Pr is ',num2str(QuaSOR_Matched_AZ_Structure.All_Boutons(AZ_num).Post_Pr)];
%     txt2 = ['Unc Norm is ', num2str(QuaSOR_Matched_AZ_Structure.All_Boutons(AZ_num).Unc_Norm)];
%     txt3 = ['Brp Norm is ', num2str(QuaSOR_Matched_AZ_Structure.All_Boutons(AZ_num).Brp_Norm)];
% 
%     text(0,1,txt,'Color','red','FontSize',8)
%     text(0,0.8,txt1, 'Color','red','FontSize',8)
%     text(NumIndices*10 - 400,1,txt2,'Color','black','FontSize',8)
%     text(NumIndices*10 - 400, 0.8, txt3, 'Color','black','FontSize',8)
    
y = QuaSOR_Matched_AZ_Structure.Is.Episode(1).AZ_Struct(AZ_num).DeltaFF0_MeanTrace;
x = [1:length(y)];
    ylabel('DeltaF/F0')
    plot(x,y)
   hold on
    peaks2 = QuaSOR_Airy_Pair_Structure_Refined(AZ_num).Modality_AiryMatch(5).Recording.Refine_Coords_Orig_byEpisodeFrame(:,3);
    peaks_by_train = zeros(length(peaks2),1);
    counter = 0;
    for ii = 1:length(peaks2)
        if ii == 1
        peaks_by_train(ii,1) = peaks2(ii,1)
        else
            if peaks2(ii)>peaks2(ii-1)
                frameadd = counter*520;
                peaks_by_train(ii) = peaks2(ii)+frameadd;
            else
                counter = counter +1;
                frameadd = counter*520;
                peaks_by_train(ii) = peaks2(ii)+frameadd;
            end
        end
    end
%     peaks2 = zeros(length(peaks2),2);
%     peaks2(:,1) = peaks2;
%     peaks2(:,1) = peaks2(:,1)*5;
    ytemp = y.';
    peakstrain1_ind = peaks_by_train(:,1)<=520;
    peakstrain1 = peakstrain1_ind.*peaks_by_train;
    peakstrain1 = nonzeros(peakstrain1);
%     
%     peaksy = ytemp(peakstrain1,1);
%     scatter(peakstrain1,peaksy);
%      box on
%     
%     subplot(3,4,9)
%     title('Traces Pre')w
%     for episode_num = 1:NumIndices
%         x = 1:10; 
%         y = QuaSOR_Matched_AZ_Structure.Ib.Episode(episode_num).AZ_Struct(AZ_num).DeltaFF0_MeanTrace;
%        % y = QuaSOR_Matched_AZ_Structure.Is.Episode(episode_num).AZ_Struct(AZ_num).DeltaFF0_MeanTrace;
%         hold on 
%         plot(x,y)
%     end
%     ylim([0 1.2])
%     box on
%     
%     subplot(3,4,10)
%     title('Traces Post')
%     for episode_num = MarkerSetInfo.Markers.MarkerStart: NumIndices
%         x = 1:10; 
%         y = QuaSOR_Matched_AZ_Structure.Ib.Episode(episode_num).AZ_Struct(AZ_num).DeltaFF0_MeanTrace;
%        % y = QuaSOR_Matched_AZ_Structure.Is.Episode(episode_num).AZ_Struct(AZ_num).DeltaFF0_MeanTrace;
%         hold on 
%         plot(x,y)
%     end
%     ylim([0 1.2])
%     box on
%     
%     subplot(3,4,11)
%     x = [1:10];
%     plot(x, mean_pre_trace)
%     hold on
%     plot(x, mean_post_trace)
%     title('Average Trace Pre and Post')
%     %legend('Pre','Post')
%     hold off
%     box on 
%     
    set(gcf, 'PaperUnits', 'inches');
    x_width=10 ;y_width=6
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
    
    baseFileName = sprintf('AZnumber_%d.png',AZ_num);
    fullFileName = fullfile('AZs',baseFileName);
    saveas(gcf,fullFileName); 
    
    baseFileName = sprintf('AZnumber_%d.emf',AZ_num);
    fullFileName = fullfile('AZs',baseFileName);
    saveas(gcf,fullFileName); 
    
    baseFileName = sprintf('AZnumber_%d.eps',AZ_num);
    fullFileName = fullfile('AZs',baseFileName);
    saveas(gcf,fullFileName); 

end
