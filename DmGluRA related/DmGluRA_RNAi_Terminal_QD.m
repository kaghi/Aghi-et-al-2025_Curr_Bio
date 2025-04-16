%%Event detection data PixelMax_Struct.All_Location_Coords_byOverallStim  
%%Quasor detection data QuaSOR_Data.All_Location_Coords_byOverallStim
%Also need BorderLine from the _Analysis_Setup.mat
clear all
%If 0.2 Hz data
All_Events_By_Stim = PixelMax_Struct.All_Location_Coords_byOverallStim;
%If 5 Hz data
All_Events_By_Stim = PixelMax_Output.All_Location_Coords_byStim;  

%Once files loaded start here.
%filter right away on which events you want to analyze by stim number
EventsOfInterest_Ind = find(All_Events_By_Stim(:,3)>9 & All_Events_By_Stim(:,3)<101);
%%%

Events_of_Interest_x = All_Events_By_Stim(EventsOfInterest_Ind,1);
Events_of_Interest_y = All_Events_By_Stim(EventsOfInterest_Ind,2);
Events_of_Interest = [Events_of_Interest_x Events_of_Interest_y]

hFig = figure;
hPlot = scatter(Events_of_Interest(:,1),Events_of_Interest(:,2),10,'Filled')
axis equal

    %Now to select events in bouton. Hold shift to select multiple areas.
    clear('brush')
     % create and enable the brush object
        hBrush = brush(hFig);
        hBrush.ActionPostCallback = @OnBrushActionPostCallback;
        hBrush.Enable = 'on';

%Run this after selecting events    
xd = get(hPlot, 'XData');
yd = get(hPlot, 'YData');
brush_data = get(hPlot, 'BrushData');
brush_data = brush_data.';
brush_data = im2double(brush_data);
selected_ind_event = find(brush_data>0);

Terminal_events = zeros(length(selected_ind_event),2);
for ll = 1:2
Terminal_events(:,ll) = Events_of_Interest(selected_ind_event,ll);
end

scatter(Terminal_events(:,1),Terminal_events(:,2),10,'Filled')
axis equal

close all
%%%%%%%%%%%%%%%%%%%%%%%%%Find area of terminal bouton
%Press enter in command window once youve selected bouton so it can loop
%through
clear('brush')
no_borders = length(BorderLine);
areas = [];
count = 1;
%Run this section multiple times for multiple terminal boutons across the
%border lines. ONLY run the above portion at the start of a NEW NMJ. 
for kk = 1:no_borders
hFig = figure;
hPlot = scatter(BorderLine{1,kk}.BorderLine(:,1),BorderLine{1,kk}.BorderLine(:,2))
axis equal
 % create and enable the brush object
    hBrush = brush(hFig);
    hBrush.ActionPostCallback = @OnBrushActionPostCallback;
    hBrush.Enable = 'on';
    pause
   
xd = get(hPlot, 'XData');
yd = get(hPlot, 'YData');
brush_data = get(hPlot, 'BrushData');
brush_data = brush_data.';
brush_data = im2double(brush_data);

selected_ind = find(brush_data>0);

pg_points = zeros(length(selected_ind),2);
for nn = 1:2
    pg_points(:,nn) = (BorderLine{1,kk}.BorderLine(selected_ind,nn))
end

pgon = polyshape(pg_points);
term_bout_area = area(pgon);
%Change the pix_2_area value based on the pixel to um conversion from
%original microscope.
pix_2_area = 0.2112;
%%%%%%
true_term_bout_area = term_bout_area*pix_2_area;
if true_term_bout_area  == 0;
areas{kk,count} = [];
else
areas{kk,count}.polygon = pg_points;
areas{kk,count}.area = true_term_bout_area;

end
count = count+1;
end

%Only run once done
areas = areas(~cellfun('isempty',areas));

%%%%Calculate QD
total_area = 0;
for nnn = 1:length(areas);
selected_area = areas{nnn}.area;
total_area = total_area + selected_area;
end

%Only want to get events from last 90 stims (last 90);
QuantalDensity = (length(Terminal_events))/total_area;
QuantalDensity

%Plot selected boutons
for nn  =1: length(areas)
 x_val = areas{1, nn}.polygon(:,1);
 y_val = areas{1, nn}.polygon(:,2);
  pgon = polyshape(x_val,y_val); 
  plot(pgon)
 axis equal
  hold on
 
end
hold on
hPlot = scatter(BorderLine{1,kk}.BorderLine(:,1),BorderLine{1,kk}.BorderLine(:,2),5)
axis equal
% 
% 
% scatter(BorderLine_Upscale{1,1}.BorderLine_Crop(:,1),BorderLine_Upscale{1, 1}.BorderLine_Crop(:,2))
% 
%     
% xd = get(hPlot, 'XData');
% yd = get(hPlot, 'YData');
% brush_data = get(hPlot, 'BrushData');
% 
% %find area of selected region
% pgon = polyshape(BorderLine_Upscale{1,2}.BorderLine_Crop(:,1) , BorderLine_Upscale{1, 2}.BorderLine_Crop(:,2) ) 
% 
% axis equal
% hold on 
% scatter(All_Events_By_Stim(:,1),All_Events_By_Stim(:,2),10,'Filled')
% axis equal


