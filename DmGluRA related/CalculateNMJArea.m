%Requires _Quasor_Data.mat and the BoutonArray structure
% %Ib
% blb_l = size(BoutonArray(1).BorderLine_Orig);
% areab = 0;
% for ii = 1:bls_l(2)
% xbvalues = BoutonArray(1).BorderLine_Orig{1, 2}.BorderLine(:,1);
% ybvalues = BoutonArray(1).BorderLine_Orig{1, 2}.BorderLine(:,2);; 
% scatter(xbvalues,ybvalues)
% areab = polyarea(xbvalues, ybvalues);
% end
% area2b = areab*(0.253968);%The value we are using here is the measurement of each pixel, depending on the camera detector. Scale in um
% 
% %Is
% bls_l = size(BoutonArray(2).BorderLine);
% areas = 0;
% for ii = 1:bls_l(2)
% xsvalues = BoutonArray(2).BorderLine_Orig{1,ii}.BorderLine(:,1);
% ysvalues = BoutonArray(2).BorderLine_Orig{1,ii}.BorderLine(:,2); 
% areas1 = polyarea(xsvalues, ysvalues);
% areas = areas1 + areas;
% end
% area2s = areas*(0.253968);%The value we are using here is the measurement of each pixel, depending on the camera detector. Scale in um
% 
% scatter(xbvalues,ybvalues)
% hold on
% scatter(xsvalues,ysvalues);

imagesc(BoutonMerge.BoutonArray(1).AllBoutonsRegion)
hold on 
imagesc(BoutonMerge.BoutonArray(2).AllBoutonsRegion)

area1b = sum(BoutonMerge.BoutonArray(1).AllBoutonsRegion(:) == 1);
area1b_um = area1b*(0.253968);
area1s = sum(BoutonMerge.BoutonArray(2).AllBoutonsRegion(:) == 1);
area1s_um = area1s*(0.253968);
