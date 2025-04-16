NormalizeImage = ones(332,1219)*1000;
uintNorm = uint16(NormalizeImage);
NormalizedImage = sampleFrame - uintNorm;
save('demon.mat', 'demon')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save(file_name);
tracktiff = Tiff(file_name, 'w');

file_name = 'ok6-gc6s-animal2-2_16bit_registerNMJ_4'

demontiff = 'UAS-GC6S_Control_161012_P01_VNC_M01_16bit.tif';
tracktiff = imwrite(track(:,:,1), strcat(file_name),'tif')
for n=2:2000
    imwrite(track(:,:,n),strcat(file_name),'tif','WriteMode' ,'append');
end 


movefile('tracktiff.tiff',file_name);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
obj = Tiff(file_name,'a');
imwrite(demon, 'demonstack.tif');
imwrite(demon(:,:,1), 'demonstack.tiff');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = imread('demonstack.tiff');
[L,N] = superpixels(A(:,:,1),500);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Normalize traces
frame_number_T = 1:2000;
frame_number = frame_number_T.';
for n=1:2000
    fluorvalues = cells_mean(:,:);
end
for n=1:2000
    fittedcurve = fit(frame_number,cells_mean(:,n),'poly2'); 
end

%%%%%%%%%%%%%%%%
t = imwrite(demon,strcat(file_name),'tif');
t.write(demon);

%%%%%%%%%
demond = double(demon);
for i=1:2000
  tiff = demond(:, :, i);
  file_name = sprintf('UAS-GC6S_Control_161012_P01_VNC_M01_16bit.tif', i);
  imwrite(tiff,file_name,'WriteMode', 'append')
end

%%%%%%%
for i=1:2000
  tiff = demon(:, :, i);
  file_name = sprintf('UAS-GC6S_Control_161012_P01_VNC_M01_16bit.tif', i);
  imwrite(tiff,file_name,'WriteMode', 'append')
end
%%%%%%%
%Convert raw traces into deltaF/F
for ii =1:2000