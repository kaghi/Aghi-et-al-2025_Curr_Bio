imagefiles = dir('*.tif');      
nfiles = length(imagefiles); 

for ii=1:nfiles
   currentfilename = imagefiles(ii).name;
   currentimage = imread(currentfilename);
   images{ii} = currentimage;
end

movieReg = imregister(images{1,2}, images{1,1}, 'affine');
for kk = 1:ii
    matrix = images{1,kk};
    meanval(kk) = mean2(matrix);
end

plot(meanval(1:100));

fftmat = fft(meanval);
emat = 1:length(fftmat);
plot(emat, fftmat);