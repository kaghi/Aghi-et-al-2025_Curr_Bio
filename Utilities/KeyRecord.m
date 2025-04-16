
%Load in file
file_name = ''
ImageArray = bfopen(file_name)
for i=1:2000
    movie(:,:,i) = ImageArray{1,1}{i,1};
end 
r2cmovie = im2double(movie);

%Convert your tif file into an AVI so the movie reader can parse out the
%frames
v = VideoWriter(strcat(file_name,'.avi','Uncompressed AVI');
 open(v);
for k=1:length(movie)      % assumes 10 images to write to file
     writeVideo(v,randi(255,100,200,'uint8'));
 end
%Record points in which the nerve terminals are on and off
offframes = ones(1,length(numFrames));
onframes = ones(1,length(numFrames));

while ~isDone(VideoReader)
  videoFrame = step(VideoReader);
  step(VideoPlayer, videoFrame);
  switch event.Key
    case 'leftarrow'
       disp(str("On state"));
       for ii = 1:length(numFrames);
           currFrame = get(videoFrame);
           if currFrame = ii;
               onframes(1,ii) = 1;
           end
       end 
        
    case 'rightarrow'
        disp(str("Off state"));
        for ii = 1:length(numFrames);
           currFrame = get(videoFrame); 
           if currFrame = ii;
               offframes(1,ii) = 2;
           end
        end 
  end
end 


%Calculate percent active time for each presynaptic terminal you have
offframes + onframes = totalrecord;

%specify if it is Ib focused or Is focused. Make sure there are no spaces
%between the prompt and your response (I.e., answer will look like "....as
%IsY"

prompt = 'Did you focus on Ib inputs? Type Y for yes, otherwise it will automatically label as Is';
str = input(prompt,'s');
if str == 'Y';
    filename_new = strcat(filename, '_Ib_rec');
else
    filename_new = strcat(filename, '_Is_rec');
end

prompt = 'What number nerve terminal is this?';
str2 = input(prompt);
newfilename = strcat(filename_new, str2)

function KeyRecord(src,evnt)
   persistent record 
   record = evnt.Key;
   disp(evnt.Key);
end

