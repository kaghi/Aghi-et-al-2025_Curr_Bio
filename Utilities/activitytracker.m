
h_fig = figure;
set(h_fig,'KeyPressFcn',@myfun);

function KeyRecord(src,event)
   record = event.Key;
   disp(event.Key);
end

while ~isDone(VideoReader)
  videoFrame = step(VideoReader);
  step(videoPlayer, videoFrame);
  set(VideoPlayer,'KeyPressFcn',KeyRecord);
end

function myfun(src,event)
   record = event.Key;
   disp(event.Key);
end

function cb(src,evnt)
     persistent ticStartTime;       % creates a persistent variable
     out = sprintf('Key: %s\n',evnt.Key);
     disp(out)
     if isempty(ticStartTime)
         % use last the tic from the key_pressFcn 
         time=toc;
     else
         % use tic that was called within this callback
         time = toc(ticStartTime);
     end
     disp(time)  
     % restart tic for this function
     ticStartTime = tic;
end 
