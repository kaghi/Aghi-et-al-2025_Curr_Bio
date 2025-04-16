imadjust(tifanimal(:,:,1));

[optimizer, metric] = imregconfig('multimodal');

moving_reg = imregister(tifanimal(:,:,2),ans,'affine',optimizer,metric );

for i = 1:384
    moving_reg(:,:,i) = imregister(tifanimal(:,:,i), ans, 'affine', optimizer, metric);
end