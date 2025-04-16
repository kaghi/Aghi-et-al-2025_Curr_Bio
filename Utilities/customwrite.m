for x = 1 : 207
    imwrite(moving_reg(:,:,x),strcat(file_name,'.tif'),'WriteMode','append');
end