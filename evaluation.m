clear;

% this script is written for evaluation mertics
% read the original and result video files
obj = VideoReader('D:/Download/Source and Result Videos/Video_02.mp4');
frame_num = obj.NumFrames;
obj2 = VideoReader('D:/Download/Source and Result Videos/testResultEV/Video_02-FIRWindowBP-band0.40-1.00-sr30-alpha15-mp0-sigma20-scale0.50-frames1-127-octave.avi');
frame_num_2 = obj2.NumFrames;

% Loop to extract frames
for i = 1:frame_num
 frames = read(obj,i);
 imwrite(frames,['D:/Download/Source and Result Videos/data/Image' int2str(i), '.jpg']);
end
for j = 1:frame_num_2
 frames2 = read(obj2,j);
 imwrite(frames2,['D:/Download/Source and Result Videos/data/Image_out' int2str(j), '.jpg']);
end
 
 % set the initial evaluation value to 0 
value = 0;
pvalue = 0;
svalue = 0;
 
for i = 1:frame_num_2 % normally frames are same for two videos
original = imread(['D:/Download/Source and Result Videos/data/Image' int2str(i), '.jpg']);
magnified = imread(['D:/Download/Source and Result Videos/data/Image_out' int2str(i), '.jpg']);
% resize the frames to same size
o = imresize(original,[320,320]);
m = imresize(magnified,[320,320]);
% get the ssim, psnr value using the matlab function
[ssimval,ssimmap] = ssim(m,o);
[peaksnr, snr] = psnr(m, o);
imshow(ssimmap,[])
title(['SSIM Map witth SSIM Value: ',num2str(ssimval)])
% add values for all the frames
value = value + ssimval;
pvalue = pvalue + peaksnr;
svalue = svalue + snr;

end
% calculate the final evaluation value
ssimvalue = value/frame_num;
psnrvalue = pvalue/frame_num;
snrvalue = svalue/frame_num;

