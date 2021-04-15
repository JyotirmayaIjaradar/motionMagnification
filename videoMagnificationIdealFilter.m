% videoMagnificationIdealFilter(vidFile, outDir, alpha, lambda_c,
%                                     wl, wh, samplingRate, chromAttenuation)
% Ideal bandpass filter based on Laplacian pyramid

function videoMagnificationIdealFilter(vidFile, outDir ...
    ,alpha, lambda_c, fl, fh ...
    ,samplingRate, chromAttenuation)
    

    [~,vidName] = fileparts(vidFile);

    outName = fullfile(outDir,[vidName '-ideal-from-' num2str(fl) ...
                       '-to-' num2str(fh) '-alpha-' num2str(alpha) ...
                       '-lambda_c-' num2str(lambda_c) '-chromAtn-' ...
                       num2str(chromAttenuation) '.avi']);

    % Read the input video
    vid = VideoReader(vidFile);
    % Extract and save the video information (video height, width, frame number etc.)
    vidHeight = vid.Height;
    vidWidth = vid.Width;
    nChannels = 3;
    fr = vid.FrameRate;
    len = vid.NumberOfFrames;
    temp = struct('cdata', zeros(vidHeight, vidWidth, nChannels, 'uint8'), 'colormap', []);

    startIndex = 1;
    endIndex = len-10;

    vidOut = VideoWriter(outName);
    vidOut.FrameRate = fr;

    open(vidOut)


    % calculate the Laplacian pyramid for each frame of the video
    [pyr_stack, pind] = build_Lpyr_stack(vidFile, startIndex, endIndex);

    filtered_stack = ideal_bandpassing(pyr_stack, 3, fl, fh, samplingRate);


    % amplify each spatial frequency bands according to the theory of my report
    ind = size(pyr_stack(:,1,1),1);
    nLevels = size(pind,1);
    
    delta = lambda_c/8/(1+alpha);
    
    % the factor to boost alpha above the bound we have in the report
    exaggeration_factor = 2;
    
    % compute the representative wavelength lambda for the lowest spatial 
    % freqency band of Laplacian pyramid
    % here 3 is experimental constant,  you can adjust it based on your result 
    lambda = (vidHeight^2 + vidWidth^2).^0.5/3; 

    for l = nLevels:-1:1
      indices = ind-prod(pind(l,:))+1:ind;
      % Calculate the updated alpha for this level
      currAlpha = lambda/delta/8 - 1;
      currAlpha = currAlpha*exaggeration_factor;
       
      % ignore the highest and lowest frequency band
      if (l == nLevels || l == 1) 
          filtered_stack(indices,:,:) = 0;
      elseif (currAlpha > alpha)  % representative lambda exceeds lambda_c
          filtered_stack(indices,:,:) = alpha*filtered_stack(indices,:,:);
      else
          filtered_stack(indices,:,:) = currAlpha*filtered_stack(indices,:,:);
      end
      
      ind = ind - prod(pind(l,:));
      % now go one level down on pyramid, 
      % representative lambda will reduce by factor of 2
      lambda = lambda/2; 
    end
    
    %% Render the video

    % reconstruct the output video
    k = 0;
    for i=startIndex+1:endIndex
        i
        k = k+1;
        temp.cdata = read(vid, i);
        [rgbframe,~] = frame2im(temp);
        rgbframe = im2double(rgbframe);
        frame = rgb2ntsc(rgbframe);

        filtered = zeros(vidHeight,vidWidth,3);

        filtered(:,:,1) = reconLpyr(filtered_stack(:,1,k),pind);
        filtered(:,:,2) = reconLpyr(filtered_stack(:,2,k),pind)*chromAttenuation;
        filtered(:,:,3) = reconLpyr(filtered_stack(:,3,k),pind)*chromAttenuation;

        filtered = filtered+frame;

        frame = ntsc2rgb(filtered);

        frame(frame > 1) = 1;
        frame(frame < 0) = 0;


        writeVideo(vidOut,im2uint8(frame));
    end


    close(vidOut);

end
