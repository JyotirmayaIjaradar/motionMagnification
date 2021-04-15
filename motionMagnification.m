clear;

% select the source file directory
dataDir = 'D:/sourceVideos';
% set the output file directory 
resultsDir = 'D:/resultVideos/';
mkdir(resultsDir);
scaleAndClipLargeVideos = true;

%% Motion Magnification based on Linear Eulerian motion magnification (LVMM)
% Here I have given the optimised parameter for all our example to reproduce
% the result. For using other videos you need to adjust the amplification
% factor , cutoff frequency etc.
% amplify_spatial_lpyr_temporal_butter(inFile, resultsDir, amplification factor, lammda_c, low cutoff frequency, high cutoff frequency, frame rate, chrome attenuation(we kept deafult 0.1 for all));
% amplify_spatial_Gdown_temporal_ideal(inFile,resultsDir, amplification factor, levels, low cutoff frequency, high cutoff frequency, frame rate, chrome attenuation);
% Based on the above example you can set your parameter. 
% the temporal butter filter is for motionn magnification.
% the temporal ideal filter is for color magnification
% tic-toc feature is used for computation time calculation

%% Test Variance 1: (motion)
tic
inFile = fullfile(dataDir,'Video_08.mp4');
fprintf('Start working for video file  %s\n', inFile);
% Select the magnification factor
alpha = 20;
% Select the lambda value, you can get idea from  the theory
lambda_c = 20;
% Select the low cutoff frquency in Hz
wl = 0.6;
% Select the high cutoff frquency in Hz
wh = 1;
% Select the frame rate in fps
frameRate = 30;
amplify_spatial_lpyr_temporal_butter(inFile, resultsDir, alpha, lambda_c, wl, wh, frameRate, 0.1);
toc

%% Test Variance 2: (color)
tic
inFile = fullfile(dataDir,'Video_08.mp4');
fprintf('Start working for video file  %s\n', inFile);
% Select the magnification factor
alpha = 80;
% Select the low cutoff frquency in Hz
wl = 0.6;
% Select the high cutoff frquency in Hz
wh = 0.8;
% Select the frame rate in fps
frameRate = 30;
% Select the levl for ideal filter
level = 6;
amplify_spatial_Gdown_temporal_ideal(inFile,resultsDir, alpha, level, wl, wh, frameRate, 1);
toc


%% Test Variance 3:(motion)
tic
inFile = fullfile(dataDir,'Sample-01.avi');
fprintf('Start working for video file  %s\n', inFile);
% Select the magnification factor
alpha = 20;
% Select the lambda value, you can get idea from  the theory
lambda_c = 20;
% Select the low cutoff frquency in Hz
wl = 0.8;
% Select the high cutoff frquency in Hz
wh = 0.9;
% Select the frame rate in fps
frameRate = 30;
amplify_spatial_lpyr_temporal_butter(inFile, resultsDir, alpha, lambda_c, wl, wh, frameRate, 0.1);
toc

%% Test Variance 4:(color)
tic
inFile = fullfile(dataDir,'Sample-01.avi');
fprintf('Start working for video file  %s\n', inFile);
% Select the magnification factor
alpha = 50;
% Select the low cutoff frquency in Hz
wl = 0.8;
% Select the high cutoff frquency in Hz
wh = 1;
% Select the frame rate in fps
frameRate = 30;
% Select the levl for ideal filter
level = 6;
amplify_spatial_Gdown_temporal_ideal(inFile,resultsDir, alpha, level, wl, wh, frameRate, 1);
toc

%% MOtion Magnification based on Phase Based Motion Magnification (PBVMM)
% Again here I have fix the parameters to reproduce our result
% For new video here again you need to set the parameters.

%% Case1
tic
inFile = fullfile(dataDir, 'Video_08.mp4');
% Select the sampling rate in Hz
samplingRate = 16; 
% Select the low cutoff frquency in Hz
loCutoff = 0.5;    
% Select the low cutoff frquency in Hz
hiCutoff = 1; 
% Select the magnification factor
alpha = 10;   
% Select the sigma value
sigma = 8;  
% select the pyramid type for decomposition (octave, halfOctave,)
pyrType = 'octave'; 
phaseAmplify(inFile, alpha, loCutoff, hiCutoff, samplingRate, resultsDir,'sigma', sigma,'pyrType', pyrType,'scaleVideo', 1);
toc

%% Case2
tic
inFile = fullfile(dataDir, 'Sample-01.avi');
% Select the sampling rate in Hz
samplingRate = 40; 
% Select the low cutoff frquency in Hz
loCutoff = 0.9;    
% Select the low cutoff frquency in Hz
hiCutoff = 1; 
% Select the magnification factor
alpha = 20;   
% Select the sigma value
sigma = 20;  
% select the pyramid type for decomposition (octave, halfOctave,)
pyrType = 'octave'; 
phaseAmplify(inFile, alpha, loCutoff, hiCutoff, samplingRate, resultsDir,'sigma', sigma,'pyrType', pyrType,'scaleVideo', 1);
toc

%% Case3
tic
inFile = fullfile(dataDir, 'Video_05.mp4');
% Select the sampling rate in Hz
samplingRate = 20; 
% Select the low cutoff frquency in Hz
loCutoff = 0.4;    
% Select the low cutoff frquency in Hz
hiCutoff = 0.8; 
% Select the magnification factor
alpha = 20;   
% Select the sigma value
sigma = 10;  
% select the pyramid type for decomposition (octave, halfOctave,)
pyrType = 'octave'; 
phaseAmplify(inFile, alpha, loCutoff, hiCutoff, samplingRate, resultsDir,'sigma', sigma,'pyrType', pyrType,'scaleVideo', 1);
toc

%% Case4
tic
inFile = fullfile(dataDir, 'Sample-03.avi');
% Select the sampling rate in Hz
samplingRate = 75; 
% Select the low cutoff frquency in Hz
loCutoff = 0.6;    
% Select the low cutoff frquency in Hz
hiCutoff = 0.9; 
% Select the magnification factor
alpha = 30;   
% Select the sigma value
sigma = 50;  
% select the pyramid type for decomposition (octave, halfOctave,)
pyrType = 'octave'; 
phaseAmplify(inFile, alpha, loCutoff, hiCutoff, samplingRate, resultsDir,'sigma', sigma,'pyrType', pyrType,'scaleVideo', 1);
toc

