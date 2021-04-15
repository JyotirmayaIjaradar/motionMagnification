
import cv2
import os
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal
import warnings

# local imports
from common import *

# skip warnings; warning are expected to happen for some division by zero which is handled in the code
warnings.filterwarnings("ignore")

# Compute the complex steerable pyramid. This function is created based the work at Wadhia et al. and their upgraded version of Riesz pyramid.
def complexSteerablePyramid(frame, levels):
    # First the basic construction is based on the LVMM laplacian pyramid
    laplacianPyramid, sizes = buildLaplacian(frame, levels=levels)
    kernelX = np.array([[-0.12,0,0.12],[-0.34, 0, 0.34],[-0.12,0,0.12]], dtype="float32")
    kernelY = np.transpose(kernelX)    
    rieszX = []
    rieszY = []
    sizes = []
    # Based on the frame and level we construct the pyramid
    for level in range(levels-1):
        rieszX.append(cv2.filter2D(laplacianPyramid[level], -1, kernelX))
        rieszY.append(cv2.filter2D(laplacianPyramid[level], -1, kernelY))
        sizes.append(rieszX[level].shape)
    return laplacianPyramid, rieszX, rieszY, sizes          

# As stated in the report, the PBVMM works on the phase differences and amplification applied on phases. This function is to calculate the phase differences and amplitude in the phase
def computePhaseDifferenceAndAmplitude(currentReal, currentX, currentY, previousReal, previousX, previousY):
    
    qConjProdReal = np.multiply(currentReal , previousReal) + np.multiply(currentX , previousX) + np.multiply(currentY , previousY)
    
    qConjProdX = np.multiply(-1*currentReal , previousX) + np.multiply(previousReal , currentX)
    qConjProdY = np.multiply(-1*currentReal , previousY) + np.multiply(previousReal , currentY)
    
    qConjProdAmplitudeTemp = np.power(qConjProdReal, 2) + np.power(qConjProdX, 2) + np.power(qConjProdY, 2)
    qConjProdAmplitude = np.sqrt(qConjProdAmplitudeTemp)    
    
    temp1 = np.divide(qConjProdReal , qConjProdAmplitude)
    temp1 = np.nan_to_num(temp1)
    phaseDifference = np.arccos(temp1)

    temp2Sum = np.power(qConjProdX, 2) + np.power(qConjProdY,2)
    temp2 = np.sqrt(temp2Sum)
    cosOrientation = np.divide(qConjProdX , temp2)
    sinOrientation = np.divide(qConjProdY , temp2)
    cosOrientation =  np.nan_to_num(cosOrientation)
    sinOrientation = np.nan_to_num(sinOrientation)

    phaseDifferenceCos = np.multiply(phaseDifference , cosOrientation)
    phaseDifferenceSin = np.multiply(phaseDifference , sinOrientation)
    
    amplitude = np.sqrt(qConjProdAmplitude)
    
    return phaseDifferenceCos, phaseDifferenceSin, amplitude

# compute the ideal temporal filter. here our main argument is the phase
def IIRTemporalFilter(B,A,phase,register0, register1):
    temporallyFilteredPhase = (B[0] * phase) + register0
    newRegister0 = (B[1] * phase) + register1 - (A[1]*temporallyFilteredPhase)
    newRegister1 = (B[2] * phase)             - (A[2]*temporallyFilteredPhase)
    return temporallyFilteredPhase, newRegister0, newRegister1

# The optiontional amplitude wighting blur finction. Mainly we add the weighted blur to the temporal filtered images
def amplitudeWeightedBlur(temporallyFilteredPhase, amplitude, blurKernel):
    denom = cv2.filter2D(amplitude, -1, blurKernel)
    numer = cv2.filter2D(np.multiply(temporallyFilteredPhase,amplitude), -1, blurKernel)
    spatiallySmoothedTemporallyFilteredPhase = np.divide(numer, denom)
    spatiallySmoothedTemporallyFilteredPhase = np.nan_to_num(spatiallySmoothedTemporallyFilteredPhase)
    return spatiallySmoothedTemporallyFilteredPhase

#  Function to calculate phase shift coefficient
def phaseShiftCoefficientRealPart(rieszReal, rieszX, rieszY, phaseCos, phaseSin):
    temp1 = np.power(phaseCos, 2) + np.power(phaseSin,2)
    phaseMagnitude = np.sqrt(temp1)
    
    expPhaseReal = np.cos(phaseMagnitude)
    
    temp2 = np.divide(phaseCos, phaseMagnitude)
    temp3 = np.divide(phaseSin, phaseMagnitude)
    temp2 = np.nan_to_num(temp2)
    temp3 = np.nan_to_num(temp3)
    
    expPhaseX = np.multiply(temp2 , np.sin(phaseMagnitude))
    expPhaseY = np.multiply(temp3 , np.sin(phaseMagnitude))
    
    result = np.multiply(expPhaseReal, rieszReal) - np.multiply(expPhaseX, rieszX) - np.multiply(expPhaseY, rieszY)
    
    return result

# Function to apply magnification to phases. The argument for this are amplification factor, cutoff frequency etc.
def phaseBasedMagnification(fileToProcess, alpha, samplingRate, lowFreq, highFreq):

    # Select the source the and result directories
    datadir = os.getcwd()
    # Source video directory.
    samplesFn = datadir + "/Source"
    # Result video directory
    resultsFn = datadir + "/results/phaseBased"
    # Split video file to frames
    vidFn = fileToProcess.split(".")[0]
    resParentFn = resultsFn + "/" + vidFn
    # Original video frame temporary directory
    vidFramesFn = resParentFn + "/videoFramesOriginal"
    # Magnified video frame temporary directory
    reconstructedFn = resParentFn + "/magnifiedFrames"

    # create results folder according to processed file
    if os.path.exists(resParentFn) == False:
        os.makedirs(resParentFn)
    if os.path.exists(vidFramesFn) == False:
        os.makedirs(vidFramesFn)
    if os.path.exists(reconstructedFn) == False:
        os.makedirs(reconstructedFn)

    # load video file into frames np array
    # save every frames to local storage as well
    # calculate the video height, width and framenumber
    videoInputFn    = cv2.VideoCapture(samplesFn + "/" + fileToProcess)
    frameCount      = int(videoInputFn.get(cv2.CAP_PROP_FRAME_COUNT))
    frameHeight     = int(videoInputFn.get(cv2.CAP_PROP_FRAME_HEIGHT))
    frameWidth      = int(videoInputFn.get(cv2.CAP_PROP_FRAME_WIDTH))
    frameChannels   = 3
    # Calling the cv2 video codec function for mp4
    four_cc = cv2.VideoWriter_fourcc('m', 'p', '4', 'v')
    # Write the output video
    writer = cv2.VideoWriter(resultsFn + "/" + vidFn + "_phase_magnified_alpha"+ str(alpha) +".mp4", four_cc, 30, (frameWidth, frameHeight), 1)

    # get the video frames based on the video data
    frames = np.zeros((frameCount, frameHeight, frameWidth, frameChannels), dtype="float32")

    success, image = videoInputFn.read()
    # set initial frame index
    frameIndex = 1
    print("Processing the video", fileToProcess , "...")
    print("Extracting original video frames...")
    # Loop to extract all the original video frames and save that to temporary folder
    while success:
        frames[frameIndex-1,:,:,:] = np.float32(image/255.0)
        cv2.imwrite(vidFramesFn + "/frame%04d.png" %(frameIndex), image)
        success, image = videoInputFn.read()
        frameIndex += 1
    print("Extracted %d frames..." %(frameIndex-1))

    # Get butterworth bandpass filter parameters
    # Calculate the Nyquist frequency based on the theory in wadhia et al.
    nyQuistFrequency = samplingRate / 2 
    temporalFilterOrder = 1 
    wl = lowFreq / nyQuistFrequency
    wh = highFreq / nyQuistFrequency
    B, A = signal.butter(temporalFilterOrder, [wl, wh], btype='bandpass')

    # Prepare gaussianKernel used later for smoothing 
    gaussianKernel = cv2.getGaussianKernel(7, 2)
    gaussianKernel2D = gaussianKernel * np.transpose(gaussianKernel)

    # Assumed number of levels (here we assumed to 7 you can assume different level based on your video)
    pyrLevels = 7

    # First frame, convert to YIQ space & use only first channel
    previousFrame = convertBGR2YIQ(frames[0])[:,:,0]

    # sizeGuidance is just needed to ease setup of upcoming arrays
    _, sizeGuidance = buildLaplacian(previousFrame, levels=pyrLevels)

    # Pyramid construction for the first frame
    previousLaplacian, previousRieszX, previousRieszY, sizeGuidance = complexSteerablePyramid(previousFrame, levels=pyrLevels)

    # we use the term Riesz to keep it different and as stated the implementation is influenced by Riesz upgradation of Phase based based motion magnification of MIT
    rieszPyrLevels = pyrLevels - 1

    # Setup variables needed for the main for loop for magnification
    phaseCos = [0] * rieszPyrLevels
    phaseSin = [0] * rieszPyrLevels
    register0Cos = [0] * rieszPyrLevels
    register1Cos = [0] * rieszPyrLevels
    register0Sin = [0] * rieszPyrLevels
    register1Sin = [0] * rieszPyrLevels
    motionMagnifiedLaplacianPyr = [0] * (rieszPyrLevels+1)
    # loop to get every pyramid level coefficients
    for k in range(rieszPyrLevels):
        phaseCos[k] = np.zeros(sizeGuidance[k])
        phaseSin[k] = np.zeros(sizeGuidance[k])
        register0Cos[k] = np.zeros(sizeGuidance[k])
        register1Cos[k] = np.zeros(sizeGuidance[k])
        register0Sin[k] = np.zeros(sizeGuidance[k])
        register1Sin[k] = np.zeros(sizeGuidance[k])
        motionMagnifiedLaplacianPyr[k] = np.zeros(sizeGuidance[k])

    # Save first frame as original frame, no magnification could be done here
    saveVideo(frames[0],reconstructedFn, writer, 0)

    print("Started phase based magnification...")
    # progress is used for cosmetic display of how much of the video is done
    progress = np.linspace(0,frameCount-1, 11, dtype="int32")

    # main for loop to magnifiy frames starting from the second frame
    for frameIdx in range(1,frameCount):
        currentFrame = convertBGR2YIQ(frames[frameIdx])[:,:,0]
        
        currentLaplacian, currRieszX, currRieszY, _ = complexSteerablePyramid(currentFrame, levels=pyrLevels)
        for k in range(rieszPyrLevels):
            phaseDifferenceCos, phaseDifferenceSin, amplitude = computePhaseDifferenceAndAmplitude(currentLaplacian[k], currRieszX[k], currRieszY[k], previousLaplacian[k], previousRieszX[k], previousRieszY[k])
            
            phaseCos[k] = phaseCos[k] + phaseDifferenceCos
            phaseSin[k] = phaseSin[k] + phaseDifferenceSin
            
            phaseFilteredCos, register0Cos[k], register1Cos[k] = IIRTemporalFilter(B,A,phaseCos[k], register0Cos[k], register1Cos[k])
            phaseFilteredSin, register0Sin[k], register1Sin[k] = IIRTemporalFilter(B,A,phaseSin[k], register0Sin[k], register1Sin[k])
            
            phaseFilteredCos = amplitudeWeightedBlur(phaseFilteredCos, amplitude, gaussianKernel2D)
            phaseFilteredSin = amplitudeWeightedBlur(phaseFilteredSin, amplitude, gaussianKernel2D)
            
            phaseMagnifiedFilteredCos = alpha * phaseFilteredCos
            phaseMagnifiedFilteredSin = alpha * phaseFilteredSin
            
            motionMagnifiedLaplacianPyr[k] = phaseShiftCoefficientRealPart(currentLaplacian[k], currRieszX[k], currRieszY[k], phaseMagnifiedFilteredCos, phaseMagnifiedFilteredSin)
        
        motionMagnifiedLaplacianPyr[-1] = currentLaplacian[-1]
        
        previousLaplacian, previousRieszX, previousRieszY = np.copy(currentLaplacian), np.copy(currRieszX), np.copy(currRieszY)

        # reconstruct the magnified video
        reconstructedChannel = reconstructVideoFrames(motionMagnifiedLaplacianPyr)
        reconstructedFrame = np.copy(convertBGR2YIQ(frames[frameIdx]))
        reconstructedFrame[:,:,0] = reconstructedChannel
        reconstructedFrameRGB = convertYIQ2RGB(reconstructedFrame)
        saveVideo(reconstructedFrameRGB, reconstructedFn, writer, frameIdx)
        if frameIdx in progress:
            print("Completed the PBVMM " , np.where(progress==frameIdx)[0][0] * 10 ,"%" , "of the video")

    # release the writer used for video after finishing last frame
    writer.release()
    print("Done...")