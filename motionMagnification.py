from linearMotionMagnification import *
from phaseBasedMotionMagnification import *
import time

# Select the video file
fileName = 'Video_08.mp4'

# Select the arguments
# Magnification factor
alpha = 30
# Value of lambda_c
lambda_c = 16
# Sampling rate in fps
samplingRate = 30
chromeAttenuation = 0.1
# Low cutoff frequency in Hz
lowFreq = 0.4
# High cutoff frequency in Hz
highFreq = 1
#  Select the filter for LVMM (butter for motion and ideal for color magnification)
filter = "butter"

# LVMM
print("working on Eulerian based magnification")
# We use timer function to calculate computation time
start = time.time()
if filter == "butter":
    # Call the butterworth filter based magnification function
    try:
        videoMagnificationButterWorthFilter(fileName, alpha, lambda_c, samplingRate, chromeAttenuation, lowFreq, highFreq)
    except Exception as e:
        print("Failed in processing", fileName, "due to", e)

elif filter == "ideal":
    # Call the ideal filter based magnification function
    try:
        videoMagnificationIdealFilter(fileName, alpha, samplingRate, chromeAttenuation, lowFreq, highFreq)
    except Exception as e:
        print("Failed in processing", fileName, "due to", e)

end = time.time()
x = end - start
print('Total time: ', x)

# PBVMM
print("working on Phase Based magnification")
start = time.time()
# Call the phase based magnification function
try:
    phaseBasedMagnification(fileName, alpha, samplingRate, lowFreq, highFreq)
except Exception as e:
    print("Failed in processing", fileName, "due to", e)

end = time.time()
x = end - start

print('Total time: ', x)



