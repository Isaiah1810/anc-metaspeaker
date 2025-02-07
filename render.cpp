
#include <array>
#include <Bela.h>
#include <libraries/Fft/Fft.h>
#include <libraries/Scope/Scope.h>
#include <cmath>
#include <random>

// Constants
//const float carrierFreq = 40000.0f;
const float carrierFreq = 0.0f;
const float carrierAmp = 2.0f;
const float stepSize = 0.05f;
const unsigned int primaryFilterSize = 256;
const unsigned int secondaryFilterSize = 256;
const float convergenceThreshold = 0.01f;

// Global Variables
unsigned int currSample = 0;
std::vector<float> refBlock, errorBlock, antiNoiseBlock, trainingNoiseBlock;
std::vector<float> primaryFilter(primaryFilterSize, 0.01f);
std::vector<float> secondaryFilter(secondaryFilterSize, 0.01f);
std::vector<float> prevSecondaryFilter(secondaryFilterSize, 0.0f);
std::vector<float> output;
Scope scope;
Fft fft;
bool doNoiseControl = false;


// Constants for simulation
float gFrequency1 = 440.0;	// Frequency of the sine wave in Hz
float gAmplitude = .30f;		// Amplitude of the sine wave (1.0 is maximum)
float sampleRate;
float gPhase1 = 0;
// Simulated filter
std::vector<float> testFilter = {1.0, 0.0, 0.0, 0.0, 0.5}; 

// Precomputed modulation buffer
std::vector<float> modulationBuffer;

// Function to initialize modulation buffer
void initializeModulationBuffer(size_t blockSize, float sampleRate) {
    modulationBuffer.resize(blockSize);
    for (size_t n = 0; n < blockSize; ++n) {
        modulationBuffer[n] = carrierAmp * sinf(2.0f * M_PI * carrierFreq * n / sampleRate);
    }
}

/**
 * Initializes vectors, fft module, osciliscope and modulation buffer 
 * 
 **/
bool setup(BelaContext *context, void *userData) {
    size_t blockSize = context->audioFrames;
    
    refBlock.resize(blockSize, 0.0f);
    errorBlock.resize(blockSize, 0.0f);
    
    antiNoiseBlock.resize(blockSize, 0.0f);
    trainingNoiseBlock.resize(blockSize, 0.0f);
    output.resize(blockSize, 0.0f);
    initializeModulationBuffer(blockSize, context->audioSampleRate);
    
	sampleRate = context->audioSampleRate;

    fft.setup(blockSize);
    scope.setup(4, context->audioSampleRate);
    return true;
}

/**
 * Given the specific block since start, generate some known random noise
 * 
 * &noise: Reference to a vector to where the noise will be stored into 
 * currSample: an integer representing the number of samples passed since begining of the program
 **/
void generateNoise(std::vector<float> &noise, unsigned int currSample){
    static std::default_random_engine generator;
    static std::uniform_real_distribution<float> distribution(-1.0, 1.0);
    for (size_t n = 0; n < noise.size(); ++n){
        noise[n] = distribution(generator);
    }
    return;
}


/**
 *  Applies FxLMS adaptive filtering algorithm on a specific filter vector
 * 
 *&reference: A vector containing the reference signal to be compared to
 *&error: A vector containing the error mic signal to find the difference from 
 *&filter: A vector that represents filter coefficients to be learned/updated
 *stepSize: Learning rate of FxLMS
 *&output: A vector to where to copy output of filter with reference
 *&copyOut: Whether to write to the output vector
 **/
void applyFxLMS(const std::vector<float> &reference, const std::vector<float> &error,
                std::vector<float> &filter, float stepSize, std::vector<float> &output, bool copyOut) {
    size_t filterOrder = filter.size();
    size_t numSamples = reference.size();

    for (size_t n = 0; n < numSamples; ++n) {
        float filterOutput = 0.0f;
        for (size_t k = 0; k < filterOrder && n >= k; ++k) {
            filterOutput += filter[k] * reference[n - k];
        }

        float errorValue = error[n] - filterOutput;
        for (size_t k = 0; k < filterOrder && n >= k; ++k) {
            filter[k] += stepSize * errorValue * reference[n - k];
        }

        if (copyOut) output[n] = filterOutput;
    }
}

void processAntiNoise(const std::vector<float> &antiNoise, std::vector<float> &output, BelaContext *context) {
    fft.fft(antiNoise);
    size_t blockSize = output.size();
    
	std::vector<float> reBuffer(context->audioFrames, 0.0f);
	std::vector<float> imBuffer(context->audioFrames, 0.0f);

	
    // Remove negative frequencies
    for (size_t n = blockSize / 2; n < blockSize; ++n) {
        reBuffer[n] = fft.fdr(n);
        imBuffer[n] = fft.fdi(n);
    }

    fft.ifft(reBuffer, imBuffer);
    for (size_t n = 0; n < blockSize; ++n) {
        output[n] = fft.td(n) * modulationBuffer[n];
    }
}

bool checkConvergence(const std::vector<float> &filter, const std::vector<float> &prevFilter) {
    float delta = 0.0f;
    for (size_t n = 0; n < filter.size(); ++n) {
        delta += (filter[n] - prevFilter[n]) * (filter[n] - prevFilter[n]);
    }
    return delta < convergenceThreshold;
}


void convolve(const std::vector<float>& signal, 
              const std::vector<float>& filter, 
              std::vector<float>& output) {
    int signalSize = signal.size();
    int filterSize = filter.size();
    
    // Ensure the output vector has the correct size
    output.assign(signalSize, 0.0f);

    // Perform FIR filtering (convolution with fixed output size)
    for (int i = 0; i < signalSize; i++) {
        for (int j = 0; j < filterSize; j++) {
            if (i - j >= 0) {  // Prevent out-of-bounds access
                output[i] += signal[i - j] * filter[j];
            }
        }
    }
}




void render(BelaContext *context, void *userData) {
    size_t blockSize = context->audioFrames;
    // for (size_t n = 0; n < blockSize; ++n) {
    //     refBlock[n] = audioRead(context, n, 0);  // Reference mic
    //     errorBlock[n] = audioRead(context, n, 2); // Error mic
    // }

	for (unsigned int n = 0; n < context->audioFrames; n++) {
	    // TODO: Calculate a sample of the sine wave
	    //       Start by copying the calculation from the 'sine-generator' project
	    gPhase1 += 2.0 * M_PI * gFrequency1 / sampleRate;
	    if (gPhase1 > 2.0 * M_PI) {
	    	gPhase1 -= 2.0 * M_PI;
		}
		refBlock[n] = gAmplitude * sin(gPhase1);
		}
		
	std::vector<float> tempErrorBlock(blockSize, 0.0f);
	for (size_t i = 0; i < blockSize; ++i) {
	    for (size_t j = 0; j < testFilter.size() && (i + j) < blockSize; ++j) {
	        tempErrorBlock[i + j] += refBlock[i] * testFilter[j];
	    }
	}
errorBlock = tempErrorBlock;
		
    if (!doNoiseControl) {
        prevSecondaryFilter = secondaryFilter;
        generateNoise(trainingNoiseBlock, currSample);
        applyFxLMS(trainingNoiseBlock, errorBlock, secondaryFilter, stepSize, antiNoiseBlock, true);
        doNoiseControl = checkConvergence(secondaryFilter, prevSecondaryFilter);
        rt_printf("Secondary Path Learned \n");
        
    } else {
    	//convolve(refBlock, secondaryFilter, output);
        applyFxLMS(refBlock, errorBlock, primaryFilter, stepSize, antiNoiseBlock, true);
    }

    processAntiNoise(antiNoiseBlock, output, context);
    for (size_t n = 0; n < blockSize; ++n) {
        audioWrite(context, n, 4, -antiNoiseBlock[n]);
        scope.log(refBlock[n], errorBlock[n], -antiNoiseBlock[n], errorBlock[n]- antiNoiseBlock[n]);
    }

    currSample += blockSize;
}


void cleanup(BelaContext *context, void *userData)
{

}
