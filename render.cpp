#include <array>
#include <Bela.h>
#include <libraries/Fft/Fft.h>
#include <libraries/Scope/Scope.h>
#include <cmath>
#include <random>

// Constants
const float carrierFreq = 40000.0f;
const float carrierAmp = 1.0f;
const float stepSize = 0.2f;
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

// Precomputed modulation buffer
std::vector<float> modulationBuffer;

// Function to initialize modulation buffer
void initializeModulationBuffer(size_t blockSize, float sampleRate) {
    modulationBuffer.resize(blockSize);
    for (size_t n = 0; n < blockSize; ++n) {
        modulationBuffer[n] = carrierAmp * sinf(2.0f * M_PI * carrierFreq * n / sampleRate);
    }
}

bool setup(BelaContext *context, void *userData) {
    size_t blockSize = context->audioFrames;
    refBlock.resize(blockSize, 0.0f);
    errorBlock.resize(blockSize, 0.0f);
    antiNoiseBlock.resize(blockSize, 0.0f);
    trainingNoiseBlock.resize(blockSize, 0.0f);
    output.resize(blockSize, 0.0f);
    initializeModulationBuffer(blockSize, context->audioSampleRate);

    fft.setup(blockSize);
    scope.setup(3, context->audioSampleRate);
    return true;
}

void generateNoise(std::vector<float> &noise, unsigned int currSample){
    static std::default_random_engine generator;
    static std::uniform_real_distribution<float> distribution(-1.0, 1.0);
    for (size_t n = 0; n < noise.size(); ++n){
        noise[n] = distribution(generator);
    }
    return;
}


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

void render(BelaContext *context, void *userData) {
    size_t blockSize = context->audioFrames;
    for (size_t n = 0; n < blockSize; ++n) {
        refBlock[n] = audioRead(context, n, 0);  // Reference mic
        errorBlock[n] = audioRead(context, n, 2); // Error mic
    }

    if (!doNoiseControl) {
        prevSecondaryFilter = secondaryFilter;
        generateNoise(trainingNoiseBlock, currSample);
        applyFxLMS(trainingNoiseBlock, errorBlock, secondaryFilter, stepSize, antiNoiseBlock, true);
        doNoiseControl = checkConvergence(secondaryFilter, prevSecondaryFilter);
    } else {
        applyFxLMS(refBlock, errorBlock, primaryFilter, stepSize, antiNoiseBlock, true);
    }

    processAntiNoise(antiNoiseBlock, output, context);
    for (size_t n = 0; n < blockSize; ++n) {
        audioWrite(context, n, 4, output[n]);
        scope.log(refBlock[n], errorBlock[n], antiNoiseBlock[n]);
    }

    currSample += blockSize;
}


void cleanup(BelaContext *context, void *userData)
{

}
