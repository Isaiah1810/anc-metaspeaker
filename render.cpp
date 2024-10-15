#include <memory>
#include <random>
#include <Bela.h>
#include <libraries/Fft/Fft.h>
#include <vector>
#include <libraries/Convolver/Convolver.h>
#include <cmath>

/** @todo 
 * -ACOUSTIC FEEDBACK (antiNoise->Reference needs to be estimated and subtracted)
 * I'm hoping for now it isn't an issue due to directionality of meta-speaker
 * -REDUCE AMOUNT OF VECTORS. Theres a bunch of vectors that could be shared to reduce size**/

const float carrierFreq = 40000;
const float carrierAmp = 0.8;

const float sampleRate = 44100;

const float stepSize = 0.2;

// Number of samples that have passed
unsigned int currSample;

// Define Audio Channels //
const unsigned int refChannel = 0;
const unsigned int errorChannel = 2;
const unsigned int speakerChannel = 4;
// Define microphone input vectors //
std::vector<float> refBlock;
std::vector<float> errorBlock;

// Define antiNoise(output) vector //
std::vector<float> antiNoiseBlock;

// Define training noise for secondary path //
std::vector<float> trainingNoiseBlock;

// Define filter constants //
unsigned int primaryFilterSize = 256;
unsigned int secondaryFilterSize = 256;
const float threshold = 0.001;

// Define Filters //
std::vector<float> primaryFilter;
std::vector<float> secondaryFilter;
std::vector<float> prevSecondaryFilter;

// Define Output Vector 
std::vector<float> output;
/** Boolean representing whether to train primary path and run noise control.
 *	is false when secondary path estimator filter is being trained */
bool doNoiseControl;

bool setup(BelaContext *context, void *userData)
{
	currSample = 0;

	//Initialize microphone input vectors to 0
	refBlock.resize(context->audioFrames, 0.0f);
	errorBlock.resize(context->audioFrames, 0.0f);

	//Initialize Primary and Secondary Path Filters to 
	primaryFilter.resize(primaryFilterSize, 0.0f);
	secondaryFilter.resize(secondaryFilterSize, 0.0f);
	prevSecondaryFilter.resize(secondaryFilterSize, 0.0f);

	//Initialize antiNoise vector
	antiNoiseBlock.resize(context->audioFrames, 0.0f);

	// Initialize training noise vector 
	trainingNoiseBlock.resize(context->audioFrames, 0.0f);

	// Initialize Output Vector
	output.resize(context->audioFrames, 0.0f);

	doNoiseControl = false;

	return true;
}

/** 
 * @brief Updates a filter using FxLMS, and applies the filter to reference signal. 
 * @param reference the audio block from the reference signal
 * @param error the audio block from the error mic
 * @param filter the filter coefficients
 * @param stepsSize step size for FxLMS
 * @param filterOrder the order of the filter
 * @param output the vector to which the output of the reference and the filter is stored
*/
void applyFxLMS(const std::vector<float> &reference,
				const std::vector<float> &error, std::vector<float> &filter,
				float stepSize, std::vector<float> &output, bool copyOut)
{
	//Ensure reference and error signal are same size
	if (reference.size() != error.size()) {
		throw std::invalid_argument("Reference and error signal must be of the same size.");
	}	

	// Set filter order to size of the filter
	size_t filterOrder = filter.size();

	// Initialize output to correct size
	size_t numSamples = reference.size();
	output.resize(numSamples);

	for (size_t n = 0; n < numSamples; ++n){
		// Convolve the filter with the reference signal
		float filterOutput = 0;
		for (size_t k = 0; k < filterOrder; ++k){
			if (n >= k){
				filterOutput += filter[k] * reference[n - k];
			}
		}

		// Get error value from error signal and filter output
		float errorValue = error[n] - filterOutput;

		// Update filter coefficients 
		for (size_t k = 0; k < filterOrder; ++k){
			if (n >= k){
				filter[k] += stepSize * errorValue * reference[n - k];
			}
		}
		// Store filter output
		if (copyOut)
			output[n] = filterOutput;
	}
	
}

/** 
 * @brief Generates training noise for use in secondary path estimation learning
 * @param noise Vector where noise is writen to
 * @param currSample The current elapsed sample from boot
*/
void generateNoise(std::vector<float> &noise, unsigned int currSample){
	std::default_random_engine generator;
	std::uniform_real_distribution<float> distribution(-1.0, 1.0);
	for (size_t n = 0; n < noise.size(); ++n){
		noise[n] = distribution(generator);
	}
	return;
}

/**
 * @brief Modulates audio to carrier frequency and removes negative baseband using hilber transform
 * @param context Bela context
 * @param antiNoise The audio to be processed
 * @param output A vector to where the processed audio is written
 * @param currSample The number of samples that have passed since the begining of the program
 * 
*/

// SHOULD BE HILBERT THEN MODULATE NOT VIS VERSA
void processAntiNoise(BelaContext *context, std::vector<float> &antiNoise, std::vector<float> &output, unsigned int currSample){

	static Fft fft(context->audioFrames);
	fft.setup(context->audioFrames);
	fft.fft(antiNoise);
	
	// Remove Negative Basebands with Hilbert Transform (Set buffers to zero, copy only first half)
	std::vector<float> reBuffer(context->audioFrames, 0.0f);
	std::vector<float> imBuffer(context->audioFrames, 0.0f);
	for (size_t n = 0; n < context->audioFrames / 2; n++) {
		reBuffer[n] = fft.fdr(n);
		imBuffer[n] = fft.fdi(n);
	} 
	
	fft.ifft(reBuffer, imBuffer);
	for (size_t n = 0; n < context->audioFrames; n++) {
		output[n] = fft.td(n);
	}

	// Modulate antiNoise by carrier 
	for (size_t n = 0; n < context->audioFrames; n++){
		float wave = sin(2 * M_PI * (n + currSample) * carrierFreq / context->audioSampleRate);
		output[n] = output[n] * wave;
	}
	return;
}

/// @brief Determines whether a given filter has converged using coefficient stability
/// @param filter The currently updated filter
/// @param prevFilter The previous filter
/// @param threshold The threshold for convergance 
/// @return True if absolutue difference is less than threshold, otherwise, false
bool checkConvergence(std::vector<float> &filter, std::vector<float> &prevFilter, float threshold){
	float delta = 0;
	for (size_t n = 0; n < filter.size(); ++n){
		delta += std::abs(filter[n] - prevFilter[n]);
	}
	return delta < threshold;
}


void render(BelaContext *context, void *userData)
{
	float ref;
	float error;
	
	//Read input from microphones and store them into vectors
	for (unsigned int n = 0; n < context->audioFrames; n++){
		// Read analog frames for both microphones
		ref = audioRead(context, n, refChannel);
		error = audioRead(context, n, errorChannel);
		errorBlock[n] = error;
		// Subtract out Antinoise from reference signal
		refBlock[n] = ref;
	}


	// Train secondary path filter if not converged
	if (!doNoiseControl){
		prevSecondaryFilter = secondaryFilter;
		generateNoise(trainingNoiseBlock, currSample);
		applyFxLMS(trainingNoiseBlock, errorBlock, secondaryFilter, stepSize, 
				   antiNoiseBlock, true);
		
		/** TODO: CHECK CONVERGANCE CONDITION */
		doNoiseControl = checkConvergence(secondaryFilter, prevSecondaryFilter, threshold);
	}
	// Perform Active Noise Control and FxLMS on primary path
	else{
		// Convolve reference signal with primary path 
		for (size_t n = 0; n < context->audioFrames; n++){
			float primarySignal = 0;
			for (size_t k = 0; k < primaryFilter.size(); k++){
				if (n >= k){
					primarySignal += primaryFilter[k] * refBlock[n - k];
				}
			}
			// Invert primary path signal 
			antiNoiseBlock[n] = -primarySignal;
		}

		// Convolve reference with secondary path
		for (size_t n = 0; n < context->audioFrames; n++){
			float secondarySignal = 0;
			for (size_t k = 0; k < secondaryFilter.size(); k++){
				if (n >= k){
					secondarySignal += secondaryFilter[k] * refBlock[n - k];
				}
			}
			// NOTE: Reusing trainingNoiseBlock since it does nothing on this branch
			trainingNoiseBlock[n] = secondarySignal;
		}	
		// Perform FxLMS to update weights (don't use output from this)	
		applyFxLMS(trainingNoiseBlock, errorBlock, primaryFilter, stepSize, output, false);

	}

	// Process output for speaker (hilbert+modulate)
	processAntiNoise(context, antiNoiseBlock, output, currSample);
	
	// Output antiNoise to Speaker
	for (size_t n = 0; n < context->audioFrames; n++){
		audioWrite(context, n, speakerChannel, output[n]);
	}

	currSample = currSample + context->audioFrames;
}

void cleanup(BelaContext *context, void *userData)
{

}
