#include <memory>
#include <Bela.h>
#include <libraries/Fft/Fft.h>
#include <vector>
#include <libraries/Convolver/Convolver.h>
#include <cmath>

/** @todo 
 * -Make FxLMS use fft_convolve
 * -Be smart about timing for secondary path training
 * -Implement sacfolds	 **/

const float carrierFreq = 40000;
const float carrierAmp = 0.8;

const float sampleRate = 44100;

const float stepSize = 0.2;

int gAudioFramesPerAnalogFrame;

// Number of samples that have passed
unsigned int currSample;

// Define Microphone Channels //
const unsigned int refChannel = 0;
const unsigned int errorChannel = 2;

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

// Define Filters //
std::vector<float> primaryFilter;
std::vector<float> secondaryFilter;

/** Boolean representing whether to train primary path and run noise control.
 *	is false when secondary path estimator filter is being trained */
bool doNoiseControl;

bool setup(BelaContext *context, void *userData)
{
	gAudioFramesPerAnalogFrame = context->audioFrames / context->analogFrames;

	currSample = 0;

	//Initialize microphone input vectors to 0
	refBlock.resize(context->audioFrames, 0.0f);
	errorBlock.resize(context->audioFrames, 0.0f);

	//Initialize Primary and Secondary Path Filters to 
	primaryFilter.resize(primaryFilterSize, 0.0f);
	secondaryFilter.resize(secondaryFilterSize, 0.0f);
	
	//Initialize antiNoise vector
	antiNoiseBlock.resize(context->audioFrames, 0.0f);

	//Initialize training noise vector 
	trainingNoiseBlock.resize(context->audioFrames, 0.0f);

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
				float stepSize, std::vector<float> &output)
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
		output[n] = filterOutput;
	}
	
}

/** 
 * @brief Generates training noise for use in secondary path estimation learning
 * @param noise Vector where noise is writen to
 * @param currSample The current elapsed sample from boot
*/
void generateNoise(std::vector<float> &noise, unsigned int currSample);

/**
 * @brief Modulates audio to carrier frequency and removes negative baseband using hilber transform
 * @param context Bela context
 * @param antiNoise The audio to be processed
 * @param output A vector to where the processed audio is written
 * @param currSample The number of samples that have passed since the begining of the program
 * 
*/
void processAntiNoise(BelaContext *context, std::vector<float> &antiNoise, std::vector<float> &output, unsigned int currSample){
	// Modulate antiNoise by carrier 
	for (size_t n = 0; n < context->audioFrames; n++){
		float wave = sin(2 * M_PI * (n + currSample) * carrierFreq / context->audioSampleRate);
		output[n] = antiNoise[n] * wave;
	}
	static Fft fft(context->audioFrames);
	fft.setup(context->audioFrames);
	fft.fft(output);
	
	// Remove Negative Basebands with Hilbert Transform (Set buffers to zero, copy only first half)
	std::vector<float> reBuffer(context->audioFrames, 0.0f);
	std::vector<float> imBuffer(context->audioFrames, 0.0f);
	for (size_t n = 0; n < context->audioFrames / 2; n++) {
		reBuffer[n] = fft.fdr(n);
		imBuffer[n] = fft.fdi(n);
	} 
	
	fft.ifft(&reBuffer, &imBuffer);
	for (size_t n = 0; n < context->audioFrames; n++) {
		output[n] = fft.td(n);
	}

}


void render(BelaContext *context, void *userData)
{
	float ref;
	float error;
	
	//Read input from microphones and store them into vectors
	for (unsigned int n = 0; n < context->audioFrames; n++){
		if(gAudioFramesPerAnalogFrame && !(n % gAudioFramesPerAnalogFrame)) {
			//Read analog frames for both microphones
			ref = analogRead(context, n/gAudioFramesPerAnalogFrame, refChannel);
			error = analogRead(context, n/gAudioFramesPerAnalogFrame, errorChannel);
			refBlock[n] = ref;
			errorBlock[n] = error;
		}
	}

	// Train secondary path filter if not converged
	if (!doNoiseControl){
		generateNoise(&trainingNoiseBlock, currSample);
		applyFxLMS(&trainingNoiseBlock, &errorBlock, &secondaryFilter, stepSize, 
				   &antiNoiseBlock);
	}
	// Perform Active Noise Control and FxLMS on primary path
	else{
		// Convolve reference signal with primary path 
		
		// Invert primary path signal 

		// Convolve reference with secondary path

		// Compute error signal for FxLMS as errorBlock - (ref conv secondary)

		// Perform FxLMS to update weights (don't use output from this)

	}

	// Process output for speaker (hilbert+modulate)
 
	//
	currSample = currSample + context->audioFrames;
}

void cleanup(BelaContext *context, void *userData)
{

}
