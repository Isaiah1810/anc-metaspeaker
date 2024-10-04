#include <Bela.h>
#include <libraries/Fft/Fft.h>
#include <vector>

float carrierFreq = 40000;
float sampleRate = 44100;
const unsigned int blockSize = 16; //This should be changed from the IDE

int audio_frames_per_analog_frame;

/** Define Microphone Channels **/
const unsigned int refChannel = 0;
const unsigned int errorChannel = 2;


std::vector<float> *refMic = NULL;
std::vector<float> *errorMic = NULL;

bool setup(BelaContext *context, void *userData)
{
	refMic = new std::vector<float>;
	errorMic = new std::vector<float>;
	return true;
}

void render(BelaContext *context, void *userData)
{
//Get samples for current audio block
for (unsigned int n = 0; n < context->audioFrames; n++){
	

}


}

void cleanup(BelaContext *context, void *userData)
{

}
