#include <Bela.h>
#include <libraries/Fft/Fft.h>
#include <vector>

float carrier_freq = 40000;
float sample_rate = 44100;
const unsigned int block_size = 16; //This should be changed from the IDE

std::vector<float> *ref_mic = NULL;
std::vector<float> *error_mic = NULL;


bool setup(BelaContext *context, void *userData)
{
	
	return true;
}

void render(BelaContext *context, void *userData)
{

}

void cleanup(BelaContext *context, void *userData)
{

}
