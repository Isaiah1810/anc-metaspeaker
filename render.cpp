#include <array>
#include <Bela.h>
#include <libraries/Fft/Fft.h>
#include <libraries/Scope/Scope.h>
#include <cmath>
#include <random>
#include <memory>

/** OSCILLOSCOPE LEGEND
 * Secondary Path Estimation:
 * red: error mic signal
 * blue: lms output (estimated secondary path with excitation signal)
 * green: excitation signal
 * 
 * Primary Path Estimation:
 * red: reference mic signal
 * blue: error mic signal
 * green: adaptive filter OUTPUT_CHANNEL
 * orange: lms error
 * **/


// Configuration
constexpr bool SIMULATION = true;
constexpr int IMPULSE_SIZE = 64;
constexpr int REFERENCE_CHANNEL = 0;
constexpr int ERROR_CHANNEL = 2;
constexpr int OUTPUT_CHANNEL = 4;
constexpr float STEP_SIZE = 0.0001f;
constexpr float CONVERGENCE_THRESHOLD = 0.001f;
constexpr float SMALL_VALUE = 1e-6f; // Avoid division by zero

// Class to handle Adaptive Noise Cancellation
class AdaptiveNoiseCanceller {
private:
    // Bela signal buffers 
    std::vector<float> referenceSignal;
    std::vector<float> errorSignal;
    
    // Filter coefficients
    std::array<float, IMPULSE_SIZE> primaryPath;
    std::array<float, IMPULSE_SIZE> secondaryPath;
    
    // Circular buffers
    std::array<float, IMPULSE_SIZE> referenceBuffer;
    std::array<float, IMPULSE_SIZE> excitationBuffer;
    std::array<float, IMPULSE_SIZE> outputBuffer;
    std::array<float, IMPULSE_SIZE> filteredReference;
    
    // Buffer indices
    int referenceIdx = 0;
    int excitationIdx = 0;
    int outputIdx = 0;
    
    // Simulation parameters
    std::array<float, IMPULSE_SIZE> primaryGroundTruth;
    std::array<float, IMPULSE_SIZE> secondaryGroundTruth;
    
    // Sine wave generation for simulation
    float frequency = 440.0f;
    float amplitude = 0.2f;
    float phase = 0.0f;
    float sampleRate = 44100.0f;
    
    // State tracking
    bool secondaryPathLearned = false;
    int totalSamples = 0;
    
    // Random number generator
    std::mt19937 rng;
    std::uniform_real_distribution<float> dist{-0.5f, 0.5f};
    
    Scope& scope;

public:
    AdaptiveNoiseCanceller(BelaContext* context, Scope& scopeRef) 
        : scope(scopeRef), rng(std::random_device{}()) {
        
        // Initialize buffers
        referenceSignal.resize(context->audioFrames, 0.0f);
        errorSignal.resize(context->audioFrames, 0.0f);
        
        // Initialize filters
        initializeGroundTruthFilters();
        initializeFiltersRandomly();
        
        // Set sample rate
        sampleRate = context->audioSampleRate;
    }
    
    // Generate white noise with specified amplitude
    inline float generateWhiteNoise(float amp) {
        return amp * dist(rng);
    }
    
    // Initialize filters with ground truth values for simulation
    void initializeGroundTruthFilters() {
        // Primary path: 10ms delay, low pass filter
        std::fill(primaryGroundTruth.begin(), primaryGroundTruth.end(), 0.0f);
        primaryGroundTruth[20] = 0.8f;
        primaryGroundTruth[25] = -0.4f;
        primaryGroundTruth[30] = 0.2f;
        
        // Secondary path: 5ms delay, low pass filter
        std::fill(secondaryGroundTruth.begin(), secondaryGroundTruth.end(), 0.0f);
        secondaryGroundTruth[10] = 0.6f;
        secondaryGroundTruth[12] = -0.3f;
        secondaryGroundTruth[15] = 0.1f;
    }
    
    // Initialize filters with small random values
    void initializeFiltersRandomly() {
        // Random distribution for filter initialization
        std::uniform_real_distribution<float> initDist{-0.1f, 0.1f};
        
        // Initialize filters with small random values
        for (int i = 0; i < IMPULSE_SIZE; ++i) {
            primaryPath[i] = initDist(rng);
            secondaryPath[i] = initDist(rng);
            outputBuffer[i] = 0.01f * initDist(rng);
            referenceBuffer[i] = 0.0f;
            excitationBuffer[i] = 0.0f;
            filteredReference[i] = 0.0f;
        }
    }
    
    // Process a block of audio samples
    void process(BelaContext* context) {
        // Read or generate input signals
        readInputSignals(context);
        
        // Choose algorithm based on state
        if (!secondaryPathLearned) {
            learnSecondaryPath(context);
        } else {
            applyFxNLMS(context);
        }
    }
    
    // Read input signals from microphones or generate simulated signals
    void readInputSignals(BelaContext* context) {
        if (!SIMULATION) {
            // Read microphone inputs
            for (int n = 0; n < context->audioFrames; ++n) {
                totalSamples++;
                referenceSignal[n] = audioRead(context, n, REFERENCE_CHANNEL);
                referenceBuffer[referenceIdx] = referenceSignal[n];
                referenceIdx = (referenceIdx + 1) % IMPULSE_SIZE;
                
                errorSignal[n] = audioRead(context, n, ERROR_CHANNEL);
            }
        } else {
            // Generate simulated inputs
            for (unsigned int n = 0; n < context->audioFrames; n++) {
                totalSamples++;
                
                // Generate sine wave reference
                phase += 2.0f * M_PI * frequency / sampleRate;
                if (phase > 2.0f * M_PI) {
                    phase -= 2.0f * M_PI;
                }
                
                referenceSignal[n] = amplitude * sinf(phase);
                referenceBuffer[referenceIdx] = referenceSignal[n];
                
                // Convolve reference with primary path and output with secondary path
                float primaryOutput = 0.0f;
                float secondaryOutput = 0.0f;
                
                // Convolve primary and secondary paths with reference and output respectivly 
                for (int k = 0; k < IMPULSE_SIZE; ++k) {
                    int refIdx = (referenceIdx - k + IMPULSE_SIZE) % IMPULSE_SIZE;
                    primaryOutput += primaryGroundTruth[k] * referenceBuffer[refIdx];
                    
                    int outIdx = (outputIdx - k - (context->audioFrames - n) + IMPULSE_SIZE) % IMPULSE_SIZE;
                    secondaryOutput += secondaryGroundTruth[k] * outputBuffer[outIdx];
                }
                
                referenceIdx = (referenceIdx + 1) % IMPULSE_SIZE;
                errorSignal[n] = primaryOutput + secondaryOutput;
            }
        }
    }
    
    // Learn the secondary path using white noise excitation
    void learnSecondaryPath(BelaContext* context) {
        float blockError = 0.0f;
        
        // Process each sample
        for (int n = 0; n < context->audioFrames; ++n) {
            // Convolve excitation with estimated secondary path
            float filterOutput = 0.0f;
            
            for (int k = 0; k < IMPULSE_SIZE; ++k) {
                int excitIdx = (excitationIdx - k - (context->audioFrames - n) + IMPULSE_SIZE) % IMPULSE_SIZE;
                filterOutput += excitationBuffer[excitIdx] * secondaryPath[k];
            }
            
            // Calculate LMS error
            float lmsError = errorSignal[n] - filterOutput;
            blockError += lmsError * lmsError;
            
            // Calculate excitation signal power for normalization
            float excitationPower = SMALL_VALUE;
            for (int k = 0; k < IMPULSE_SIZE; ++k) {
                int excitIdx = (excitationIdx - k - (context->audioFrames - n) + IMPULSE_SIZE) % IMPULSE_SIZE;
                float sample = excitationBuffer[excitIdx];
                excitationPower += sample * sample;
            }
            
            // Update filter coefficients using NLMS
            float normalizedStepSize = STEP_SIZE / excitationPower;
            for (int k = 0; k < IMPULSE_SIZE; ++k) {
                int excitIdx = (excitationIdx - k - (context->audioFrames - n) + IMPULSE_SIZE) % IMPULSE_SIZE;
                secondaryPath[k] += normalizedStepSize * lmsError * excitationBuffer[excitIdx];
            }
            
            // Log to scope for debugging
            int currentExcitIdx = (excitationIdx - (context->audioFrames - n) + IMPULSE_SIZE) % IMPULSE_SIZE;
            scope.log(errorSignal[n], filterOutput, excitationBuffer[currentExcitIdx]);
        }
        
        // Generate and output excitation signal
        for (int n = 0; n < context->audioFrames; ++n) {
            float excitation = generateWhiteNoise(0.2f);
            excitationBuffer[excitationIdx] = excitation;
            excitationIdx = (excitationIdx + 1) % IMPULSE_SIZE;
            
            // Output the excitation signal
            if (!SIMULATION) {
                audioWrite(context, n, OUTPUT_CHANNEL, excitation);
            } else {
                outputIdx = (outputIdx + 1) % IMPULSE_SIZE;
                outputBuffer[outputIdx] = excitation;
            }
        }
        
        // Check convergence
        if ((totalSamples > (context->audioFrames + 1)) && 
            (blockError / context->audioFrames) < CONVERGENCE_THRESHOLD) {
            secondaryPathLearned = true;
            rt_printf("Secondary path learned at sample %i\n", totalSamples);
        }
    }
    
    // Apply the Filtered-X Normalized LMS algorithm
    void applyFxNLMS(BelaContext* context) {
        // Pre-compute filtered reference signal through secondary path
        for (int n = 0; n < context->audioFrames; ++n) {
            float xFiltered = 0.0f;
            for (int k = 0; k < IMPULSE_SIZE; ++k) {
                int refIdx = (referenceIdx - k + IMPULSE_SIZE - (context->audioFrames - n)) % IMPULSE_SIZE;
                xFiltered += secondaryPath[k] * referenceBuffer[refIdx];
            }
            
            int filtIdx = (referenceIdx - (context->audioFrames - n) + IMPULSE_SIZE + 1) % IMPULSE_SIZE;
            filteredReference[filtIdx] = xFiltered;
        }
        
        // Process each sample with FxNLMS
        for (int n = 0; n < context->audioFrames; ++n) {
            // Convolve filtered reference with primary path
            float filterOutput = 0.0f;
            for (int k = 0; k < IMPULSE_SIZE; ++k) {
                int filtIdx = (referenceIdx - k - (context->audioFrames - n) + IMPULSE_SIZE) % IMPULSE_SIZE;
                filterOutput += primaryPath[k] * filteredReference[filtIdx];
            }
            
            float lmsError = errorSignal[n] - filterOutput;
            
            // Calculate power of filtered reference for normalization
            float filteredReferencePower = SMALL_VALUE;
            for (int k = 0; k < IMPULSE_SIZE; ++k) {
                int filtIdx = (referenceIdx - k - (context->audioFrames - n) + IMPULSE_SIZE) % IMPULSE_SIZE;
                float sample = filteredReference[filtIdx];
                filteredReferencePower += sample * sample;
            }
            
            // Update filter coefficients using FxNLMS
            float normalizedStepSize = STEP_SIZE / filteredReferencePower;
            for (int k = 0; k < IMPULSE_SIZE; ++k) {
                int filtIdx = (referenceIdx - k - (context->audioFrames - n) + IMPULSE_SIZE) % IMPULSE_SIZE;
                primaryPath[k] += normalizedStepSize * lmsError * filteredReference[filtIdx];
            }
            
            // Log to scope for debugging
            scope.log(referenceSignal[n], errorSignal[n], -filterOutput, lmsError);
            
            // Output the anti-noise
            if (!SIMULATION) {
                audioWrite(context, n, OUTPUT_CHANNEL, -filterOutput);
            } else {
                outputIdx = (outputIdx + 1) % IMPULSE_SIZE;
                outputBuffer[outputIdx] = -filterOutput;
            }
        }
    }
    
    // Print filter coefficients for debugging
    void printFilters() const {
        rt_printf("Secondary Path:\n");
        for (int k = 0; k < IMPULSE_SIZE; ++k) {
            rt_printf("%f\n", secondaryPath[k]);
        }
        
        rt_printf("Primary Path:\n");
        for (int k = 0; k < IMPULSE_SIZE; ++k) {
            rt_printf("%f\n", primaryPath[k]);
        }
    }
};

// Global objects
Scope scope;
AdaptiveNoiseCanceller* ancSystem = nullptr;

bool setup(BelaContext *context, void *userData)
{
    // Setup scope
    scope.setup(4, context->audioSampleRate);
    
    // Create ANC system
    ancSystem = new AdaptiveNoiseCanceller(context, scope);
    
    return true;
}

void render(BelaContext *context, void *userData)
{
    // Process audio block
    ancSystem->process(context);
}

void cleanup(BelaContext *context, void *userData)
{
    // Print final filter coefficients
    if (ancSystem) {
        ancSystem->printFilters();
        delete ancSystem;
        ancSystem = nullptr;
    }
}
