// LPFilter.h
#ifndef LPFILTER_H
#define LPFILTER_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

// Structure to store LPF state

typedef struct {
    float x_prev;  // Previous input sample
    float y_prev;  // Previous output sample
    float alpha;   // Coefficient for the filter
} LPFState;

/**
 * @brief Initializes a low-pass filter state with cutoff frequency and sample rate
 * @param state Pointer to LPFState structure
 * @param cutoff_freq Cutoff frequency in Hz
 * @param sample_rate Sampling rate in Hz
 */
void LPF_Init(LPFState *state, float cutoff_freq, float sample_rate);

/**
 * @brief Applies the low-pass filter to the current input sample
 * @param x_now Current input value
 * @param state Pointer to LPFState structure
 * @return Filtered output value
 */
float LPFilter(float x_now, LPFState *state);

/**
 * @brief Initializes all motor-related LPF states with default parameters
 */
void LPF_Init_Motor(void);

// External filter states (defined in LPFilter.c)
extern LPFState lpfState_pris_velocity;
extern LPFState lpfState_pris_current;
extern LPFState lpfState_rev_velocity;
extern LPFState lpfState_rev_current;

#ifdef __cplusplus
}
#endif

#endif // LPFILTER_H
