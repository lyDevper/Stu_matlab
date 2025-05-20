// LPFilter.c
#include "LPFilter.h"
#include <math.h>

// Define PI for alpha calculation
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Initializes a filter state and calculates alpha coefficient
void LPF_Init(LPFState *state, float cutoff_freq, float sample_rate) {
    state->x_prev = 0.0f;
    state->y_prev = 0.0f;
    //   alpha = (2 * pi * f_c) / (2 * pi * f_c + f_s)
    state->alpha = (2.0f * M_PI * cutoff_freq) / (2.0f * M_PI * cutoff_freq + sample_rate);
}

// To be called in the main loop
float LPFilter(float x_now, LPFState *state) {
    //   y[n] = alpha * x[n] + alpha * x[n-1] + (1 - 2 * alpha) * y[n-1]
    float y = state->alpha * x_now
            + state->alpha * state->x_prev
            + (1.0f - 2.0f * state->alpha) * state->y_prev;

    state->x_prev = x_now;
    state->y_prev = y;

    return y;
}

// implement to motor sensors ------------------------
LPFState lpfState_pris_velocity;
LPFState lpfState_pris_current;

LPFState lpfState_rev_velocity;
LPFState lpfState_rev_current;

// To be called in main begin
void LPF_Init_Motor() {
    // sampling rate = 1000 Hz    
    // cutoff frequency = 25 Hz
    LPF_Init(&lpfState_pris_velocity, 15.0f, 1000.0f);
    LPF_Init(&lpfState_rev_velocity, 15.0f, 1000.0f);

    LPF_Init(&lpfState_pris_current, 50.0f, 1000.0f);
    LPF_Init(&lpfState_rev_current, 50.0f, 1000.0f);
}

/*
--------- example usage ------------
  pris_velocity_filtered = LPFilter(pris_velocity_raw, &lpfState_pris_velocity);
  rev_velocity_filtered = LPFilter(rev_velocity_raw, &lpfState_rev_velocity);
*/
