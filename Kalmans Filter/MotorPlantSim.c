#include "MotorPlantSim.h"

#define N 4

// State vector
static float x_sim[N] = {0};

// Discrete-time state-space matrices (from Kalman filter model)
static const float A_sim[N][N] = {
    {1.0f, 9.02071e-5f, -1.0601f, 0.0633f},
    {0.0f, -0.00443f, -1144.76f, 68.04f},
    {0.0f, 0.0f, 1.0f, 0.0f},
    {0.0f, -6.1875e-5f, 0.7307f, 0.9459f}
};

static const float B_sim[N] = {3.33e-4f, 0.7213f, 0.0f, 0.0111f};

void MotorPlant_Init(void) {
    for (int i = 0; i < N; ++i)
        x_sim[i] = 0.0f;
}

void MotorPlant_Step(float Volt, float* theta, float* omega, float* current) {
    float x_new[N] = {0};

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            x_new[i] += A_sim[i][j] * x_sim[j];
        }
        x_new[i] += B_sim[i] * Volt;
    }

    for (int i = 0; i < N; ++i)
        x_sim[i] = x_new[i];

    if (theta) *theta = x_sim[0];
    if (omega) *omega = x_sim[1];
    if (current) *current = x_sim[3];
}
