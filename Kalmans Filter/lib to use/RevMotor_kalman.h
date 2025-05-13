
#ifndef REVMOTOR_KALMAN_DSP_H
#define REVMOTOR_KALMAN_DSP_H

#include "arm_math.h"
#include "Kalman_Data.h"

#ifdef __cplusplus
extern "C" {
#endif

void RevKalman_Init(void);
void RevKalman_Predict(float u);
void RevKalman_Update(float y[2]);
float* RevKalman_GetState(void);
MotorState RevKalman_Step(float volt, float position, float current);

extern MotorState revMotor_state;

#ifdef __cplusplus
}
#endif

#endif // REVMOTOR_KALMAN_DSP_H
