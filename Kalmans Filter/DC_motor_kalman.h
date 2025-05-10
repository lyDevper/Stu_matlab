/* DC_motor_kalman_dsp.h */
#ifndef DC_MOTOR_KALMAN_DSP_H
#define DC_MOTOR_KALMAN_DSP_H

#include "arm_math.h"

#ifdef __cplusplus
extern "C" {
#endif

void Kalman_Init(void);
void Kalman_Predict(float u);
void Kalman_Update(float y[2]);
float* Kalman_GetState(void);

#ifdef __cplusplus
}
#endif

#endif // DC_MOTOR_KALMAN_DSP_H
