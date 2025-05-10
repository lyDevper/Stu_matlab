/* DC_motor_kalman_nodsp.h */
#ifndef DC_MOTOR_KALMAN_NODSP_H
#define DC_MOTOR_KALMAN_NODSP_H

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

#endif // DC_MOTOR_KALMAN_NODSP_H
