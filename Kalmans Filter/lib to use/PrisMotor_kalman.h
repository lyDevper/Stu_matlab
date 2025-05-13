#ifndef PRISMOTOR_KALMAN_DSP_H
#define PRISMOTOR_KALMAN_DSP_H

#include "arm_math.h"
#include "Kalman_Data.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Initialize the Kalman filter for the prismatic motor using CMSIS-DSP.
 */
void PrisKalman_Init(void);

/**
 * @brief Perform the prediction step of the Kalman filter.
 *
 * @param u Input voltage to the motor.
 */
void PrisKalman_Predict(float u);

/**
 * @brief Perform the update step of the Kalman filter with new measurements.
 *
 * @param y Array of measurements: {position, current}.
 */
void PrisKalman_Update(float y[2]);

/**
 * @brief Get the current estimated state vector.
 *
 * @return Pointer to state array {position, velocity, load, current}.
 */
float* PrisKalman_GetState(void);

/**
 * @brief Perform one full Kalman filter step and return the MotorState struct.
 *
 * @param volt Input voltage to the motor.
 * @param position Measured position.
 * @param current Measured current.
 * @return MotorState Struct containing {position, velocity, load, current}.
 */
MotorState PrisKalman_Step(float volt, float position, float current);

/**
 * @brief Latest estimated motor states as struct (position, velocity, load, current).
 */
extern MotorState prisMotor_state;

#ifdef __cplusplus
}
#endif

#endif // PRISMOTOR_KALMAN_DSP_H
