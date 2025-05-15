#ifndef MOTOR_PLANT_SIM_H
#define MOTOR_PLANT_SIM_H

#ifdef __cplusplus
extern "C" {
#endif

// Initialize the motor plant simulation (resets internal state)
void MotorPlant_Init(void);

// Simulate one time step of the motor model
// Inputs:  Volt   - input voltage to the motor
// Outputs: theta  - position in radians
//          omega  - angular velocity in rad/s
//          current - motor current in Amperes
void MotorPlant_Step(float Volt, float* theta, float* omega, float* current);

float AddNoise(float y, float noise);

#ifdef __cplusplus
}
#endif

#endif // MOTOR_PLANT_SIM_H
