#ifndef KALMAN_DATA_H
#define KALMAN_DATA_H

typedef struct {
    float position;  // rad or m (depends on context)
    float velocity;  // rad/s or m/s
    float load;      // Nm or N (disturbance torque/force)
    float current;   // A (motor current)
} MotorState;

#endif // KALMAN_DATA_H
