/* DC_motor_kalman.c - CMSIS DSP version */
#include "DC_motor_kalman.h"
#include "arm_math.h"

#define N 4   // State size
#define M 2   // Measurement size

// State vector and covariance
static arm_matrix_instance_f32 x;       // (N x 1)
static arm_matrix_instance_f32 P;       // (N x N)

// Matrices
static arm_matrix_instance_f32 A, B, C, Q, R;
static arm_matrix_instance_f32 At, Ct;
static arm_matrix_instance_f32 tempNxN, tempNxM, tempMxM;

static float x_data[N] = {0};
static float P_data[N * N] = {0};

static float A_data[N * N] = {
    1.0f, 9.02071e-5f, -1.0601f, 0.0633f,
    0.0f, -0.00443f, -1144.76f, 68.04f,
    0.0f, 0.0f, 1.0f, 0.0f,
    0.0f, -6.1875e-5f, 0.7307f, 0.9459f
};

static float B_data[N] = {3.33e-4f, 0.7213f, 0.0f, 0.0111f};
static float C_data[M * N] = {
    1, 0, 0, 0,
    0, 0, 0, 1
};

static float Q_data[N * N] = {
    0.0f,     0.0f,     0.0f,     0.0f,
    0.0f,     0.0f,     0.0f,     0.0f,
    0.0f,     0.0f,     1e-7f,     0.0f,
    0.0f,     0.0f,     0.0f,     0.0f
};

static float R_data[M * M] = {
    1.32e-9f, 0.0f,
    0.0f,     1.3e-6f
};

static float tempNxN_data[N * N];
static float tempNxM_data[N * M];
static float tempMxM_data[M * M];
static float K_data[N * M];
static float y_data[M];

static arm_matrix_instance_f32 K, y_vec;

static int errFlag = 0;

void Kalman_Init(void) {
    arm_mat_init_f32(&x, N, 1, x_data);
    arm_mat_init_f32(&P, N, N, P_data);
    arm_mat_init_f32(&A, N, N, A_data);
    arm_mat_init_f32(&B, N, 1, B_data);
    arm_mat_init_f32(&C, M, N, C_data);
    arm_mat_init_f32(&Q, N, N, Q_data);
    arm_mat_init_f32(&R, M, M, R_data);

    arm_mat_init_f32(&At, N, N, tempNxN_data);
    arm_mat_init_f32(&Ct, N, M, tempNxM_data);
    arm_mat_init_f32(&tempNxN, N, N, tempNxN_data);
    arm_mat_init_f32(&tempNxM, N, M, tempNxM_data);
    arm_mat_init_f32(&tempMxM, M, M, tempMxM_data);
    arm_mat_init_f32(&K, N, M, K_data);
    arm_mat_init_f32(&y_vec, M, 1, y_data);

    for (int i = 0; i < N; i++)
        P_data[i * N + i] = 1.0f;
}

void Kalman_Predict(float u) {
    // x = A*x + B*u
    arm_mat_mult_f32(&A, &x, &tempNxM); // A*x
    for (int i = 0; i < N; ++i)
        x_data[i] = tempNxM_data[i] + B_data[i] * u;

    // P = A*P*A^T + Q
    arm_mat_mult_f32(&A, &P, &tempNxN);
    arm_mat_trans_f32(&A, &At);
    arm_mat_mult_f32(&tempNxN, &At, &P);
    for (int i = 0; i < N * N; ++i)
        P_data[i] += Q_data[i];
}

void Kalman_Update(float y[2]) {
    // y - Cx
    arm_mat_trans_f32(&C, &Ct);
    arm_mat_mult_f32(&C, &x, &y_vec);
    y_vec.pData[0] = y[0] - y_vec.pData[0];
    y_vec.pData[1] = y[1] - y_vec.pData[1];

    // K = P*C^T * inv(C*P*C^T + R)
    arm_mat_mult_f32(&P, &Ct, &tempNxM); // P*C^T
    arm_mat_mult_f32(&C, &tempNxM, &tempMxM); // C*P*C^T
    for (int i = 0; i < M * M; ++i)
        tempMxM_data[i] += R_data[i];

    arm_status status = arm_mat_inverse_f32(&tempMxM, &tempMxM);
    if (status == ARM_MATH_SUCCESS)
        arm_mat_mult_f32(&tempNxM, &tempMxM, &K); // K = PC^T * inv(...)
    else
    	errFlag = 1;

    // x = x + K*(y - Cx)
    arm_mat_mult_f32(&K, &y_vec, &tempNxM);
    for (int i = 0; i < N; ++i)
        x_data[i] += tempNxM_data[i];

    // P = (I - K*C) * P
    float KC_data[N * N] = {0};
    arm_matrix_instance_f32 KC;
    arm_mat_init_f32(&KC, N, N, KC_data);
    arm_mat_mult_f32(&K, &C, &KC);
    for (int i = 0; i < N; ++i)
        KC_data[i * N + i] = 1.0f - KC_data[i * N + i]; // I - KC (diagonal approximation)
    arm_mat_mult_f32(&KC, &P, &P);
}

float* Kalman_GetState(void) {
    return x_data;
}
