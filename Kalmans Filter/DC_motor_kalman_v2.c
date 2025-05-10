#include "DC_motor_kalman.h"
#include "arm_math.h"

#define N 4   // State size
#define M 2   // Measurement size

// State vector and covariance
static arm_matrix_instance_f32 x;       // (N x 1)
static arm_matrix_instance_f32 P;       // (N x N)

// Kalman matrices
static arm_matrix_instance_f32 A, B, C, Q, R;
static arm_matrix_instance_f32 At, Ct;
static arm_matrix_instance_f32 K, y_vec, KC;

// Temporary buffers
static arm_matrix_instance_f32 temp1;   // (N x 1) for A*x
static arm_matrix_instance_f32 temp2;   // (M x M) for S matrix
static arm_matrix_instance_f32 temp3;   // (N x N) for P update
static arm_matrix_instance_f32 temp4;   // (N x M) for P*C^T
static arm_matrix_instance_f32 temp5;   // (M x 1) for y - Cx
static arm_matrix_instance_f32 temp6;   // (M x 1) for K*(y - Cx)

// Buffers
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
    1e-5f, 0.0f, 0.0f, 0.0f,
    0.0f, 1e-5f, 0.0f, 0.0f,
    0.0f, 0.0f, 1e-5f, 0.0f,
    0.0f, 0.0f, 0.0f, 1e-5f
};
static float R_data[M * M] = {
    1.32e-9f, 0.0f,
    0.0f, 1.3e-6f
};

static float At_data[N * N];
static float Ct_data[N * M];
static float K_data[N * M];
static float y_data[M];
static float KC_data[N * N];
static float temp1_data[N];   // (N x 1)
static float temp2_data[M * M];   // (M x M)
static float temp3_data[N * N];   // (N x N)
static float temp4_data[N * M];   // (N x M)
static float temp5_data[M];   // (M x 1)
static float temp6_data[M];   // (M x 1)

static arm_matrix_instance_f32 Bu;
static float Bu_data[N];

static int errFlag = 0;

void Kalman_Init(void) {
    arm_mat_init_f32(&x, N, 1, x_data);
    arm_mat_init_f32(&P, N, N, P_data);
    arm_mat_init_f32(&A, N, N, A_data);
    arm_mat_init_f32(&B, N, 1, B_data);
    arm_mat_init_f32(&C, M, N, C_data);
    arm_mat_init_f32(&Q, N, N, Q_data);
    arm_mat_init_f32(&R, M, M, R_data);

    arm_mat_init_f32(&At, N, N, At_data);
    arm_mat_init_f32(&Ct, N, M, Ct_data);
    arm_mat_init_f32(&K, N, M, K_data);
    arm_mat_init_f32(&y_vec, M, 1, y_data);
    arm_mat_init_f32(&KC, N, N, KC_data);
    arm_mat_init_f32(&temp1, N, 1, temp1_data);  // (N x 1)
    arm_mat_init_f32(&temp2, M, M, temp2_data);  // (M x M)
    arm_mat_init_f32(&temp3, N, N, temp3_data);  // (N x N)
    arm_mat_init_f32(&temp4, N, M, temp4_data);  // (N x M)
    arm_mat_init_f32(&temp5, M, 1, temp5_data);  // (M x 1)
    arm_mat_init_f32(&temp6, M, 1, temp6_data);  // (M x 1)

    for (int i = 0; i < N * N; ++i) P_data[i] = 0.0f;
    for (int i = 0; i < N; ++i) P_data[i * N + i] = 1.0f;
}

void Kalman_Predict(float u) {
    // x = A*x + B*u
    arm_mat_mult_f32(&A, &x, &temp1);  // temp1 = A*x

    for (int i = 0; i < N; ++i) Bu_data[i] = B_data[i] * u;
    arm_mat_init_f32(&Bu, N, 1, Bu_data);
    arm_mat_add_f32(&temp1, &Bu, &x);  // x = temp1 + Bu

    // P = A*P*A^T + Q
    arm_mat_mult_f32(&A, &P, &temp3);  // temp3 = A*P
    arm_mat_trans_f32(&A, &At);  // At = A^T
    arm_mat_mult_f32(&temp3, &At, &P);  // P = temp3 * At
    arm_mat_add_f32(&P, &Q, &P);  // P = P + Q
}

void Kalman_Update(float y[2]) {
    // Innovation y - Cx
    arm_mat_trans_f32(&C, &Ct);  // Ct = C^T
    arm_mat_mult_f32(&C, &x, &y_vec);  // y_vec = C*x
    for (int i = 0; i < M; ++i) y_vec.pData[i] = y[i] - y_vec.pData[i];  // y_vec = y - Cx

    // S = C*P*C^T + R
    arm_mat_mult_f32(&P, &Ct, &temp4);  // temp4 = P*C^T
    arm_mat_mult_f32(&C, &temp4, &temp2);  // temp2 = C*temp4
    arm_mat_add_f32(&temp2, &R, &temp2);  // temp2 = temp2 + R (S matrix)

    // K = P*C^T * inv(S)
    if (arm_mat_inverse_f32(&temp2, &temp2) != ARM_MATH_SUCCESS) {
        errFlag = 1;
        return;
    }
    arm_mat_mult_f32(&temp4, &temp2, &K);  // K = temp4 * inv(S)

    // x = x + K*(y - Cx)
    arm_mat_mult_f32(&K, &y_vec, &temp6);  // temp6 = K*(y - Cx)
    arm_mat_add_f32(&x, &temp6, &x);  // x = x + temp6

    // P = (I - K*C)*P
    arm_mat_mult_f32(&K, &C, &KC);  // KC = K*C
    arm_mat_scale_f32(&KC, -1.0f, &KC);  // KC = -K*C
    for (int i = 0; i < N; ++i) KC.pData[i * N + i] += 1.0f;  // Add identity matrix
    arm_mat_mult_f32(&KC, &P, &P);  // P = KC * P
}

float* Kalman_GetState(void) {
    return x_data;
}
