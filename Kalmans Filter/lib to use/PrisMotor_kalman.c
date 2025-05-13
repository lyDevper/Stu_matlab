/* Kalman Filter for Prismatic DC Motor using CMSIS-DSP */
#include "PrisMotor_kalman.h"
#include "Kalman_Data.h"
#include "arm_math.h"

#define N 4 // State size
#define M 2 // Measurement size

// State and covariance
static arm_matrix_instance_f32 x;        // x: (N x 1)
static arm_matrix_instance_f32 P;        // P: (N x N)

// System matrices
static arm_matrix_instance_f32 A;        // A: (N x N)
static arm_matrix_instance_f32 B;        // B: (N x 1)
static arm_matrix_instance_f32 C;        // C: (M x N)
static arm_matrix_instance_f32 Q;        // Q: (N x N)
static arm_matrix_instance_f32 R;        // R: (M x M)

// Transposes and identity
static arm_matrix_instance_f32 At;       // At: (N x N)
static arm_matrix_instance_f32 Ct;       // Ct: (N x M)
static arm_matrix_instance_f32 I;        // I: (N x N)

// Kalman gain and intermediate values
static arm_matrix_instance_f32 K;        // K: (N x M)
static arm_matrix_instance_f32 Ax;       // Ax = A*x: (N x 1)
static arm_matrix_instance_f32 Bu;       // Bu = B*u: (N x 1)
static arm_matrix_instance_f32 AP;       // AP = A*P: (N x N)
static arm_matrix_instance_f32 APAt;     // APAt = AP*At: (N x N)
static arm_matrix_instance_f32 Cx;       // Cx = C*x: (M x 1)
static arm_matrix_instance_f32 y;        // y = measurement: (M x 1)
static arm_matrix_instance_f32 y_res;    // y_res = y - Cx: (M x 1)
static arm_matrix_instance_f32 CP;       // CP = C*P: (M x N)
static arm_matrix_instance_f32 CPCt;     // CPCt = CP*Ct: (M x M)
static arm_matrix_instance_f32 S;        // S = CPCt + R: (M x M)
static arm_matrix_instance_f32 S_inv;    // S_inv = inv(S): (M x M)
static arm_matrix_instance_f32 PCt;      // PCt = P*Ct: (N x M)
static arm_matrix_instance_f32 Ky_res;   // Ky_res = K*y_res: (N x 1)
static arm_matrix_instance_f32 KC;       // KC = K*C: (N x N)
static arm_matrix_instance_f32 IsKC;     // IsKC = I - KC: (N x N)

// Buffers
static float x_data[N] = {0};
static float P_data[N * N] = {0};
static float A_data[N * N] = {
    1.0f, 9.02071e-5f, -1.0601f, 0.0633f,
    0.0f, -0.00443f, -1144.76f, 68.04f,
    0.0f, 0.0f, 1.0f, 0.0f,
    0.0f, -6.1875e-5f, 0.7307f, 0.94591f
};
static float B_data[N] = {3.33e-4f, 0.7213f, 0.0f, 0.0111f};
static float C_data[M * N] = {
    1.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 0.0f, 1.0f
};
static float Q_data[N * N] = {
    0, 0, 0, 0,
    0, 0, 0, 0,
    0, 0, 1e-10f, 0,
    0, 0, 0, 0
};
static float R_data[M * M] = {
    1.32e-6f, 0,
    0, 1.3e-5f
};
static float I_data[N * N] = {
    1.0f, 0, 0, 0,
    0, 	1.0f, 0, 0,
    0, 	0, 1.0f, 0,
    0, 	0, 0, 1.0f
};
static float At_data[N * N], Ct_data[N * M];
static float Ax_data[N], Bu_data[N], AP_data[N * N], APAt_data[N * N];
static float Cx_data[M], y_data[M], y_res_data[M], CP_data[M * N];
static float CPCt_data[M * M], S_data[M * M], S_inv_data[M * M];
static float PCt_data[N * M], K_data[N * M], Ky_res_data[N];
static float KC_data[N * N], IsKC_data[N * N];

MotorState prisMotor_state;

void PrisKalman_Init() {
    arm_mat_init_f32(&x, N, 1, x_data);
    arm_mat_init_f32(&P, N, N, P_data);

    arm_mat_init_f32(&A, N, N, A_data);
    arm_mat_init_f32(&B, N, 1, B_data);
    arm_mat_init_f32(&C, M, N, C_data);
    arm_mat_init_f32(&Q, N, N, Q_data);
    arm_mat_init_f32(&R, M, M, R_data);

    arm_mat_init_f32(&At, N, N, At_data);
    arm_mat_trans_f32(&A, &At);

    arm_mat_init_f32(&Ct, N, M, Ct_data);
    arm_mat_trans_f32(&C, &Ct);

    arm_mat_init_f32(&I, N, N, I_data);

    arm_mat_init_f32(&Ax, N, 1, Ax_data);
    arm_mat_init_f32(&Bu, N, 1, Bu_data);
    arm_mat_init_f32(&AP, N, N, AP_data);
    arm_mat_init_f32(&APAt, N, N, APAt_data);

    arm_mat_init_f32(&Cx, M, 1, Cx_data);
    arm_mat_init_f32(&y, M, 1, y_data);
    arm_mat_init_f32(&y_res, M, 1, y_res_data);

    arm_mat_init_f32(&CP, M, N, CP_data);
    arm_mat_init_f32(&CPCt, M, M, CPCt_data);
    arm_mat_init_f32(&S, M, M, S_data);
    arm_mat_init_f32(&S_inv, M, M, S_inv_data);

    arm_mat_init_f32(&PCt, N, M, PCt_data);
    arm_mat_init_f32(&K, N, M, K_data);
    arm_mat_init_f32(&Ky_res, N, 1, Ky_res_data);
    arm_mat_init_f32(&KC, N, N, KC_data);
    arm_mat_init_f32(&IsKC, N, N, IsKC_data);
}

void PrisKalman_Predict(float u) {
    // x = A*x + B*u
    arm_mat_mult_f32(&A, &x, &Ax);
    arm_mat_scale_f32(&B, u, &Bu);
    arm_mat_add_f32(&Ax, &Bu, &x);

    // P = A*P*A^T + Q
    arm_mat_trans_f32(&A, &At);
    arm_mat_mult_f32(&A, &P, &AP);
    arm_mat_mult_f32(&AP, &At, &APAt);
    arm_mat_add_f32(&APAt, &Q, &P);
}

void PrisKalman_Update(float y_meas[2]) {
    memcpy(y_data, y_meas, sizeof(float) * M);

    // y_res = y - C*x
    arm_mat_mult_f32(&C, &x, &Cx);
    arm_mat_sub_f32(&y, &Cx, &y_res);

    // S = C*P*C^T + R
    arm_mat_trans_f32(&C, &Ct);
    arm_mat_mult_f32(&C, &P, &CP);
    arm_mat_mult_f32(&CP, &Ct, &CPCt);
    arm_mat_add_f32(&CPCt, &R, &S);

    // K = P*C^T * inv(S)
    arm_mat_mult_f32(&P, &Ct, &PCt);
    arm_mat_inverse_f32(&S, &S_inv);
    arm_mat_mult_f32(&PCt, &S_inv, &K);

    // x = x + K*y_res
    arm_mat_mult_f32(&K, &y_res, &Ky_res);
    arm_mat_add_f32(&x, &Ky_res, &x);

    // P = (I - K*C)*P
    arm_mat_mult_f32(&K, &C, &KC);
    arm_mat_sub_f32(&I, &KC, &IsKC);
    arm_mat_mult_f32(&IsKC, &P, &P);
}

float* PrisKalman_GetState() {
	prisMotor_state.position = x_data[0];
	prisMotor_state.velocity = x_data[1];
	prisMotor_state.load = x_data[2];
	prisMotor_state.current = x_data[3];

    return x_data;
}

MotorState PrisKalman_Step(float volt, float position, float current) {
	PrisKalman_Predict(volt);
	float y_meas[2] = {position, current};
	PrisKalman_Update(y_meas);
	PrisKalman_GetState();

	return prisMotor_state;
}
