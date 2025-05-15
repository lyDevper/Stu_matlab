
#include "RevMotor_kalman.h"
#include "Kalman_Data.h"
#include "arm_math.h"

#define N 4
#define M 2

// State and covariance
static arm_matrix_instance_f32 x, P;

// System matrices
static arm_matrix_instance_f32 A, B, C, Q, R;

// Transposes and identity
static arm_matrix_instance_f32 At, Ct, I;

// Kalman gain and intermediate variables
static arm_matrix_instance_f32 K, Ax, Bu, AP, APAt;
static arm_matrix_instance_f32 Cx, y, y_res, CP, CPCt, S, S_inv, PCt, Ky_res, KC, IsKC;

// Buffers
static float x_data[N] = {0};
static float P_data[N * N] = {0};

// including load inertia
static float A_data[N * N] = {
    1.0f, 0.0009995f, -0.00000293f, 0.00000139f,
    0.0f, 0.9989f, -0.00586f, 0.00275f,
    0.0f, 0.0f, 1.0f, 0.0f,
    0.0f, -0.01899f, 0.00005599f, 0.96201f
};

static float B_data[N] = {1.825e-8f, 5.4569e-5f, 0.0f, 3.8618e-2f};

// for y_measure = {theta, omega}
static float C_data[M * N] = {
    1.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 1.0f, 0.0f, 0.0f
};

/*
// for y_measure = {theta, current}
static float C_data[M * N] = {
    1.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 0.0f, 1.0f
};
*/

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
    0, 1.0f, 0, 0,
    0, 0, 1.0f, 0,
    0, 0, 0, 1.0f
};

static float At_data[N * N], Ct_data[N * M];
static float Ax_data[N], Bu_data[N], AP_data[N * N], APAt_data[N * N];
static float Cx_data[M], y_data[M], y_res_data[M], CP_data[M * N];
static float CPCt_data[M * M], S_data[M * M], S_inv_data[M * M];
static float PCt_data[N * M], K_data[N * M], Ky_res_data[N];
static float KC_data[N * N], IsKC_data[N * N];

MotorState revMotor_state;

void RevKalman_Init() {
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

void RevKalman_Predict(float u) {
    arm_mat_mult_f32(&A, &x, &Ax);
    arm_mat_scale_f32(&B, u, &Bu);
    arm_mat_add_f32(&Ax, &Bu, &x);

    arm_mat_trans_f32(&A, &At);
    arm_mat_mult_f32(&A, &P, &AP);
    arm_mat_mult_f32(&AP, &At, &APAt);
    arm_mat_add_f32(&APAt, &Q, &P);
}

void RevKalman_Update(float y_meas[2]) {
    memcpy(y_data, y_meas, sizeof(float) * M);

    arm_mat_mult_f32(&C, &x, &Cx);
    arm_mat_sub_f32(&y, &Cx, &y_res);

    arm_mat_trans_f32(&C, &Ct);
    arm_mat_mult_f32(&C, &P, &CP);
    arm_mat_mult_f32(&CP, &Ct, &CPCt);
    arm_mat_add_f32(&CPCt, &R, &S);

    arm_mat_mult_f32(&P, &Ct, &PCt);
    arm_mat_inverse_f32(&S, &S_inv);
    arm_mat_mult_f32(&PCt, &S_inv, &K);

    arm_mat_mult_f32(&K, &y_res, &Ky_res);
    arm_mat_add_f32(&x, &Ky_res, &x);

    arm_mat_mult_f32(&K, &C, &KC);
    arm_mat_sub_f32(&I, &KC, &IsKC);
    arm_mat_mult_f32(&IsKC, &P, &P);
}

float* RevKalman_GetState() {
    revMotor_state.position = x_data[0];
    revMotor_state.velocity = x_data[1];
    revMotor_state.load = x_data[2];
    revMotor_state.current = x_data[3];

    return x_data;
}

MotorState RevKalman_Step(float volt, float position, float velocity, float current) {
    RevKalman_Predict(volt);
    // use velocity instead of current
    float y_meas[2] = {position, velocity};
    RevKalman_Update(y_meas);
    RevKalman_GetState();

    return revMotor_state;
}
