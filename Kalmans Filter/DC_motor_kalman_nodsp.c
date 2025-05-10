/* DC_motor_kalman_nodsp.c - No DSP version */
#include "DC_motor_kalman.h"
#include <string.h> // for memset

#define N 4   // state size
#define M 2   // measurement size

static float x[N];          // state estimate
static float P[N][N];       // covariance

static float A[N][N] = {
    {1.0f, 9.02071e-5f, -1.0601f, 0.0633f},
    {0.0f, -0.00443f, -1144.76f, 68.04f},
    {0.0f, 0.0f, 1.0f, 0.0f},
    {0.0f, -6.1875e-5f, 0.7307f, 0.9459f}
};

static float B[N] = {3.33e-4f, 0.7213f, 0.0f, 0.0111f};

static float C[M][N] = {
    {1, 0, 0, 0},
    {0, 0, 0, 1}
};

static float Q[N][N]; // process noise
static float R[M][M]; // measurement noise

void Kalman_Init(void) {
    memset(x, 0, sizeof(x));
    memset(P, 0, sizeof(P));
    for (int i = 0; i < N; ++i) {
        P[i][i] = 1.0f;
        Q[i][i] = 1e-3f;
    }
    for (int i = 0; i < M; ++i)
        R[i][i] = 1e-2f;
}

void Kalman_Predict(float u) {
    float x_new[N] = {0};
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j)
            x_new[i] += A[i][j] * x[j];
        x_new[i] += B[i] * u;
    }
    memcpy(x, x_new, sizeof(x));

    float AP[N][N] = {0}, APA[N][N] = {0};
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            for (int k = 0; k < N; ++k)
                AP[i][j] += A[i][k] * P[k][j];

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            for (int k = 0; k < N; ++k)
                APA[i][j] += AP[i][k] * A[j][k];

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            P[i][j] = APA[i][j] + Q[i][j];
}

void Kalman_Update(float y[2]) {
    float y_est[M] = {0};
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < N; ++j)
            y_est[i] += C[i][j] * x[j];

    float y_diff[M];
    for (int i = 0; i < M; ++i)
        y_diff[i] = y[i] - y_est[i];

    float CP[N][M] = {0}, CPC[M][M] = {0};
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < M; ++j)
            for (int k = 0; k < N; ++k)
                CP[i][j] += P[i][k] * C[j][k];

    for (int i = 0; i < M; ++i)
        for (int j = 0; j < M; ++j)
            for (int k = 0; k < N; ++k)
                CPC[i][j] += C[i][k] * CP[k][j];

    for (int i = 0; i < M; ++i)
        for (int j = 0; j < M; ++j)
            CPC[i][j] += R[i][j];

    float det = CPC[0][0]*CPC[1][1] - CPC[0][1]*CPC[1][0];
    float invCPC[M][M];
    invCPC[0][0] =  CPC[1][1] / det;
    invCPC[1][1] =  CPC[0][0] / det;
    invCPC[0][1] = -CPC[0][1] / det;
    invCPC[1][0] = -CPC[1][0] / det;

    float K[N][M] = {0};
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < M; ++j)
            for (int k = 0; k < M; ++k)
                K[i][j] += CP[i][k] * invCPC[k][j];

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < M; ++j)
            x[i] += K[i][j] * y_diff[j];

    float KC[N][N] = {0};
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            for (int k = 0; k < M; ++k)
                KC[i][j] += K[i][k] * C[k][j];

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            P[i][j] -= KC[i][j] * P[j][j];
}

float* Kalman_GetState(void) {
    return x;
}
