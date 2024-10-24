// Copyright 2023 The gumy Authors. All rights reserved.
// Q运算
// Authors: gumy <ojbk@live.com>

#include "iqmath.h"

// sincos index
#define SIN_RAD     0x0300
#define U0_90       0x0000
#define U90_180     0x0100
#define U180_270    0x0200
#define U270_360    0x0300

const uint16_t SinCosTable[256] = {
    0x0000, 0x00c9, 0x0192, 0x025b, 0x0324, 0x03ed, 0x04b6, 0x057f, 0x0648,
    0x0711, 0x07d9, 0x08a2, 0x096b, 0x0a33, 0x0afb, 0x0bc4, 0x0c8c, 0x0d54,
    0x0e1c, 0x0ee4, 0x0fab, 0x1073, 0x113a, 0x1201, 0x12c8, 0x138f, 0x1455,
    0x151c, 0x15e2, 0x16a8, 0x176e, 0x1833, 0x18f9, 0x19be, 0x1a83, 0x1b47,
    0x1c0c, 0x1cd0, 0x1d93, 0x1e57, 0x1f1a, 0x1fdd, 0x209f, 0x2162, 0x2224,
    0x22e5, 0x23a7, 0x2467, 0x2528, 0x25e8, 0x26a8, 0x2768, 0x2827, 0x28e5,
    0x29a4, 0x2a62, 0x2b1f, 0x2bdc, 0x2c99, 0x2d55, 0x2e11, 0x2ecc, 0x2f87,
    0x3042, 0x30fc, 0x31b5, 0x326e, 0x3327, 0x33df, 0x3497, 0x354e, 0x3604,
    0x36ba, 0x3770, 0x3825, 0x38d9, 0x398d, 0x3a40, 0x3af3, 0x3ba5, 0x3c57,
    0x3d08, 0x3db8, 0x3e68, 0x3f17, 0x3fc6, 0x4074, 0x4121, 0x41ce, 0x427a,
    0x4326, 0x43d1, 0x447b, 0x4524, 0x45cd, 0x4675, 0x471d, 0x47c4, 0x486a,
    0x490f, 0x49b4, 0x4a58, 0x4afb, 0x4b9e, 0x4c40, 0x4ce1, 0x4d81, 0x4e21,
    0x4ec0, 0x4f5e, 0x4ffb, 0x5098, 0x5134, 0x51cf, 0x5269, 0x5303, 0x539b,
    0x5433, 0x54ca, 0x5560, 0x55f6, 0x568a, 0x571e, 0x57b1, 0x5843, 0x58d4,
    0x5964, 0x59f4, 0x5a82, 0x5b10, 0x5b9d, 0x5c29, 0x5cb4, 0x5d3e, 0x5dc8,
    0x5e50, 0x5ed7, 0x5f5e, 0x5fe4, 0x6068, 0x60ec, 0x616f, 0x61f1, 0x6272,
    0x62f2, 0x6371, 0x63ef, 0x646c, 0x64e9, 0x6564, 0x65de, 0x6657, 0x66d0,
    0x6747, 0x67bd, 0x6832, 0x68a7, 0x691a, 0x698c, 0x69fd, 0x6a6e, 0x6add,
    0x6b4b, 0x6bb8, 0x6c24, 0x6c8f, 0x6cf9, 0x6d62, 0x6dca, 0x6e31, 0x6e97,
    0x6efb, 0x6f5f, 0x6fc2, 0x7023, 0x7083, 0x70e3, 0x7141, 0x719e, 0x71fa,
    0x7255, 0x72af, 0x7308, 0x735f, 0x73b6, 0x740b, 0x7460, 0x74b3, 0x7505,
    0x7556, 0x75a6, 0x75f4, 0x7642, 0x768e, 0x76d9, 0x7723, 0x776c, 0x77b4,
    0x77fb, 0x7840, 0x7885, 0x78c8, 0x790a, 0x794a, 0x798a, 0x79c9, 0x7a06,
    0x7a42, 0x7a7d, 0x7ab7, 0x7aef, 0x7b27, 0x7b5d, 0x7b92, 0x7bc6, 0x7bf9,
    0x7c2a, 0x7c5a, 0x7c89, 0x7cb7, 0x7ce4, 0x7d0f, 0x7d3a, 0x7d63, 0x7d8a,
    0x7db1, 0x7dd6, 0x7dfb, 0x7e1e, 0x7e3f, 0x7e60, 0x7e7f, 0x7e9d, 0x7eba,
    0x7ed6, 0x7ef0, 0x7f0a, 0x7f22, 0x7f38, 0x7f4e, 0x7f62, 0x7f75, 0x7f87,
    0x7f98, 0x7fa7, 0x7fb5, 0x7fc2, 0x7fce, 0x7fd9, 0x7fe2, 0x7fea, 0x7ff1,
    0x7ff6, 0x7ffa, 0x7ffe, 0x7fff};

const int32_t AtanDiv[13] = {4096, 2418, 1278, 649, 326,
                             163, 81, 41, 20, 10, 5, 3, 1};

// IQSinCos 正余弦Q格式计算
void IQSinCos(IQSinCosParam *phasor, int32_t theta) {
    int32_t hindex = (theta >> 5); // 32768 / 32 => 1024 / 4 = 90° = 256
    int32_t quadrant = hindex & SIN_RAD;
    if (quadrant == U0_90) {
        phasor->Sine = SinCosTable[(uint8_t)(hindex)];
        phasor->Cosine = SinCosTable[(uint8_t)(0xFF - (uint8_t)(hindex))];
    } else if (quadrant == U90_180) {
        phasor->Sine = SinCosTable[(uint8_t)(0xFF - (uint8_t)(hindex))];
        phasor->Cosine = -SinCosTable[(uint8_t)(hindex)];
    } else if (quadrant == U180_270) {
        phasor->Sine = -SinCosTable[(uint8_t)(hindex)];
        phasor->Cosine = -SinCosTable[(uint8_t)(0xFF - (uint8_t)(hindex))];
    } else if (quadrant == U270_360) {
        phasor->Sine = -SinCosTable[(uint8_t)(0xFF - (uint8_t)(hindex))];
        phasor->Cosine = SinCosTable[(uint8_t)(hindex)];
    }

    return;
}

// IQSin 正弦查表法计算
int32_t IQSin(int32_t theta) {
    int32_t hindex = (theta >> 5); // 32768 / 32 => 1024 / 4 = 90° = 256
    int32_t quadrant = hindex & SIN_RAD;
    if (quadrant == U0_90) {
        return SinCosTable[(uint8_t)(hindex)];
    } else if (quadrant == U90_180) {
        return SinCosTable[(uint8_t)(0xFF - (uint8_t)(hindex))];
    } else if (quadrant == U180_270) {
        return -SinCosTable[(uint8_t)(hindex)];
    } else if (quadrant == U270_360) {
        return -SinCosTable[(uint8_t)(0xFF - (uint8_t)(hindex))];
    }

    return 0;
}

// IQCos 余弦查表法计算
int32_t IQCos(int32_t theta) {
    int32_t hindex = (theta >> 5); // 32768 / 32 => 1024 / 4 = 90° = 256
    int32_t quadrant = hindex & SIN_RAD;
    if (quadrant == U0_90) {
        return SinCosTable[(uint8_t)(0xFF - (uint8_t)(hindex))];
    } else if (quadrant == U90_180) {
        return -SinCosTable[(uint8_t)(hindex)];
    } else if (quadrant == U180_270) {
        return -SinCosTable[(uint8_t)(0xFF - (uint8_t)(hindex))];
    } else if (quadrant == U270_360) {
        return SinCosTable[(uint8_t)(hindex)];
    }

    return 0;
}

// IQAtan2 反正切Q格式计算
int32_t IQAtan2(int32_t y, int32_t x) {
    unsigned char quadrant;
    int32_t angle, newx;
    // 确认象限
    if ((x >= 0) && (y >= 0)) // 象限 0, 0 < angle < 90
    {
        quadrant = 0;
    } else if ((x < 0) && (y >= 0)) // 象限 1, 90 < angle < 180
    {
        quadrant = 1;
        x = -x;
    } else if ((x < 0) && (y < 0)) // 象限 2, 180 < angle < 270
    {
        quadrant = 2;
        x = -x;
        y = -y;
    } else if ((x >= 0) && (y < 0)) // 象限 3, 270 < angle < 360
    {
        quadrant = 3;
        y = -y;
    }
    angle = 0;

    for (int i = 0; i < 13; i++) {
        if (y < 0) {
            newx = x - (y >> i);
            y += (x >> i);
            x = newx;
            angle += AtanDiv[i];
        } else {
            newx = x + (y >> i);
            y -= (x >> i);
            x = newx;
            angle -= AtanDiv[i];
        }
    }

    // 0 ~ 2pi
    switch (quadrant) {
    case 0:
        angle = -angle;
        break;

    case 1:
        angle = _IQ15(0.5f) + angle;
        break;

    case 2:
        angle = _IQ15(0.5f) - angle;
        break;

    case 3:
        angle = _IQ15(1.0f) + angle;
        break;
    }

    return angle;
}

// IQSqrt 平方根计算
int32_t IQSqrt(uint32_t num) {
    uint32_t input = num;
    uint32_t sqrtX = 0;
    uint32_t temp = 0;

    for (int i = 30; i >= 0;) {
        temp = sqrtX + (1 << i);
        sqrtX >>= 1;
        if (temp <= input) {
            input -= temp;
            sqrtX += (1 << i);
        }
        i -= 2;
    }

    return sqrtX;
}
