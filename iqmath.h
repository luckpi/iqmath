// Copyright 2023 The gumy Authors. All rights reserved.
// Q运算
// Authors: gumy <ojbk@live.com>

#ifndef IQMATH_H
#define IQMATH_H

#include "stdint.h"

#define _IQabs(A)           (int32_t)((A) < 0 ? -(A) : (A))
#define _IQsign(A)          (int32_t)((A) < 0 ? -1 : 1)
#define _IQmin(A, B)        (int32_t)((A) > (B) ? (B) : (A))
#define _IQmax(A, B)        (int32_t)((A) < (B) ? (B) : (A))
#define _IQsat(A, Max, Min) (int32_t)((A) > (Max) ? (Max) : ((A) < (Min) ? (Min) : (A)))

// Q格式转换
#define _IQ15(A)            (int32_t)((A) * 32768.0f)
#define _IQ15toF(A)         ((A) / 32768.0f)

#define _IQmpy(A, B, C)     (int32_t)(((A) * (B)) >> (C))
#define _IQdiv(A, B, C)     (int32_t)(((A) << (C)) / (B))
#define _IQdiv2(A)          (int32_t)((A) >> 1)
#define _IQmpy2(A)          (int32_t)((A) << 1)
#define _IQ15mpy(A, B)      (int32_t)(((A) * (B)) >> 15)
#define _IQ15div(A, B)      (int32_t)(((A) << 15) / (B))

// 2pi_pu
#define TWO_PI_PU           _IQ15(1.0f)
#define THETA_PU            TWO_PI_PU
// 1分钟 (s)
#define ONE_MIN             60.0f
// pi
#define ONE_PI              3.141592654f
// 2pi
#define TWO_PI              6.283185307f
// 1 / √3
#define ONE_BY_SQRT3        0.577350269f
// √3 / 2
#define SQRT3_BY_TWO        0.866025404f
// √3
#define SQRT3               1.732050808f

#define IQ_THETA_WARP(A)    (int32_t)((A) & 0x7FFF)

typedef struct
{
    int32_t Sine;
    int32_t Cosine;
} IQSinCosParam;

// IQSinCos 正余弦Q格式计算
void IQSinCos(IQSinCosParam *phasor, int32_t theta);

// IQSine 正弦Q格式计算
int32_t IQSin(int32_t theta);

// IQCos 余弦Q格式计算
int32_t IQCos(int32_t theta);

// IQAtan2 反正切Q格式计算
int32_t IQAtan2(int32_t y, int32_t x);

// IQSqrt 平方根Q格式计算
int32_t IQSqrt(uint32_t x);

#endif
