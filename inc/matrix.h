#pragma once
#define _USE_MATH_DEFINES
#include <math.h>

namespace matrix
{
	//typedef float Mat4x4[16];
	using Mat4x4 = float[16];

	void Identity(Mat4x4 m)
	{
		m[0] = m[5] = m[10] = m[15] = 1;
		m[1] = m[2] = m[3] = m[4] = m[6] = m[7] = m[8] = m[9] = m[11] = m[12] = m[13] = m[14] = 0;
	}

	void Scale(Mat4x4 m, float x, float y, float z)
	{
		m[0] = x;
		m[1] = 0;
		m[2] = 0;
		m[3] = 0;

		m[4] = 0;
		m[5] = y;
		m[6] = 0;
		m[7] = 0;

		m[8] = 0;
		m[9] = 0;
		m[10] = z;
		m[11] = 0;

		m[12] = 0;
		m[13] = 0;
		m[14] = 0;
		m[15] = 1;
	}

	void Translate(float* m, float x, float y, float z)
	{
		m[0] = 1;
		m[1] = 0;
		m[2] = 0;
		m[3] = x;

		m[4] = 0;
		m[5] = 1;
		m[6] = 0;
		m[7] = y;

		m[8] = 0;
		m[9] = 0;
		m[10] = 1;
		m[11] = z;

		m[12] = 0;
		m[13] = 0;
		m[14] = 0;
		m[15] = 1;
	}

	void Rotate(float* m, double a, double x, double y, double z)
	{
		if (!a || (!x && !y && !z))
		{
			Identity(m);
			return;
		}

		auto d = sqrt((x * x) + (y * y) + (z * z));
		a *= M_PI / 180; x /= d; y /= d; z /= d;
		auto c = cos(a), s = sin(a), t = 1 - c;

		m[0] = x * x * t + c;
		m[1] = x * y * t - z * s;
		m[2] = x * z * t + y * s;
		m[3] = 0;

		m[4] = y * x * t + z * s;
		m[5] = y * y * t + c;
		m[6] = y * z * t - x * s;
		m[7] = 0;

		m[8] = z * x * t - y * s;
		m[9] = z * y * t + x * s;
		m[10] = z * z * t + c;
		m[11] = 0;

		m[12] = 0;
		m[13] = 0;
		m[14] = 0;
		m[15] = 1;
	}

	void Transpose(float* const m, float const* const src)
	{
		m[0] = src[0]; m[1] = src[4]; m[2] = src[8]; m[3] = src[12];
		m[4] = src[1]; m[5] = src[5]; m[6] = src[9]; m[7] = src[13];
		m[8] = src[2]; m[9] = src[6]; m[10] = src[10]; m[11] = src[14];
		m[12] = src[3]; m[13] = src[7]; m[14] = src[11]; m[15] = src[15];
	}

	void Multiply(float* const m, float const* const l, float const* const r)
	{
		m[0] = l[0] * r[0] + l[1] * r[4] + l[2] * r[8] + l[3] * r[12];
		m[1] = l[0] * r[1] + l[1] * r[5] + l[2] * r[9] + l[3] * r[13];
		m[2] = l[0] * r[2] + l[1] * r[6] + l[2] * r[10] + l[3] * r[14];
		m[3] = l[0] * r[3] + l[1] * r[7] + l[2] * r[11] + l[3] * r[15];

		m[4] = l[4] * r[0] + l[5] * r[4] + l[6] * r[8] + l[7] * r[12];
		m[5] = l[4] * r[1] + l[5] * r[5] + l[6] * r[9] + l[7] * r[13];
		m[6] = l[4] * r[2] + l[5] * r[6] + l[6] * r[10] + l[7] * r[14];
		m[7] = l[4] * r[3] + l[5] * r[7] + l[6] * r[11] + l[7] * r[15];

		m[8] = l[8] * r[0] + l[9] * r[4] + l[10] * r[8] + l[11] * r[12];
		m[9] = l[8] * r[1] + l[9] * r[5] + l[10] * r[9] + l[11] * r[13];
		m[10] = l[8] * r[2] + l[9] * r[6] + l[10] * r[10] + l[11] * r[14];
		m[11] = l[8] * r[3] + l[9] * r[7] + l[10] * r[11] + l[11] * r[15];

		m[12] = l[12] * r[0] + l[13] * r[4] + l[14] * r[8] + l[15] * r[12];
		m[13] = l[12] * r[1] + l[13] * r[5] + l[14] * r[9] + l[15] * r[13];
		m[14] = l[12] * r[2] + l[13] * r[6] + l[14] * r[10] + l[15] * r[14];
		m[15] = l[12] * r[3] + l[13] * r[7] + l[14] * r[11] + l[15] * r[15];
	}

	void Ortho(float* m, float width, float height, float near, float far)
	{
		float left = -width / 2, right = width / 2, bottom = -height / 2, top = height / 2;

		m[0] = 2 / (right - left);
		m[1] = 0;
		m[2] = 0;
		m[3] = -(right + left) / (right - left);

		m[4] = 0;
		m[5] = 2 / (top - bottom);
		m[6] = 0;
		m[7] = -(top + bottom) / (top - bottom);

		m[8] = 0;
		m[9] = 0;
		m[10] = -2 / (far - near);
		m[11] = -(far + near) / (far - near);

		m[12] = 0;
		m[13] = 0;
		m[14] = 0;
		m[15] = 1;
	}
}
