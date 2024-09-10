#include <smmintrin.h> // For _mm_stream_load_si128
#include <emmintrin.h> // For _mm_mul_ps
#include <assert.h>
#include <stdint.h>

extern void saxpySerial(int N,
			float scale,
			float X[],
			float Y[],
			float result[]);


void saxpyStreaming(int N,
                    float scale,
                    float X[],
                    float Y[],
                    float result[])
{
        for (int i = 0; i < N; i += 8) {
            // Load inputs into registers
            __m128 x_vals = _mm_loadu_ps(&X[i]);
            __m128 y_vals = _mm_loadu_ps(&Y[i]);
            __m128 scale_vals = _mm_set1_ps(scale);

            // saxpy calculation
            __m128 result_vals = _mm_add_ps(_mm_mul_ps(scale_vals, x_vals), y_vals);

            // Write result directly to memory
            _mm_stream_ps(&result[i], result_vals);
        }

        // Make sure everything is committed before returning
        _mm_sfence();
}

