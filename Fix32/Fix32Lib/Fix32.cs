#define CHECKMATH
//#define USE_DOUBLES // Used for testing

using System;
using System.Runtime.CompilerServices;

/// <summary>
/// Represents a Q1.17.14 fixed-point number.
/// </summary>
public enum Fix32 : int {
	MaxValue = Fix32Ext.MAX_VALUE,
	MinValue = Fix32Ext.MIN_VALUE,

	Three = One * 3,
	Two = One * 2,
	One = Fix32Ext.ONE,
	Zero = 0,
	MinusOne = -One,
	MinusTwo = -One * 2,
	MinusThree = -One * 3,

	Pi = Fix32Ext.PI,
	PiOver2 = Fix32Ext.PI_OVER_2,
	PiOver4 = Fix32Ext.PI_OVER_4,
	PiTimes2 = Fix32Ext.PI_TIMES_2,
	PiInv = Fix32Ext.PI_INV,
	PiOver2Inv = Fix32Ext.PI_OVER_2_INV,
	E = Fix32Ext.E_RAW,
	EPow4 = Fix32Ext.EPOW4,
	Ln2 = Fix32Ext.LN2,
	Log2Max = Fix32Ext.LOG2MAX,
	Log2Min = Fix32Ext.LOG2MIN,
}

/// <summary>
/// Operations for <see cref="Fix32"/>.
/// </summary>
public static partial class Fix32Ext {

	// Precision of this type is 2^-14, that is 6.103515625E-5
	public static readonly decimal Precision = ((Fix32) 1).ToDecimal();

	internal const int MAX_VALUE = int.MaxValue;
	internal const int MIN_VALUE = int.MinValue;
	internal const int MAX_INT_VALUE = int.MaxValue >> FRACTIONAL_BITS;
	internal const int MIN_INT_VALUE = int.MinValue >> FRACTIONAL_BITS;
	public const int NUM_BITS = 32;
	public const int FRACTIONAL_BITS = 14;
	internal const int SIGNED_INTEGER_PLACES = NUM_BITS - FRACTIONAL_BITS;
	internal const int SIGN_MASK = unchecked((int) (1U << NUM_BITS_MINUS_ONE));

	internal const int NUM_BITS_MINUS_ONE = NUM_BITS - 1;
	internal const int ONE = 1 << FRACTIONAL_BITS;
	internal const int FRACTIONAL_MASK = ONE - 1;
	internal const int INTEGER_MASK = ~FRACTIONAL_MASK;
	internal const int LOG2MAX = (NUM_BITS - 1) << FRACTIONAL_BITS;
	internal const int LOG2MIN = -(NUM_BITS << FRACTIONAL_BITS);
	internal const int LUT_SIZE_RS = FRACTIONAL_BITS / 2 - 1;
	internal const int LUT_SIZE = PI_OVER_2 >> LUT_SIZE_RS;

	static readonly Fix32 LutInterval = (LUT_SIZE - 1).ToFix32().Div(Fix32.PiOver2);
	static readonly Fix32 C0p28 = 0.28m.ToFix32();

	// Const before rounding
	const decimal D_PI = ONE * (3.1415926535897932384626433832795028841971693993751058209749445923078164m);
	const decimal D_PI_TIMES_2 = ONE * (2m * 3.1415926535897932384626433832795028841971693993751058209749445923078164m);
	const decimal D_PI_OVER_2 = ONE * (3.1415926535897932384626433832795028841971693993751058209749445923078164m / 2m);
	const decimal D_PI_OVER_4 = ONE * (3.1415926535897932384626433832795028841971693993751058209749445923078164m / 4m);
	const decimal D_PI_INV = ONE * (1m / 3.1415926535897932384626433832795028841971693993751058209749445923078164m);
	const decimal D_PI_OVER_2_INV = ONE * (1m / (3.1415926535897932384626433832795028841971693993751058209749445923078164m / 2m));
	const decimal D_E_RAW = ONE * (2.71828182845904523536028747135266249775724709369995957496696762772407663035354759457138217852516642742746639193200m);
	const decimal D_EPOW4 = ONE * (2.71828182845904523536028747135266249775724709369995957496696762772407663035354759457138217852516642742746639193200m * 2.71828182845904523536028747135266249775724709369995957496696762772407663035354759457138217852516642742746639193200m * 2.71828182845904523536028747135266249775724709369995957496696762772407663035354759457138217852516642742746639193200m * 2.71828182845904523536028747135266249775724709369995957496696762772407663035354759457138217852516642742746639193200m);
	const decimal D_LN2 = ONE * (0.693147180559945309417232121458m);

	// Const rounded to int instead of truncated
	internal const int PI = (int) (D_PI < 0 ? D_PI - 0.5m : D_PI + 0.5m);
	internal const int PI_TIMES_2 = (int) (D_PI_TIMES_2 < 0 ? D_PI_TIMES_2 - 0.5m : D_PI_TIMES_2 + 0.5m);
	internal const int PI_OVER_2 = (int) (D_PI_OVER_2 < 0 ? D_PI_OVER_2 - 0.5m : D_PI_OVER_2 + 0.5m);
	internal const int PI_OVER_4 = (int) (D_PI_OVER_4 < 0 ? D_PI_OVER_4 - 0.5m : D_PI_OVER_4 + 0.5m);
	internal const int PI_INV = (int) (D_PI_INV < 0 ? D_PI_INV - 0.5m : D_PI_INV + 0.5m);
	internal const int PI_OVER_2_INV = (int) (D_PI_OVER_2_INV < 0 ? D_PI_OVER_2_INV - 0.5m : D_PI_OVER_2_INV + 0.5m);
	internal const int E_RAW = (int) (D_E_RAW < 0 ? D_E_RAW - 0.5m : D_E_RAW + 0.5m);
	internal const int EPOW4 = (int) (D_EPOW4 < 0 ? D_EPOW4 - 0.5m : D_EPOW4 + 0.5m);
	internal const int LN2 = (int) (D_LN2 < 0 ? D_LN2 - 0.5m : D_LN2 + 0.5m);



	/// <summary>
	/// Adds x and y. Performs saturating addition, i.e. in case of overflow, 
	/// rounds to MinValue or MaxValue depending on sign of operands.
	/// </summary>
	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	public static Fix32 Add(this Fix32 x, Fix32 y) {
#if USE_DOUBLES
		return (x.ToDouble() + y.ToDouble()).ToFix32();
#endif
		// https://stackoverflow.com/questions/17580118/signed-saturated-add-of-64-bit-ints/17587197#17587197
		// determine the lower or upper bound of the result
		//int ret = (x.RawValue < 0) ? MIN_VALUE : MAX_VALUE;
		int ret = (int) ((((uint) x >> NUM_BITS_MINUS_ONE) - 1U) ^ (1U << NUM_BITS_MINUS_ONE));
		// this is always well defined:
		// if x < 0 this adds a positive value to INT64_MIN
		// if x > 0 this subtracts a positive value from INT64_MAX
		//int comp = ret - xRaw;
		// the condition is equivalent to
		// ((x < 0) && (y > comp)) || ((x >=0) && (y <= comp))
		return (Fix32) (((int) x < 0) != ((int) y > (ret - (int) x)) ? ret : (int) x + (int) y);
	}

	/// <summary>
	/// Adds x and y. Doesn't saturate, wraps around when overflows.
	/// </summary>
	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	public static Fix32 AddFast(this Fix32 x, Fix32 y) {
		return (Fix32) ((int) x + (int) y);
	}

	/// <summary>
	/// Subtracts y from x. Performs saturating substraction, i.e. in case of overflow, 
	/// rounds to MinValue or MaxValue depending on sign of operands.
	/// </summary>
	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	public static Fix32 Sub(this Fix32 x, Fix32 y) {
#if USE_DOUBLES
		return (x.ToDouble() - y.ToDouble()).ToFix32();
#endif
		long sub = (long) x - (long) y; // TO TEST: Shift and operate to check overflow
		return (Fix32) (((int) sub) != sub ? (int) ((((uint) x >> NUM_BITS_MINUS_ONE) - 1U) ^ (1U << NUM_BITS_MINUS_ONE)) : (int) sub);
	}

	/// <summary>
	/// Subtracts y from x. Doesn't saturate, wraps around when overflows.
	/// </summary>
	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	public static Fix32 SubFast(this Fix32 x, Fix32 y) {
		return (Fix32) ((int) x - (int) y);
	}

	/// <summary>
	/// Returns a number indicating the sign of a Fix64 number.
	/// Returns 1 if the value is positive, 0 if is 0, and -1 if it is negative.
	/// </summary>
	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	public static int SignI(this Fix32 x) {
		//https://stackoverflow.com/questions/14579920/fast-sign-of-integer-in-c/14612418#14612418
		return ((int) x >> NUM_BITS_MINUS_ONE) | (int) (((uint) -(int) x) >> NUM_BITS_MINUS_ONE);
	}

	/// <summary>
	/// Returns a number indicating the sign of a Fix64 number.
	/// Returns 1 if the value is positive, 0 if is 0, and -1 if it is negative.
	/// </summary>
	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	public static Fix32 Sign(this Fix32 x) {
		const int RS = NUM_BITS_MINUS_ONE - FRACTIONAL_BITS;
		return (Fix32) ((((int) x >> RS) | (int) (((uint) -(int) x) >> RS)) & INTEGER_MASK);
	}


	/// <summary>
	/// Returns the absolute value of a Fix64 number.
	/// Note: Abs(Fix64.MinValue) == Fix64.MaxValue.
	/// </summary>
	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	public static Fix32 Abs(this Fix32 x) {
#if USE_DOUBLES
		return Math.Abs(x.ToDouble()).ToFix32();
#endif
		if ((int) x == MIN_VALUE) {
			return Fix32.MaxValue;
		}

		// branchless implementation, see http://www.strchr.com/optimized_abs_function
		var mask = (int) x >> NUM_BITS_MINUS_ONE;
		return (Fix32) (((int) x + mask) ^ mask);
	}

	/// <summary>
	/// Returns the absolute value of a Fix64 number.
	/// FastAbs(Fix64.MinValue) is undefined.
	/// </summary>
	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	public static Fix32 AbsFast(this Fix32 x) {
		var mask = (int) x >> NUM_BITS_MINUS_ONE;
		return (Fix32) (((int) x + mask) ^ mask);
	}

	/// <summary>
	/// Returns the largest integer less than or equal to the specified number.
	/// </summary>
	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	public static Fix32 Floor(this Fix32 x) {
#if USE_DOUBLES
		return Math.Floor(x.ToDouble()).ToFix32();
#endif
		// Just zero out the fractional part
		return (Fix32) ((int) x & INTEGER_MASK);
	}

	/// <summary>
	/// Returns the smallest integral value that is greater than or equal to the specified number.
	/// </summary>
	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	public static Fix32 Ceiling(this Fix32 x) {
#if USE_DOUBLES
		return Math.Floor(x.ToDouble()).ToFix32();
#endif
		var hasFractionalPart = ((int) x & FRACTIONAL_MASK) != 0;
		return hasFractionalPart ? ((Fix32) ((int) x & INTEGER_MASK)).Add(Fix32.One) : x;
	}

	/// <summary>
	/// Rounds a value to the nearest integral value.
	/// If the value is halfway between an even and an uneven value, returns the even value.
	/// </summary>
	public static Fix32 Round(this Fix32 x) {
#if USE_DOUBLES
		return Math.Round(v.ToDouble()).ToFix32();
#endif
		var fractionalPart = (int) x & FRACTIONAL_MASK;
		var integralPart = (Fix32) ((int) x & INTEGER_MASK);
		if (fractionalPart < (ONE >> 1)) {
			return integralPart;
		}
		if (fractionalPart > (ONE >> 1)) {
			return integralPart.Add(Fix32.One);
		}
		// if number is halfway between two values, round to the nearest even number
		// this is the method used by System.Math.Round().
		return ((int) integralPart & ONE) == 0 ?
			integralPart :
			integralPart.Add(Fix32.One);
	}

	/// <summary>
	/// Rounds a value to the nearest integral value.
	/// If the value is halfway between an even and an uneven value, returns the even value.
	/// FastRount(Fix64.MaxValue) is undefined
	/// </summary>
	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	public static Fix32 RoundFast(this Fix32 x) {
		// https://sestevenson.wordpress.com/2009/08/19/rounding-in-fixed-point-number-conversions/
		int odd = ((int) x & ONE) >> FRACTIONAL_BITS;
		return (Fix32) (((int) x + (ONE / 2 - 1) + odd) & INTEGER_MASK);
	}

	/// <summary>
	/// Multiply. Saturates when overflows
	/// </summary>
	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	public static Fix32 Mul(this Fix32 x, Fix32 y) {
#if USE_DOUBLES
		return Math.Round(x.ToDouble() * y.ToDouble()).ToFix32();
#endif
		long multLong = ((long) x * (long) y) >> FRACTIONAL_BITS;

		int finalSign = (int) x ^ (int) y;

		return (Fix32) (((finalSign ^ multLong) & SIGN_MASK) != 0 && multLong != 0 ?
			(int) ((((uint) finalSign >> NUM_BITS_MINUS_ONE) - 1U) ^ (1U << NUM_BITS_MINUS_ONE)) :
			(int) multLong);
	}

	/// <summary>
	/// Multiply. No saturation (overflows)
	/// </summary>
	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	public static Fix32 MulFast(this Fix32 x, Fix32 y) {
		return (Fix32) (((long) x * (long) y) >> FRACTIONAL_BITS);
	}

	/// <summary>
	/// Divide. Saturates on overflow
	/// </summary>
	public static Fix32 Div(this Fix32 x, Fix32 y) {
#if USE_DOUBLES
		return (x.ToDouble() / y.ToDouble()).ToFix32();
#endif
		long xl = (int) x;

		if ((int) y == 0) {
			return (Fix32) (unchecked((int) (((((uint) xl) >> NUM_BITS_MINUS_ONE) - 1U) ^ (1U << NUM_BITS_MINUS_ONE))));
			return xl >= 0 ? Fix32.MaxValue : Fix32.MinValue; // Branched version of the previous code, for clarity. Slower
		}

		long a = xl << FRACTIONAL_BITS;
		long b = (int) y;
		long r = a / b;
		if (r > MAX_VALUE) return Fix32.MaxValue;
		if (r < MIN_VALUE) return Fix32.MinValue;
		return (Fix32) (int) r;
	}

	/// <summary>
	/// Divide. Copied from DOOM. Wrong?
	/// </summary>
	public static Fix32 DivFast(this Fix32 x, Fix32 y) {
		//x = Fix32.MinValue;
		//y = Fix32.MinusOne;
		int xAbs = (int) x.Abs();
		int yAbs = (int) y.Abs();
		if ((xAbs >> 17/*((NUM_BITS - FRACTIONAL_PLACES) - 2)*/) >= yAbs)
			return ((int) x ^ (int) y) < 0 ? Fix32.MinValue : Fix32.MaxValue;
		return x.DivFastOverflow(y);
	}

	/// <summary>
	/// Divide. Copied from DOOM. Wrong?
	/// </summary>
	public static Fix32 DivFastOverflow(this Fix32 x, Fix32 y) {
#if true
		long c = (((long) x) << FRACTIONAL_BITS) / (long) y;
		return (Fix32) c;
#endif
		/*
		// From DOOM
		double c = ((double) x) / ((double) y) * (int) Fix32.One;

		if (c >= 2147483648.0 || c < -2147483648.0)
			I_Error("FixedDiv: divide by zero");
		return c.ToFix32();
		*/
	}

	/// <summary>
	/// Returns the base-2 logarithm of a specified number.
	/// </summary>
	/// <exception cref="ArgumentOutOfRangeException">
	/// The argument was non-positive
	/// </exception>
	public static Fix32 Log2(this Fix32 x) {
#if USE_DOUBLES
		return Math.Log(x.ToDouble(), 2).ToFix32();
#endif
		if ((int) x <= 0)
			throw new ArgumentOutOfRangeException("Non-positive value passed to Ln", "x");

		// This implementation is based on Clay. S. Turner's fast binary logarithm
		// algorithm (C. S. Turner,  "A Fast Binary Logarithm Algorithm", IEEE Signal
		//     Processing Mag., pp. 124,140, Sep. 2010.)

		//https://github.com/dmoulding/log2fix/blob/master/log2fix.c

		int b = 1 << (FRACTIONAL_BITS - 1);
		int y = 0;

		int rawX = (int) x;
		while (rawX < 1 << FRACTIONAL_BITS) {
			rawX <<= 1;
			y -= 1 << FRACTIONAL_BITS;
		}

		while (rawX >= 2 << FRACTIONAL_BITS) {
			rawX >>= 1;
			y += 1 << FRACTIONAL_BITS;
		}

		ulong z = (ulong) rawX;

		for (int i = 0; i < FRACTIONAL_BITS; i++) {
			z = z * z >> FRACTIONAL_BITS;
			if (z >= 2 << FRACTIONAL_BITS) {
				z >>= 1;
				y += b;
			}
			b >>= 1;
		}

		return (Fix32) (int) y;
	}

	/// <summary>
	/// Returns the natural logarithm of a specified number.
	/// </summary>
	/// <exception cref="ArgumentOutOfRangeException">
	/// The argument was non-positive
	/// </exception>
	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	public static Fix32 Ln(this Fix32 x) {
#if USE_DOUBLES
		return Math.Log(x.ToDouble()).ToFix32();
#endif
		return Log2(x).Mul(Fix32.Ln2);
	}

	/// <summary>
	/// Returns the square root of a specified number.
	/// </summary>
	/// <exception cref="ArgumentOutOfRangeException">
	/// The argument was negative.
	/// </exception>
	public static Fix32 Sqrt(this Fix32 x) {
#if USE_DOUBLES
		return Math.Sqrt(x.ToDouble()).ToFix32();
#endif
		if ((int) x < 0) {
			// We cannot represent infinities like Single and Double, and Sqrt is
			// mathematically undefined for x < 0. So we just throw an exception.
			throw new ArgumentOutOfRangeException("Negative value passed to Sqrt", "x");
		}

		/*
		uint t, q, b, r;
		r = (uint) x;
		b = 0x40000000;
		q = 0;
		while (b > 0x40) {
			t = q + b;
			if (r >= t) {
				r -= t;
				q = t + b; // equivalent to q += 2*b
			}
			r <<= 1;
			b >>= 1;
		}
		q >>= 8;
		return (Fix32) q;
		*/

		long xx = (long) x << FRACTIONAL_BITS;
		if (xx <= 1) return x;

		int s = 1;
		long x1 = xx - 1;
		if (x1 > 4294967295) { s += 16; x1 >>= 32; }
		if (x1 > 65535) { s += 8; x1 >>= 16; }
		if (x1 > 255) { s += 4; x1 >>= 8; }
		if (x1 > 15) { s += 2; x1 >>= 4; }
		if (x1 > 3) { s += 1; }
		long g0 = 1L << s;
		long g1 = (g0 + (xx >> s)) >> 1;
		while (g1 < g0) {
			g0 = g1;
			g1 = (g0 + (xx / g0)) >> 1;
		}
		return (Fix32) (g0);

		/*
		var xl = (int) x;

		var num = (uint) xl;
		var result = 0U;

		// second-to-top bit
		var bit = 1U << (NUM_BITS - 2);

		while (bit > num) {
			bit >>= 2;
		}

		// The main part is executed twice, in order to avoid
		// using 128 bit values in computations.
		for (var i = 0; i < 2; ++i) {
			// First we get the top 48 bits of the answer.
			while (bit != 0) {
				if (num >= result + bit) {
					num -= result + bit;
					result = (result >> 1) + bit;
				}
				else {
					result = result >> 1;
				}
				bit >>= 2;
			}

			if (i == 0) {
				// Then process it again to get the lowest 16 bits.
				if (num > (1U << (NUM_BITS / 2)) - 1) {
					// The remainder 'num' is too large to be shifted left
					// by 32, so we have to add 1 to result manually and
					// adjust 'num' accordingly.
					// num = a - (result + 0.5)^2
					//       = num + result^2 - (result + 0.5)^2
					//       = num - result - 0.5
					num -= result;
					num = (num << (NUM_BITS / 2)) - (1U << NUM_BITS_MINUS_ONE);
					result = (result << (NUM_BITS / 2)) + (1U << NUM_BITS_MINUS_ONE);
				}
				else {
					num <<= (NUM_BITS / 2);
					result <<= (NUM_BITS / 2);
				}

				bit = 1U << (NUM_BITS / 2 - 2);
			}
		}
		// Finally, if next bit would have been 1, round the result upwards.
		if (num > result) {
			++result;
		}
		return (Fix32) (int) result;
		*/
	}

	/// <summary>
	/// Returns 2 raised to the specified power.
	/// </summary>
	public static Fix32 Pow2(this Fix32 xx) {
#if USE_DOUBLES
		return Math.Pow(2, x.ToDouble()).ToFix32();
#endif
		Fix32 x = xx;
		if ((int) x == 0) return Fix32.One;

		// Avoid negative arguments by exploiting that exp(-x) = 1/exp(x).
		bool neg = (int) x < 0;
		if (neg) x = x.Neg();

		if ((int) x == (int) Fix32.One)
			return neg ? (Fix32) Fix32.One.Div(Fix32.Two) : (Fix32) Fix32.Two; // Can be cached
		if ((int) x >= (int) Fix32.Log2Max) return neg ? Fix32.One.Div(Fix32.MaxValue) : Fix32.MaxValue; // Can be cached
		if ((int) x <= (int) Fix32.Log2Min) return neg ? Fix32.MaxValue : Fix32.Zero;

		/* The algorithm is based on the power series for exp(x):
            * http://en.wikipedia.org/wiki/Exponential_function#Formal_definition
            * 
            * From term n, we get term n+1 by multiplying with x/n.
            * When the sum term drops to zero, we can stop summing.
            */
		int integerPart = (int) ((Fix32) x).Floor();
		// Take fractional part of exponent
		x = (Fix32) ((uint) x & FRACTIONAL_MASK);

		Fix32 result = Fix32.One;
		Fix32 term = Fix32.One;
		int i = 1;
		while ((int) term != 0) {
			term = x.Mul(term).Mul(Fix32.Ln2.Div(i.ToFix32()));
			result = result.Add(term);
			i++;
		}

		result = (Fix32) ((int) result << integerPart);
		if (neg) result = Fix32.One.Div(result);

		return result;
	}

	/// <summary>
	/// Returns a specified number raised to the specified power.
	/// </summary>
	/// <exception cref="DivideByZeroException">
	/// The base was zero, with a negative exponent
	/// </exception>
	/// <exception cref="ArgumentOutOfRangeException">
	/// The base was negative, with a non-zero exponent
	/// </exception>
	public static Fix32 Pow(this Fix32 b, Fix32 exp) {
#if USE_DOUBLES
		return Math.Pow(b.ToDouble(), exp.ToDouble()).ToFix32();
#endif
		if ((int) b == (int) Fix32.One)
			return Fix32.One;
		if ((int) exp == 0)
			return Fix32.One;
		if ((int) b == 0) {
			if ((int) exp < 0) {
				throw new DivideByZeroException();
			}
			return Fix32.Zero;
		}

		Fix32 log2 = b.Log2();
		return exp.Mul(log2).Pow2();
	}

	/// <summary>
	/// Returns the arccos of of the specified number, calculated using Atan and Sqrt
	/// </summary>
	public static Fix32 Acos(this Fix32 x) {
#if USE_DOUBLES
        return (Fix64) Math.Acos((double) x);
#endif
		if ((int) x < -(int) Fix32.One || (int) x > (int) Fix32.One)
			throw new ArgumentOutOfRangeException(nameof(x));

		if ((int) x == 0)
			return Fix32.PiOver2;

		var result = (Fix32.One.Sub(x.Mul(x))).Sqrt().Div(x).Atan();
		return (int) x < 0 ? result.Add(Fix32.Pi) : result;
	}

	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	static int CountLeadingZeroes(ulong x) {
		int result = 0;
		while ((x & 0xF000000000000000) == 0) { result += 4; x <<= 4; }
		while ((x & 0x8000000000000000) == 0) { result += 1; x <<= 1; }
		return result;
	}

	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	public static Fix32 Mod(this Fix32 x, Fix32 y) {
#if USE_DOUBLES
		return (x.ToDouble() % y.ToDouble()).ToFix32();
#endif
		return (Fix32) (
			(int) x == MIN_VALUE & (int) y == -1 ?
			0 :
			(int) x % (int) y);
	}

	/// <summary>
	/// Performs modulo as fast as possible; throws if x == MinValue and y == -1.
	/// Use the operator (%) for a more reliable but slower modulo.
	/// </summary>
	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	public static Fix32 ModFast(this Fix32 x, Fix32 y) {
		return (Fix32) ((int) x % (int) y);
	}

	/// <summary>
	/// Negate
	/// </summary>
	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	public static Fix32 Neg(this Fix32 x) {
#if USE_DOUBLES
		return (-x.ToDouble()).ToFix32();
#endif
		//return new Fix64(-x.RawValue);
		return (Fix32) ((int) x == MIN_VALUE ? MAX_VALUE : -(int) x);
	}

	/// <summary>
	/// Returns the Sine of x.
	/// The relative error is less than 1E-10 for x in [-2PI, 2PI], and less than 1E-7 in the worst case.
	/// </summary>
	public static Fix32 Sin(this Fix32 x) {
#if USE_DOUBLES
        return (Fix64) Math.Sin((double) x);
#endif
		var clampedRaw = ClampSinValue((int) x, out var flipHorizontal, out var flipVertical);
		var clamped = (Fix32) ((int) clampedRaw);

		// Find the two closest values in the LUT and perform linear interpolation
		// This is what kills the performance of this function on x86 - x64 is fine though
		Fix32 rawIndex = clamped.Mul(LutInterval);
		Fix32 roundedIndex = rawIndex.Round();
		Fix32 indexError = rawIndex.Sub(roundedIndex);

		var nearestValue = (Fix32) ((int) SinLut[flipHorizontal ?
			SinLut.Length - 1 - (int) roundedIndex :
			(int) roundedIndex]);
		var secondNearestValue = (Fix32) ((int) SinLut[flipHorizontal ?
			SinLut.Length - 1 - (int) roundedIndex - SignI(indexError) :
			(int) roundedIndex + SignI(indexError)]);

		int delta = (int) indexError.Mul(nearestValue.Sub(secondNearestValue).Abs());
		var interpolatedValue = (int) nearestValue + (flipHorizontal ? -delta : delta);
		var finalValue = flipVertical ? -interpolatedValue : interpolatedValue;
		return (Fix32) (finalValue);
	}

	/// <summary>
	/// Returns a rough approximation of the Sine of x.
	/// This is at least 3 times faster than Sin() on x86 and slightly faster than Math.Sin(),
	/// however its accuracy is limited to 4-5 decimals, for small enough values of x.
	/// </summary>
	public static Fix32 SinFast(this Fix32 x) {
		var clampedL = ClampSinValue((int) x, out bool flipHorizontal, out bool flipVertical);

		// Here we use the fact that the SinLut table has a number of entries
		// equal to (PI_OVER_2 >> LUT_SIZE_RS) to use the angle to index directly into it
		var rawIndex = (uint) (clampedL >> LUT_SIZE_RS);
		if (rawIndex >= LUT_SIZE) {
			rawIndex = LUT_SIZE - 1;
		}
		var nearestValue = SinLut[flipHorizontal ?
			SinLut.Length - 1 - (int) rawIndex :
			(int) rawIndex];
		return (Fix32) (int) (flipVertical ? -nearestValue : nearestValue);
	}

	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	static long ClampSinValue(long angle, out bool flipHorizontal, out bool flipVertical) {
#if CHECKMATH
		var clamped2Pi = angle;
		for (int i = 0; i < LARGE_PI_TIMES; ++i) {
			clamped2Pi %= (LARGE_PI_RAW >> i);
		}
#else
		// Clamp value to 0 - 2*PI using modulo; this is very slow but there's no better way AFAIK
		var clamped2Pi = angle % PI_TIMES_2;
#endif
		if (angle < 0) {
			clamped2Pi += PI_TIMES_2;
		}

		// The LUT contains values for 0 - PiOver2; every other value must be obtained by
		// vertical or horizontal mirroring
		flipVertical = clamped2Pi >= PI;
		// obtain (angle % PI) from (angle % 2PI) - much faster than doing another modulo
		var clampedPi = clamped2Pi;
		while (clampedPi >= PI) {
			clampedPi -= PI;
		}
		flipHorizontal = clampedPi >= PI_OVER_2;
		// obtain (angle % PI_OVER_2) from (angle % PI) - much faster than doing another modulo
		var clampedPiOver2 = clampedPi;
		if (clampedPiOver2 >= PI_OVER_2) {
			clampedPiOver2 -= PI_OVER_2;
		}
		return clampedPiOver2;
	}

	/// <summary>
	/// Returns the cosine of x.
	/// See Sin() for more details.
	/// </summary>
	public static Fix32 Cos(this Fix32 x) {
#if USE_DOUBLES
		return Math.Cos(x.ToDouble()).ToFix32();
#endif
		var xl = (int) x;
		var rawAngle = xl + (xl > 0 ? -PI - PI_OVER_2 : PI_OVER_2);
		return ((Fix32) rawAngle).Sin();
	}

	/// <summary>
	/// Returns a rough approximation of the cosine of x.
	/// See FastSin for more details.
	/// </summary>
	public static Fix32 CosFast(this Fix32 x) {
		var xl = (int) x;
		var rawAngle = xl + (xl > 0 ? -PI - PI_OVER_2 : PI_OVER_2);
		return ((Fix32) rawAngle).SinFast();
	}

	/// <summary>
	/// Returns the tangent of x.
	/// </summary>
	/// <remarks>
	/// This function is not well-tested. It may be wildly inaccurate.
	/// </remarks>
	public static Fix32 Tan(this Fix32 x) {
#if USE_DOUBLES
		return Math.Tan(x.ToDouble()).ToFix32();
#endif
		var clampedPi = (int) x % PI;
		var flip = false;
		if (clampedPi < 0) {
			clampedPi = -clampedPi;
			flip = true;
		}
		if (clampedPi > PI_OVER_2) {
			flip = !flip;
			clampedPi = PI_OVER_2 - (clampedPi - PI_OVER_2);
		}

		var clamped = (Fix32) clampedPi;

		// Find the two closest values in the LUT and perform linear interpolation
		Fix32 rawIndex = clamped.Mul(LutInterval);
		Fix32 roundedIndex = rawIndex.Round();
		Fix32 indexError = rawIndex.Sub(roundedIndex);

		Fix32 nearestValue = (Fix32) ((int) TanLut[(int) roundedIndex]);
		Fix32 secondNearestValue = (Fix32) ((int) TanLut[(int) roundedIndex + SignI(indexError)]);

		var delta = (int) indexError.Mul(nearestValue.Sub(secondNearestValue).Abs());
		var interpolatedValue = (int) nearestValue + delta;
		var finalValue = flip ? -interpolatedValue : interpolatedValue;
		return (Fix32) finalValue;
	}


	/// <summary>
	/// Returns the arctan of of the specified number, calculated using Euler series
	/// </summary>
	public static Fix32 Atan(this Fix32 zz) {
#if USE_DOUBLES
		return Math.Atan(z.ToDouble()).ToFix32();
#endif
		Fix32 z = zz;
		if ((int) z == 0)
			return Fix32.Zero;

		// Force positive values for argument
		// Atan(-z) = -Atan(z).
		bool neg = ((int) z < 0);
		if (neg) z = z.Neg();

		Fix32 result;

		if ((int) z == (int) Fix32.One)
			result = Fix32.PiOver4;
		else {
			bool invert = (int) z > (int) Fix32.One;
			if (invert) z = Fix32.One.Div(z);

			result = Fix32.One;
			Fix32 term = Fix32.One;

			Fix32 zSq = z.Mul(z);
			Fix32 zSq2 = zSq.Mul(Fix32.Two);
			Fix32 zSqPlusOne = zSq.Add(Fix32.One);
			Fix32 zSq12 = zSqPlusOne.Mul(Fix32.Two);
			Fix32 dividend = zSq2;
			Fix32 divisor = zSqPlusOne.Mul(Fix32.Three);

			for (int i = 2; i < 30; i++) {
				term = term.Mul(dividend.Div(divisor));
				result = result.Add(term);

				dividend = dividend.Add(zSq2);
				divisor = divisor.Add(zSq12);

				if ((int) term == 0)
					break;
			}

			result = result.Mul(z).Div(zSqPlusOne);

			if (invert)
				result = Fix32.PiOver2.Sub(result);
		}

		if (neg) result = result.Neg();
		return result;
	}

	public static Fix32 Atan2(this Fix32 y, Fix32 x) {
#if USE_DOUBLES
		return Math.Atan2(y.ToDouble(), x.ToDouble()).ToFix32();
#endif
		var yl = (int) y;
		var xl = (int) x;
		if (xl == 0) {
			if (yl > 0)
				return Fix32.PiOver2;
			if (yl == 0)
				return Fix32.Zero;
			return Fix32.PiOver2.Neg();
		}
		Fix32 atan;
		var z = y.Div(x);

		// Deal with overflow
		if ((int) Fix32.One.Add(C0p28.Mul(z).Mul(z)) == (int) Fix32.MaxValue) {
			return (int) y < 0 ? Fix32.PiOver2.Neg() : Fix32.PiOver2;
		}

		if ((int) Abs(z) < (int) Fix32.One) {
			atan = z.Div(Fix32.One.Add(C0p28.Mul(z).Mul(z)));
			if (xl < 0) {
				if (yl < 0) {
					return atan.Sub(Fix32.Pi);
				}
				return atan.Add(Fix32.Pi);
			}
		}
		else {
			atan = Fix32.PiOver2.Sub(z.Div(z.Mul(z).Add(C0p28)));
			if (yl < 0)
				return atan.Sub(Fix32.Pi);
		}

		return atan;
	}

	public static Fix32 Atan2Fast(this Fix32 y, Fix32 x) {
#if USE_DOUBLES
		return Math.Atan2(y.ToDouble(), x.ToDouble()).ToFix32();
#endif
		var yl = (int) y;
		var xl = (int) x;
		if (xl == 0) {
			if (yl > 0) {
				return Fix32.PiOver2;
			}
			if (yl == 0) {
				return Fix32.Zero;
			}
			return Fix32.PiOver2.Neg();
		}
		Fix32 atan;
		var z = y.Div(x);

		// Deal with overflow
		if ((int) Fix32.One.Add(C0p28.Mul(z).Mul(z)) == (int) Fix32.MaxValue) {
			return (int) y < 0 ? Fix32.PiOver2.Neg() : Fix32.PiOver2;
		}

		if ((int) Abs(z) < (int) Fix32.One) {
			atan = z.Div(Fix32.One.Add(C0p28.Mul(z).Mul(z)));
			if (xl < 0) {
				if (yl < 0) {
					return atan.Sub(Fix32.Pi);
				}
				return atan.Add(Fix32.Pi);
			}
		}
		else {
			atan = Fix32.PiOver2.Sub(z.Div(z.Mul(z).Add(C0p28)));
			if (yl < 0) {
				return atan.Sub(Fix32.Pi);
			}
		}
		return atan;
	}


	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	public static Fix32 ToFix32(this float value) {
		return (Fix32) (int) Clamp(value * ONE, MIN_VALUE, MAX_VALUE);
	}
	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	public static float ToFloat(this Fix32 value) {
		return (float) value / ONE;
	}

	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	public static Fix32 ToFix32(this double value) {
		return (Fix32) (int) Clamp(value * ONE, MIN_VALUE, MAX_VALUE);
	}
	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	public static double ToDouble(this Fix32 value) {
		return (double) value / ONE;
	}

	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	public static Fix32 ToFix32(this decimal value) {
		return (Fix32) (int) Clamp(value * ONE, MIN_VALUE, MAX_VALUE);
	}
	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	public static decimal ToDecimal(this Fix32 value) {
		return (decimal) value / ONE;
	}

	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	public static Fix32 ToFix32(this int value) {
		return (Fix32) (value > MAX_INT_VALUE ? MAX_VALUE : value < MIN_INT_VALUE ? MIN_VALUE : value * ONE);
	}
	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	public static int ToInt(this Fix32 value) {
		return (int) value / ONE;
	}


	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	public static double Clamp(float value, double min, double max) {
		return value > max ? max : value < min ? min : value;
	}
	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	public static double Clamp(double value, double min, double max) {
		return value > max ? max : value < min ? min : value;
	}
	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	public static decimal Clamp(decimal value, decimal min, decimal max) {
		return value > max ? max : value < min ? min : value;
	}


	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	public static int CompareTo(this Fix32 x, Fix32 other) {
		return ((int) x).CompareTo((int) other);
	}

	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	public static string ToStringExt(this Fix32 x) {
		return x.ToDouble().ToString("0.##########");
	}

	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	public static Fix32 Saturate(long value) {
		return (Fix32) (int) (((((ulong) value) >> 63) - 1U) ^ (1U << NUM_BITS_MINUS_ONE));
	}

#if UNITY_EDITOR
	static void GenerateSinLut(string where) {
		CalculateLargePI(out Fix32 largePIF32, out int N);

		using (var writer = new System.IO.StreamWriter(where + "/Fix32SinLut.cs")) {
            writer.Write(
@"namespace FixMath.NET {
public static partial class Fix32Ext {
	const int LARGE_PI_RAW = " + (int) largePIF32 + @";
	const int LARGE_PI_TIMES = " + N + @";

	internal static readonly int[] SinLut = new int[" + LUT_SIZE + "] {");
            int lineCounter = 0;
            for (int i = 0; i < LUT_SIZE; ++i) {
                var angle = i * 3.1415926535897932384626433832795028841971693993751058209749445923078164m * 0.5m / (LUT_SIZE - 1);
                if (lineCounter++ % 8 == 0) {
                    writer.WriteLine();
                    writer.Write("			");
                }
                var sin = Math.Sin((double) angle);
                var rawValue = (int) sin.ToFix32();
                writer.Write(string.Format("0x{0:X}U, ", rawValue));
            }
            writer.Write(
@"
	};
}
}");
        }
    }

    static void GenerateTanLut(string where) {
        using (var writer = new System.IO.StreamWriter(where + "/Fix32TanLut.cs")) {
            writer.Write(
@"namespace FixMath.NET {
public static partial class Fix32Ext {
	internal static readonly int[] TanLut = new int[" + LUT_SIZE + "] {");
            int lineCounter = 0;
            for (int i = 0; i < LUT_SIZE; ++i) {
                var angle = i * 3.1415926535897932384626433832795028841971693993751058209749445923078164m * 0.5m / (LUT_SIZE - 1);
                if (lineCounter++ % 8 == 0) {
                    writer.WriteLine();
                    writer.Write("			");
                }
                var tan = Math.Tan((double) angle);
                if (tan > Fix32.MaxValue.ToDouble() || tan < 0.0) {
                    tan = Fix32.MaxValue.ToDouble();
                }
                var rawValue = (int) (((decimal)tan > Fix32.MaxValue.ToDecimal() || tan < 0.0) ? Fix32.MaxValue : tan.ToFix32());
                writer.Write(string.Format("0x{0:X}U, ", rawValue));
            }
            writer.Write(
@"
	};
}
}");
        }
    }

	// Obtained from ((Fix64)1686629713.065252369824872831112M).RawValue
	// This is (2^29)*PI, where 29 is the largest N such that (2^N)*PI < MaxValue.
	// The idea is that this number contains way more precision than PI_TIMES_2,
	// and (((x % (2^29*PI)) % (2^28*PI)) % ... (2^1*PI) = x % (2 * PI)
	// In practice this gives us an error of about 1,25e-9 in the worst case scenario (Sin(MaxValue))
	// Whereas simply doing x % PI_TIMES_2 is the 2e-3 range.
	static void CalculateLargePI(out Fix32 largePIF64, out int N) {
		int prevN = 0;
		N = 0;
		decimal largePI = 0;
		decimal prevLargePI = 0;

		do {
			prevN = N++;
			prevLargePI = largePI;

			decimal sqrd = 2;
			for (int i = 1; i < N; i++) {
				sqrd *= 2;
			}
			largePI = sqrd * (decimal) Math.PI;
		}
		while (largePI < Fix32.MaxValue.ToDecimal());
		N = prevN;
		largePIF64 = prevLargePI.ToFix32();
	}

	// turn into a Console Application and use this to generate the look-up tables
	[UnityEditor.MenuItem("Tools/Fix64 Regenerate LUT")]
    static void GenerateLut() {
		string thisFile = new System.Diagnostics.StackTrace(true).GetFrame(0).GetFileName();
		thisFile = thisFile.Replace('\\', '/');
		string path = thisFile.Substring(0, thisFile.LastIndexOf('/'));
        GenerateSinLut(path);
        GenerateTanLut(path);
		UnityEditor.AssetDatabase.Refresh();
    }
#endif
}
