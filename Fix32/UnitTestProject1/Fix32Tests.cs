using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

public class Fix32Tests {
	static int OneRaw = (int) Fix32.One;

	static List<int> TestCases = new List<int>() {
        // Small numbers
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
		-1, -2, -3, -4, -5, -6, -7, -8, -9, -10,
  
        // Integer numbers
		1 * OneRaw, -1 * OneRaw, 2 * OneRaw, -2 * OneRaw, 3 * OneRaw, -3 * OneRaw,
		4 * OneRaw, -4 * OneRaw, 5 * OneRaw, -5 * OneRaw, 6 * OneRaw, -6 * OneRaw,
		7 * OneRaw, -7 * OneRaw, 8 * OneRaw, -8 * OneRaw, 9 * OneRaw, -9 * OneRaw,
		10 * OneRaw, -10 * OneRaw, 11 * OneRaw, -11 * OneRaw, 12 * OneRaw, -12 * OneRaw,
		13 * OneRaw, -13 * OneRaw, 14 * OneRaw, -14 * OneRaw, 15 * OneRaw, -15 * OneRaw,
		16 * OneRaw, -16 * OneRaw, 17 * OneRaw, -17 * OneRaw, 18 * OneRaw, -18 * OneRaw,
  
        // Fractions (1/2, 1/4, 1/8)
		OneRaw / 2, -OneRaw / 2, OneRaw / 4, -OneRaw / 4, OneRaw / 8, -OneRaw / 8,
		OneRaw / 16, -OneRaw / 16, OneRaw / 32, -OneRaw / 32, OneRaw / 64, -OneRaw / 64,
  
        // Problematic carry
		OneRaw - 1, - (OneRaw - 1),
		(OneRaw - 1) * 2 + 1, - ((OneRaw - 1) * 2 + 1),
		(OneRaw - 1) * 4 + 3, - ((OneRaw - 1) * 4 + 3),
  
        // Smallest and largest values
        int.MaxValue, int.MinValue,
			
        // Smallest and largest values minus one
        int.MaxValue - 1, int.MinValue + 1,
			
        // Smallest and largest values minus two
        int.MaxValue - 2, int.MinValue + 2,
  
        // Large random numbers
        630801836, -819412065, 622276300, -781232034,
		848005112, -767148444, 770879513, -521836768,
		223221456, -520369580, 686049339, -906349457,
		658689347, -655205245, 670526493,

        // Small random numbers
        -4368, -2246, 3294, 25981, 3398, 13, 1060,
		-34527, 1753, 28913, 24813, 2850, 265, -1447,
		30773,

        // Tiny random numbers
        -171,
		-359, 491, 844, 158, -413, -422, -737, -575, -330,
		-376, 435, -311, 116, 715, -1024, -487, 59, 724, 993
	};

	static Fix32Tests() {
		Random r = new Random(0);
		for (int i = 0; i < byte.MaxValue; i++)
			TestCases.Add(i);
		for (int i = 0; i < byte.MaxValue; i++)
			TestCases.Add(-i);
		for (int i = 0; i < 250; i++)
			TestCases.Add(r.Next() % short.MaxValue);
		for (int i = 0; i < 250; i++)
			TestCases.Add(r.Next());

		const bool BENCHMARK = true;
		if (BENCHMARK) {
			for (int i = 0; i < 1000; i++)
				TestCases.Add(r.Next());
		}
	}

	[Test]
	public void T001_IntToFix64AndBack() {
		List<int> sources = new List<int>() {
			int.MinValue,
			-1000000,
			-100000,
			-10000,
			-1000,
			-100,
			-10,
			-1,
			0,
			1,
			10,
			100,
			1000,
			10000,
			100000,
			1000000,
			int.MaxValue,
		};

		Random r = new Random(0);

		for (int i = 0; i < 100; i++) {
			sources.Add(r.Next());
		}

		foreach (var value in sources) {
			int expected = value > (int) Fix32.MaxValue.ToInt() ? (int) Fix32.MaxValue.ToInt() :
				value < (int) Fix32.MinValue.ToInt() ? (int) Fix32.MinValue.ToInt() :
				value;
			Assert.AreEqual(expected, value.ToFix32().ToInt(), value.ToString());
		}
	}

	[Test]
	public void T002_DoubleToFix64AndBack() {
		List<double> sources = new List<double>() {
			-int.MaxValue * 100d,
			-int.MaxValue * 2d,
			-int.MaxValue,
			int.MinValue,
			-(double)Math.PI,
			-(double)Math.E,
			-1.0,
			-0.0,
			0.0,
			1.0,
			(double)Math.PI,
			(double)Math.E,
			int.MaxValue,
			int.MaxValue * 2d,
			int.MaxValue * 100d,
		};

		Random r = new Random(0);

		for (int i = 0; i < 100; i++) {
			sources.Add(r.NextDouble());
		}

		foreach (var value in sources) {
			var expected = value > Fix32.MaxValue.ToDouble() ? Fix32.MaxValue.ToDouble() :
				value < Fix32.MinValue.ToDouble() ? Fix32.MinValue.ToDouble() :
				value;
			Assert.AreEqual(expected, value.ToFix32().ToDouble(), (double) Fix32Ext.Precision);
		}
	}

	[Test]
	public void T004_Add() {
		PrepareStatistics(out var deltas, out var swF, out var swD);

		var terms1 = new Fix32[] { Fix32.MinValue, Fix32.MinusOne, Fix32.Zero, Fix32.One, Fix32.MaxValue };
		var terms2 = new Fix32[] { Fix32.MinusOne, Fix32.Two, (-1.5).ToFix32(), Fix32.MinusTwo, Fix32.One };
		var expecteds = new Fix32[] { Fix32.MinValue, Fix32.One, (-1.5).ToFix32(), Fix32.MinusOne, Fix32.MaxValue };
		for (int i = 0; i < terms1.Length; ++i) {
			var actual = terms1[i].Add(terms2[i]);
			var expected = expecteds[i];
			Assert.AreEqual(expected, actual, terms1[i].ToStringExt() + " + " + terms2[i].ToStringExt());
		}


		for (int i = 0; i < TestCases.Count; ++i) {
			for (int j = 0; j < TestCases.Count; ++j) {
				var x = (Fix32) TestCases[i];
				var y = (Fix32) TestCases[j];

				swD.Start();
				double expected = Saturate(x.ToDouble() + y.ToDouble());
				swD.Stop();
				swF.Start();
				double actual = x.Add(y).ToDouble();
				swF.Stop();

				Assert.AreEqual(expected, actual, x.ToStringExt() + " + " + y.ToStringExt());
			}
		}
		PrintStatistics(deltas, swF, swD);
	}

	[Test]
	public void T005_Sub() {
		PrepareStatistics(out var deltas, out var swF, out var swD);

		var terms1 = new Fix32[] { Fix32.MinValue, Fix32.MinusOne, Fix32.Zero, Fix32.One, Fix32.MaxValue };
		var terms2 = new Fix32[] { Fix32.One, Fix32.MinusTwo, (1.5).ToFix32(), Fix32.Two, Fix32.MinusOne };
		var expecteds = new Fix32[] { Fix32.MinValue, Fix32.One, (-1.5).ToFix32(), Fix32.MinusOne, Fix32.MaxValue };
		for (int i = 0; i < terms1.Length; ++i) {
			var actual = terms1[i].Sub(terms2[i]);
			var expected = expecteds[i];
			Assert.AreEqual(expected, actual, terms1[i].ToStringExt() + " - " + terms2[i].ToStringExt());
		}

		for (int i = 0; i < TestCases.Count; ++i) {
			for (int j = 0; j < TestCases.Count; ++j) {
				var x = (Fix32) TestCases[i];
				var y = (Fix32) TestCases[j];

				swD.Start();
				double expected = Saturate(x.ToDouble() - y.ToDouble());
				swD.Stop();
				swF.Start();
				double actual = x.Sub(y).ToDouble();
				swF.Stop();

				Assert.AreEqual(expected, actual, x.ToStringExt() + " - " + y.ToStringExt());
			}
		}
		PrintStatistics(deltas, swF, swD);
	}

	[Test]
	public void T006_Neg() {
		foreach (var operand1 in TestCases) {
			var f = (Fix32) operand1;
			if ((int) f == (int) Fix32.MinValue) {
				Assert.AreEqual(f.Neg(), Fix32.MaxValue);
			}
			else {
				var expected = -((decimal) f);
				var actual = (decimal) (f.Neg());
				Assert.AreEqual(expected, actual);
			}
		}
	}

	[Test]
	public void T007_EqualityInequalityComparisonOperators() {
		List<Fix32> sources = TestCases.Select(t => (Fix32) t).ToList();
		foreach (var op1 in sources) {
			foreach (var op2 in sources) {
				var d1 = op1.ToDouble();
				var d2 = op2.ToDouble();
				Assert.AreEqual(d1 == d2, (int) op1 == (int) op2);
				Assert.AreEqual(d1.Equals(d2), op1.Equals(op2));
				Assert.AreEqual(d1 != d2, (int) op1 != (int) op2);
				Assert.AreNotEqual((int) op1 != (int) op2, (int) op1 == (int) op2);
				Assert.AreEqual(d1 < d2, (int) op1 < (int) op2);
				Assert.AreEqual(d1 <= d2, (int) op1 <= (int) op2);
				Assert.AreEqual(d1 > d2, (int) op1 > (int) op2);
				Assert.AreEqual(d1 >= d2, (int) op1 >= (int) op2);
			}
		}
	}

	[Test]
	public void T008_CompareTo() {
		var nums = TestCases.Select(t => (Fix32) t).ToArray();
		var numsDecimal = nums.Select(t => t.ToDouble()).ToArray();
		Array.Sort(nums);
		Array.Sort(numsDecimal);
		CollectionAssert.AreEqual(numsDecimal, nums.Select(t => t.ToDouble()));
	}

	[Test]
	public void T009_Sign() {
		var sources = new Fix32[] { Fix32.MinValue, (-128).ToFix32(), (-100).ToFix32(), Fix32.MinusOne, Fix32.Zero, Fix32.One, 100.ToFix32(), 128.ToFix32(), Fix32.MaxValue };
		var expecteds = new int[] { -1, -1, -1, -1, 0, 1, 1, 1, 1 };
		for (int i = 0; i < sources.Length; ++i) {
			int expected = expecteds[i];

			var actual = sources[i].Sign();
			Assert.AreEqual((double) expected, actual.ToDouble(), sources[i].ToStringExt());
			Assert.AreEqual(expected, actual.ToInt(), sources[i].ToStringExt());

			int actual2 = sources[i].SignI();
			Assert.AreEqual(expected, actual2, sources[i].ToStringExt());
		}
	}

	[Test]
	public void T010_Abs() {
		Assert.AreEqual(Fix32.MaxValue, Fix32.MinValue.Abs());
		var sources = new Fix32[] { Fix32.MaxValue.Neg(), Fix32.MinusOne, Fix32.Zero, Fix32.One, Fix32.MaxValue.Sub(Fix32.One), Fix32.MaxValue };
		var expecteds = new Fix32[] { Fix32.MaxValue, Fix32.One, Fix32.Zero, Fix32.One, Fix32.MaxValue.Sub(Fix32.One), Fix32.MaxValue };
		for (int i = 0; i < sources.Length; ++i) {
			var actual = sources[i].Abs();
			var expected = expecteds[i];
			Assert.AreEqual(expected, actual, sources[i].ToStringExt());
		}

		for (int i = 0; i < TestCases.Count; ++i) {
			var actual = ((Fix32) TestCases[i]).Abs();
			var expected = Math.Abs(((Fix32) TestCases[i]).ToDouble()).ToFix32();

			Assert.AreEqual(expected, actual, TestCases[i].ToFix32().ToStringExt());
		}
	}

	[Test]
	public void T011_AbsFast() {
		Assert.AreEqual(Fix32.MinValue, Fix32.MinValue.AbsFast()); // Wrong result, but it is spected

		var sources = new Fix32[] { Fix32.MaxValue.Neg(), Fix32.MinusOne, Fix32.Zero, Fix32.One, Fix32.MaxValue.Sub(Fix32.One), Fix32.MaxValue };
		var expecteds = new Fix32[] { Fix32.MaxValue, Fix32.One, Fix32.Zero, Fix32.One, Fix32.MaxValue.Sub(Fix32.One), Fix32.MaxValue };
		for (int i = 0; i < sources.Length; ++i) {
			var actual = sources[i].AbsFast();
			var expected = expecteds[i];
			Assert.AreEqual(expected, actual, sources[i].ToStringExt());
		}

		for (int i = 0; i < TestCases.Count; ++i) {
			var actual = ((Fix32) TestCases[i]).AbsFast();
			var expected = (Fix32) TestCases[i];
			if ((int) expected != (int) Fix32.MinValue)
				expected = Math.Abs(((Fix32) TestCases[i]).ToDouble()).ToFix32();

			Assert.AreEqual(expected, actual, TestCases[i].ToFix32().ToStringExt());
		}
	}

	[Test]
	public void T012_Floor() {
		var sources = new[] { -5.1, -1, 0, 1, 5.1, Fix32.MaxValue.ToDouble(), Fix32.MinValue.ToDouble() };
		var expecteds = new[] { -6, -1, 0, 1, 5, Fix32.MaxValue.ToInt(), Fix32.MinValue.ToDouble() };
		for (int i = 0; i < sources.Length; ++i) {
			var actual = sources[i].ToFix32().Floor().ToDouble();
			var expected = expecteds[i];
			Assert.AreEqual(expected, actual, sources[i].ToString());
		}
	}

	[Test]
	public void T013_Ceiling() {
		var sources = new[] { -5.1, -1, 0, 1, 5.1, Fix32.MaxValue.ToDouble(), Fix32.MinValue.ToDouble() };
		var expecteds = new[] { -5, -1, 0, 1, 6, Fix32.MaxValue.ToDouble(), Fix32.MinValue.ToDouble() };
		for (int i = 0; i < sources.Length; ++i) {
			var actual = sources[i].ToFix32().Ceiling().ToDouble();
			var expected = expecteds[i];
			Assert.AreEqual(expected, actual, sources[i].ToString());
		}
	}

	[Test]
	public void T014_Round() {
		var sources = new[] { 5.5, 5.4, 4.6, 4.5, 1, 0, -1, -4.4, -4.5, -5.1, -5.5 };
		var expecteds = new[] { 6, 5, 5, 4, 1, 0, -1, -4, -4, -5, -6 };
		for (int i = 0; i < sources.Length; ++i) {
			var actual = sources[i].ToFix32().Round().ToDouble();
			var expected = expecteds[i];
			Assert.AreEqual(expected, actual, sources[i].ToString());
		}
		Assert.AreEqual(Fix32.MaxValue, Fix32.MaxValue.Round());
	}

	[Test]
	public void T014_RoundFast() {
		var sources = new[] { 5.5, 5.4, 4.6, 4.5, 1, 0, -1, -4.4, -4.5, -5.1, -5.5 };
		var expecteds = new[] { 6, 5, 5, 4, 1, 0, -1, -4, -4, -5, -6 };
		for (int i = 0; i < sources.Length; ++i) {
			var actual = sources[i].ToFix32().RoundFast().ToDouble();
			var expected = expecteds[i];
			Assert.AreEqual(expected, actual, sources[i].ToString());
		}
	}

	[Test]
	public void T015_BasicMult() {
		var term1s = new[] { 0, 1, -1, 5, -5, 0.5, -0.5, -1.0 };
		var term2s = new[] { 16, 16, 16, 16, 16, 16, 16, -1.0 };
		var expecteds = new double[] { 0, 16, -16, 80, -80, 8, -8, 1 };
		for (int i = 0; i < term1s.Length; ++i) {
			var expected = expecteds[i];
			var actual = term1s[i].ToFix32().Mul(term2s[i].ToFix32()).ToDouble();
			Assert.AreEqual(expected, actual, term1s[i] + " * " + term2s[i]);
		}
	}

	[Test]
	public void T015_BasicFastMult() {
		var term1s = new[] { 0, 1, -1, 5, -5, 0.5, -0.5, -1.0 };
		var term2s = new[] { 16, 16, 16, 16, 16, 16, 16, -1.0 };
		var expecteds = new double[] { 0, 16, -16, 80, -80, 8, -8, 1 };
		for (int i = 0; i < term1s.Length; ++i) {
			var expected = expecteds[i];
			var actual = term1s[i].ToFix32().MulFast(term2s[i].ToFix32()).ToDouble();
			Assert.AreEqual(expected, actual, term1s[i] + " * " + term2s[i]);
		}
	}

	[Test]
	public void T016_Mult() {
		PrepareStatistics(out var deltas, out var swF, out var swD);

		for (int i = 0; i < TestCases.Count; ++i) {
			for (int j = 0; j < TestCases.Count; ++j) {
				var x = (Fix32) TestCases[i];
				var y = (Fix32) TestCases[j];
				var xM = x.ToDouble();
				var yM = y.ToDouble();
				swD.Start();
				var expectedF = xM * yM;
				swD.Stop();
				expectedF =
					expectedF > Fix32.MaxValue.ToDouble() ? Fix32.MaxValue.ToDouble() :
					expectedF < Fix32.MinValue.ToDouble() ? Fix32.MinValue.ToDouble() :
					expectedF;
				swF.Start();
				var actualF = x.Mul(y);
				swF.Stop();
				var expected = expectedF.ToFix32().ToDouble();
				var actual = actualF.ToDouble();
				Assert.AreEqual(expected, actual, (double) Fix32Ext.Precision, x.ToStringExt() + " * " + y.ToStringExt() + Delta(expected, actual, deltas));
			}
		}
		PrintStatistics(deltas, swF, swD);
	}

	[Test]
	public void T016_FastMult() {
		PrepareStatistics(out var deltas, out var swF, out var swD);

		for (int i = 0; i < TestCases.Count; ++i) {
			for (int j = 0; j < TestCases.Count; ++j) {
				var x = (Fix32) TestCases[i];
				var y = (Fix32) TestCases[j];
				var xM = x.ToDouble();
				var yM = y.ToDouble();
				swD.Start();
				var expectedF = xM * yM;
				swD.Stop();
				if (expectedF > Fix32.MaxValue.ToDouble() || expectedF < Fix32.MinValue.ToDouble())
					continue; // Fast mult version doesn't saturate
				swF.Start();
				var actualF = x.MulFast(y);
				swF.Stop();
				var expected = expectedF.ToFix32().ToDouble();
				var actual = actualF.ToDouble();
				Assert.AreEqual(expected, actual, (double) Fix32Ext.Precision, x.ToStringExt() + " * " + y.ToStringExt() + Delta(expected, actual, deltas));
			}
		}
		PrintStatistics(deltas, swF, swD);
	}

	[Test]
	public void T017_Div() {
		PrepareStatistics(out var deltas, out var swF, out var swD);

		for (int i = 0; i < TestCases.Count; ++i) {
			for (int j = 0; j < TestCases.Count; ++j) {
				var x = (Fix32) TestCases[i];
				var y = (Fix32) TestCases[j];
				var xM = x.ToDouble();
				var yM = y.ToDouble();

				if (TestCases[j] == 0) {
					Assert.AreEqual(x >= 0 ? Fix32.MaxValue : Fix32.MinValue, x.Div(y), x.ToStringExt() + " / " + y.ToStringExt());
				}
				else {
					swD.Start();
					var expectedF = xM / yM;
					swD.Stop();
					expectedF =
						expectedF > Fix32.MaxValue.ToDouble()
							? Fix32.MaxValue.ToDouble()
							: expectedF < Fix32.MinValue.ToDouble()
									? Fix32.MinValue.ToDouble()
									: expectedF;
					swF.Start();
					var actualF = x.Div(y);
					swF.Stop();
					var expected = (double) expectedF;
					var actual = actualF.ToDouble();

					Assert.AreEqual(expected, actual, (double) Fix32Ext.Precision, x.ToStringExt() + " / " + y.ToStringExt() + Delta(expected, actual, deltas));
				}
			}
		}
		PrintStatistics(deltas, swF, swD);
	}

	[Test]
	public void T017_DivFast() {
		PrepareStatistics(out var deltas, out var swF, out var swD);

		for (int i = 0; i < TestCases.Count; ++i) {
			for (int j = 0; j < TestCases.Count; ++j) {
				var x = (Fix32) TestCases[i];
				var y = (Fix32) TestCases[j];
				var xM = x.ToDouble();
				var yM = y.ToDouble();

				if (TestCases[j] == 0) {
					Assert.AreEqual(x >= 0 ? Fix32.MaxValue : Fix32.MinValue, x.DivFast(y), x.ToStringExt() + " / " + y.ToStringExt());
				}
				else {
					swD.Start();
					var expectedF = xM / yM;
					swD.Stop();
					swF.Start();
					var actualF = x.DivFast(y);
					swF.Stop();

					if (expectedF > Fix32.MaxValue.ToDouble() || expectedF < Fix32.MinValue.ToDouble()) continue;

					var expected = (double) expectedF;
					var actual = actualF.ToDouble();

					Assert.AreEqual(expected, actual, (double) Fix32Ext.Precision, x.ToStringExt() + " / " + y.ToStringExt() + Delta(expected, actual, deltas));
				}
			}
		}
		PrintStatistics(deltas, swF, swD);
	}

	[Test]
	public void T018_Log2() {
		PrepareStatistics(out var deltas, out var swF, out var swD);
		double maxDelta = 1 * (double) Fix32Ext.Precision;

		for (int j = 0; j < TestCases.Count; ++j) {
			var b = (Fix32) TestCases[j];

			if (b <= Fix32.Zero) {
				Assert.Throws<ArgumentOutOfRangeException>(() => b.Log2());
			}
			else {
				swD.Start();
				//var expected = Math.Log(b.ToDouble()) / Math.Log(2);
				var expected = Math.Log(b.ToDouble(), 2);
				swD.Stop();
				swF.Start();
				var actual = b.Log2().ToDouble();
				swF.Stop();

				Assert.AreEqual(expected, actual, maxDelta, b.ToStringExt() + Delta(expected, actual, deltas));
			}
		}
		PrintStatistics(deltas, swF, swD);
	}

	[Test]
	public void T018_Log2Fast() {
		PrepareStatistics(out var deltas, out var swF, out var swD);
		double maxDelta = 4 * (double) Fix32Ext.Precision;

		for (int j = 0; j < TestCases.Count; ++j) {
			var b = (Fix32) TestCases[j];

			if (b <= Fix32.Zero) {
				Assert.Throws<ArgumentOutOfRangeException>(() => b.Log2Fast());
			}
			else {
				swD.Start();
				//var expected = Math.Log(b.ToDouble()) / Math.Log(2);
				var expected = Math.Log(b.ToDouble(), 2);
				swD.Stop();
				swF.Start();
				var actual = b.Log2Fast().ToDouble();
				swF.Stop();

				Assert.AreEqual(expected, actual, maxDelta, b.ToStringExt() + Delta(expected, actual, deltas));
			}
		}
		PrintStatistics(deltas, swF, swD);
	}

	[Test]
	public void T019_Ln() {
		PrepareStatistics(out var deltas, out var swF, out var swD);
		double maxDelta = 8 * (double) Fix32Ext.Precision;

		for (int j = 0; j < TestCases.Count; ++j) {
			var b = (Fix32) TestCases[j];

			if (b <= Fix32.Zero) {
				Assert.Throws<ArgumentOutOfRangeException>(() => b.Ln());
			}
			else {
				swD.Start();
				var expected = Math.Log(b.ToDouble());
				swD.Stop();
				swF.Start();
				var actual = b.Ln().ToDouble();
				swF.Stop();

				Assert.AreEqual(expected, actual, maxDelta, b.ToStringExt() + Delta(expected, actual, deltas));
			}
		}
		PrintStatistics(deltas, swF, swD);
	}

	[Test]
	public void T020_SqrtSlow() {
		PrepareStatistics(out var deltas, out var swF, out var swD);
		double maxDelta = 1 * (double) Fix32Ext.Precision;

		for (int i = 0; i < TestCases.Count; ++i) {
			var f = (Fix32) TestCases[i];
			if (f.SignI() < 0) {
				Assert.AreEqual(0, (int) f.SqrtSlow(), "sqrt(" + f.ToStringExt() + ")");
			}
			else {
				swD.Start();
				var expected = Math.Sqrt(f.ToDouble());
				swD.Stop();
				swF.Start();
				var actual = f.SqrtSlow().ToDouble();
				swF.Stop();
				Assert.AreEqual(expected, actual, maxDelta, "sqrt(" + f.ToStringExt() + ")" + Delta(expected, actual, deltas));
			}
		}
		PrintStatistics(deltas, swF, swD);
	}

	[Test]
	public void T020_Sqrt() {
		PrepareStatistics(out var deltas, out var swF, out var swD);
		double maxDelta = 2 * (double) Fix32Ext.Precision;

		for (int i = 0; i < TestCases.Count; ++i) {
			var f = (Fix32) TestCases[i];
			if (f.SignI() < 0) {
				Assert.AreEqual(0, (int) f.Sqrt(), "sqrt(" + f.ToStringExt() + ")");
			}
			else {
				swD.Start();
				var expected = Math.Sqrt(f.ToDouble());
				swD.Stop();
				swF.Start();
				var actual = f.Sqrt().ToDouble();
				swF.Stop();
				Assert.AreEqual(expected, actual, maxDelta, "sqrt(" + f.ToStringExt() + ")" + Delta(expected, actual, deltas));
			}
		}
		PrintStatistics(deltas, swF, swD);
	}

	[Test]
	public void T021_Modulus() {
		foreach (var operand1 in TestCases) {
			foreach (var operand2 in TestCases) {
				var f1 = (Fix32) operand1;
				var f2 = (Fix32) operand2;

				if (operand2 == 0) {
					Assert.AreEqual(f1 >= 0 ? Fix32.MaxValue : Fix32.MinValue, f1.Div(f2));
				}
				else {
					var d1 = f1.ToDouble();
					var d2 = f2.ToDouble();
					var actual = f1.Mod(f2).ToDouble();
					var expected = d1 % d2;
					Assert.AreEqual(expected, actual, f1.ToStringExt() + " % " + f2.ToStringExt());
				}
			}
		}
	}

	[Test]
	public void T022_Pow2() {
		PrepareStatistics(out var deltas, out var swF, out var swD);

		for (int i = 0; i < TestCases.Count; ++i) {
			var e = (Fix32) TestCases[i];

			swD.Start();
			var expected = Math.Min(Math.Pow(2, e.ToDouble()), Fix32.MaxValue.ToDouble());
			swD.Stop();
			swF.Start();
			var actual = e.Pow2().ToDouble();
			swF.Stop();

			// Absolute precision deteriorates with large result values, take this into account
			double maxDelta =
				expected > 1000 ? 0.05 :
				expected > 100 ? 0.05 :
				expected > 10 ? 0.002 :
				expected > 5 ? 0.0001 :
				expected > 3 ? 32 * (double) Fix32Ext.Precision :
				expected > 2 ? 24 * (double) Fix32Ext.Precision :
				expected > 1 ? 8 * (double) Fix32Ext.Precision :
				4 * (double) Fix32Ext.Precision;

			Assert.AreEqual(expected, actual, maxDelta, "Pow2(" + e.ToStringExt() + ")" + Delta(expected, actual, deltas));
		}
		PrintStatistics(deltas, swF, swD);
	}

	[Test]
	public void T022_Pow2Fast() {
		PrepareStatistics(out var deltas, out var swF, out var swD);

		for (int i = 0; i < TestCases.Count; ++i) {
			var e = (Fix32) TestCases[i];

			swD.Start();
			var expected = Math.Min(Math.Pow(2, e.ToDouble()), Fix32.MaxValue.ToDouble());
			swD.Stop();
			swF.Start();
			var actual = e.Pow2().ToDouble();
			swF.Stop();

			// Absolute precision deteriorates with large result values, take this into account
			double maxDelta =
				expected > 1000 ? 0.05 :
				expected > 100 ? 0.05 :
				expected > 10 ? 0.002 :
				expected > 5 ? 0.0001 :
				expected > 5 ? 40 * (double) Fix32Ext.Precision :
				expected > 3 ? 32 * (double) Fix32Ext.Precision :
				expected > 2 ? 24 * (double) Fix32Ext.Precision :
				expected > 1 ? 8 * (double) Fix32Ext.Precision :
				4 * (double) Fix32Ext.Precision;

			Assert.AreEqual(expected, actual, maxDelta, "Pow2(" + e.ToStringExt() + ")" + Delta(expected, actual, deltas));
		}
		PrintStatistics(deltas, swF, swD);
	}

	[Test]
	public void T023_Pow() {
		PrepareStatistics(out var deltas, out var swF, out var swD);

		for (int i = 0; i < TestCases.Count; ++i) {
			var b = (Fix32) TestCases[i];

			for (int j = 0; j < TestCases.Count; ++j) {
				var e = (Fix32) TestCases[j];

				if (b == Fix32.Zero && e < Fix32.Zero) {
					Assert.Throws<DivideByZeroException>(() => b.Pow(e));
				}
				else if (b < Fix32.Zero && e != Fix32.Zero) {
					Assert.Throws<ArgumentOutOfRangeException>(() => b.Pow(e));
				}
				else {
					swD.Start();
					var expected = e == Fix32.Zero ? 1 : b == Fix32.Zero ? 0 : Math.Min(Math.Pow(b.ToDouble(), e.ToDouble()), Fix32.MaxValue.ToDouble());
					swD.Stop();

					// Absolute precision deteriorates with large result values, take this into account
					// Similarly, large exponents reduce precision, even if result is small.
					double maxDelta = 100000000000000 * (double) Fix32Ext.Precision;
					/*double maxDelta = 
					expected > 100 ? 5 :
					expected > 10 ? 1 :
					expected > 5 ? 0.1 :
					expected > 2 ? 0.01 :
					expected > 1 ? 0.001 :
					0.0001;*/

					swF.Start();
					var actual = b.Pow(e).ToDouble();
					swF.Stop();

					Assert.AreEqual(expected, actual, maxDelta, "Pow(" + b.ToStringExt() + ", " + e.ToStringExt() + ")" + Delta(expected, actual, deltas));
				}
			}
		}
		PrintStatistics(deltas, swF, swD);
	}
	
	[Test]
	public void T024_Sin() {
		PrepareStatistics(out var deltas, out var swF, out var swD);

		Assert.AreEqual(Fix32.Zero.ToDouble(), Fix32.Zero.Sin().ToDouble(), Fix32Ext.Precision, "sin(0)");

		Assert.AreEqual(Fix32.One.ToDouble(), Fix32.PiOver2.Sin().ToDouble(), Fix32Ext.Precision, "sin(PI^2)");
		Assert.AreEqual(Fix32.Zero.ToDouble(), Fix32.Pi.Sin().ToDouble(), Fix32Ext.Precision, "sin(PI)");
		Assert.AreEqual(Fix32.MinusOne.ToDouble(), (Fix32.Pi.Add(Fix32.PiOver2)).Sin().ToDouble(), Fix32Ext.Precision, "sin(PI + PI^2)");
		Assert.AreEqual(Fix32.Zero.ToDouble(), Fix32.PiTimes2.Sin().ToDouble(), Fix32Ext.Precision, "sin(2 * PI)");

		Assert.AreEqual(Fix32.MinusOne.ToDouble(), Fix32.PiOver2.Neg().Sin().ToDouble(), Fix32Ext.Precision, "sin(-PI^2)");
		Assert.AreEqual(Fix32.Zero.ToDouble(), Fix32.Pi.Neg().Sin().ToDouble(), Fix32Ext.Precision, "sin(-PI)");
		Assert.AreEqual(Fix32.One.ToDouble(), (Fix32.Pi.Neg().Sub(Fix32.PiOver2)).Sin().ToDouble(), Fix32Ext.Precision, "sin(-PI - PI^2)");
		Assert.AreEqual(Fix32.Zero.ToDouble(), Fix32.PiTimes2.Neg().Sin().ToDouble(), Fix32Ext.Precision, "sin(-2 * PI)"); // This doesn't return exactly 0


		for (double angle = -2 * Math.PI; angle <= 2 * Math.PI; angle += 0.0001) {
			var f = angle.ToFix32();
			swF.Start();
			var actualF = f.Sin();
			swF.Stop();
			swD.Start();
			var expected = Math.Sin(angle);
			swD.Stop();
			Assert.AreEqual(expected, actualF.ToDouble(), 4 * Fix32Ext.Precision, "Sin(" + angle + ")" + Delta(expected, actualF.ToDouble(), deltas));
		}
		
		foreach (var val in TestCases) {
			var f = (Fix32) val;
			swF.Start();
			var actualF = f.Sin();
			swF.Stop();
			swD.Start();
			var expected = Math.Sin(f.ToDouble());
			swD.Stop();
			Assert.AreEqual(expected, actualF.ToDouble(), 4 * Fix32Ext.Precision, "Sin(" + f.ToDouble() + ")" + Delta(expected, actualF.ToDouble(), deltas));
		}
		PrintStatistics(deltas, swF, swD);
	}

	[Test]
	public void T024_SinFast() {
		PrepareStatistics(out var deltas, out var swF, out var swD);

		for (double angle = -2 * Math.PI; angle <= 2 * Math.PI; angle += 0.0001) {
			var f = angle.ToFix32();
			swF.Start();
			var actualF = f.SinFast();
			swF.Stop();
			swD.Start();
			var expected = Math.Sin(angle);
			swD.Stop();
			Assert.AreEqual(expected, actualF.ToDouble(), 4 * Fix32Ext.Precision, "Sin(" + angle + ")" + Delta(expected, actualF.ToDouble(), deltas));
		}

		foreach (var val in TestCases) {
			var f = (Fix32) val;
			swF.Start();
			var actualF = f.SinFast();
			swF.Stop();
			swD.Start();
			var expected = Math.Sin(f.ToDouble());
			swD.Stop();
			Assert.AreEqual(expected, actualF.ToDouble(), 4 * Fix32Ext.Precision, "Sin(" + f.ToDouble() + ")" + Delta(expected, actualF.ToDouble(), deltas));
		}
		PrintStatistics(deltas, swF, swD);
	}

	[Test]
	public void T024_SinFastest() {
		PrepareStatistics(out var deltas, out var swF, out var swD);

		for (double angle = -2 * Math.PI; angle <= 2 * Math.PI; angle += 0.0001) {
			var f = angle.ToFix32();
			swF.Start();
			var actualF = f.SinFastest();
			swF.Stop();
			swD.Start();
			var expected = Math.Sin(angle);
			swD.Stop();
			Assert.AreEqual(expected, actualF.ToDouble(), 16 * Fix32Ext.Precision, "Sin(" + angle + ")" + Delta(expected, actualF.ToDouble(), deltas));
		}

		foreach (var val in TestCases) {
			var f = (Fix32) val;
			swF.Start();
			var actualF = f.SinFastest();
			swF.Stop();
			swD.Start();
			var expected = Math.Sin(f.ToDouble());
			swD.Stop();
			Assert.AreEqual(expected, actualF.ToDouble(), 16 * Fix32Ext.Precision, "Sin(" + f.ToDouble() + ")" + Delta(expected, actualF.ToDouble(), deltas));
		}
		PrintStatistics(deltas, swF, swD);
	}
	
	[Test]
	public void T025_Cos() {
		PrepareStatistics(out var deltas, out var swF, out var swD);

		Assert.AreEqual(Fix32.Zero.Cos().ToDouble(), Fix32.One.ToDouble(), Fix32Ext.Precision);

		Assert.AreEqual(Fix32.PiOver2.Cos().ToDouble(), Fix32.Zero.ToDouble(), Fix32Ext.Precision);
		Assert.AreEqual(Fix32.Pi.Cos().ToDouble(), Fix32.MinusOne.ToDouble(), Fix32Ext.Precision);
		Assert.AreEqual((Fix32.Pi.Add(Fix32.PiOver2)).Cos().ToDouble(), Fix32.Zero.ToDouble(), Fix32Ext.Precision);
		Assert.AreEqual(Fix32.PiTimes2.Cos().ToDouble(), Fix32.One.ToDouble(), Fix32Ext.Precision);

		Assert.AreEqual(Fix32.PiOver2.Neg().Cos().ToDouble(), Fix32.Zero.ToDouble(), Fix32Ext.Precision);
		Assert.AreEqual(Fix32.Pi.Neg().Cos().ToDouble(), Fix32.MinusOne.ToDouble(), Fix32Ext.Precision);
		Assert.AreEqual((Fix32.Pi).Neg().Sub(Fix32.PiOver2).Cos().ToDouble(), Fix32.Zero.ToDouble(), Fix32Ext.Precision);
		Assert.AreEqual(Fix32.PiTimes2.Neg().Cos().ToDouble(), Fix32.One.ToDouble(), Fix32Ext.Precision);


		for (double angle = -2 * Math.PI; angle <= 2 * Math.PI; angle += 0.0001) {
			var f = angle.ToFix32();
			swF.Start();
			var actualF = f.Cos();
			swF.Stop();
			swD.Start();
			var expected = Math.Cos(angle);
			swD.Stop();
			Assert.AreEqual(expected, actualF.ToDouble(), 4 * Fix32Ext.Precision, "Cos(" + angle + ")" + Delta(expected, actualF.ToDouble(), deltas));
		}

		foreach (var val in TestCases) {
			var f = (Fix32) val;
			swF.Start();
			var actualF = f.Cos();
			swF.Stop();
			swD.Start();
			var expected = Math.Cos(f.ToDouble());
			swD.Stop();
			Assert.AreEqual(expected, actualF.ToDouble(), 4 * Fix32Ext.Precision, "Cos(" + f.ToDouble() + ")" + Delta(expected, actualF.ToDouble(), deltas));
		}
		PrintStatistics(deltas, swF, swD);
	}

	[Test]
	public void T025_CosFast() {
		PrepareStatistics(out var deltas, out var swF, out var swD);

		for (double angle = -2 * Math.PI; angle <= 2 * Math.PI; angle += 0.0001) {
			var f = angle.ToFix32();
			swF.Start();
			var actualF = f.CosFast();
			swF.Stop();
			swD.Start();
			var expected = Math.Cos(angle);
			swD.Stop();
			Assert.AreEqual(expected, actualF.ToDouble(), 4 * Fix32Ext.Precision, "Cos(" + angle + ")" + Delta(expected, actualF.ToDouble(), deltas));
		}

		foreach (var val in TestCases) {
			var f = (Fix32) val;
			swF.Start();
			var actualF = f.CosFast();
			swF.Stop();
			swD.Start();
			var expected = Math.Cos(f.ToDouble());
			swD.Stop();
			Assert.AreEqual(expected, actualF.ToDouble(), 4 * Fix32Ext.Precision, "Cos(" + f.ToDouble() + ")" + Delta(expected, actualF.ToDouble(), deltas));
		}
		PrintStatistics(deltas, swF, swD);
	}

	[Test]
	public void T025_CosFastest() {
		PrepareStatistics(out var deltas, out var swF, out var swD);

		for (double angle = -2 * Math.PI; angle <= 2 * Math.PI; angle += 0.0001) {
			var f = angle.ToFix32();
			swF.Start();
			var actualF = f.CosFastest();
			swF.Stop();
			swD.Start();
			var expected = Math.Cos(angle);
			swD.Stop();
			Assert.AreEqual(expected, actualF.ToDouble(), 16 * Fix32Ext.Precision, "Cos(" + angle + ")" + Delta(expected, actualF.ToDouble(), deltas));
		}

		foreach (var val in TestCases) {
			var f = (Fix32) val;
			swF.Start();
			var actualF = f.CosFastest();
			swF.Stop();
			swD.Start();
			var expected = Math.Cos(f.ToDouble());
			swD.Stop();
			Assert.AreEqual(expected, actualF.ToDouble(), 16 * Fix32Ext.Precision, "Cos(" + f.ToDouble() + ")" + Delta(expected, actualF.ToDouble(), deltas));
		}
		PrintStatistics(deltas, swF, swD);
	}
	/*
	[Test]
	public void Atan() {
		var deltas = new List<decimal>();

		Assert.AreEqual(Fix64.Zero, Fix64.Atan(Fix64.Zero));

		// Precision
		for (var x = -1.0; x < 1.0; x += 0.0001) {
			var xf = (Fix64) x;
			var actual = (decimal) Fix64.Atan(xf);
			var expected = (decimal) Math.Atan((double) xf);
			var delta = Math.Abs(actual - expected);
			deltas.Add(delta);
			Assert.LessOrEqual(delta, 10 * Fix64.Precision, string.Format("Precision: Atan({0}): expected {1} but got {2}", xf, expected, actual));
		}

		// Scalability and edge cases
		foreach (var x in TestCases) {
			var xf = Fix64.FromRaw(x);
			var actual = (decimal) Fix64.Atan(xf);
			var expected = (decimal) Math.Atan((double) xf);
			var delta = Math.Abs(actual - expected);
			deltas.Add(delta);
			Assert.LessOrEqual(delta, 10 * Fix64.Precision, string.Format("Scalability: Atan({0}): expected {1} but got {2}", xf, expected, actual));
		}
		Console.WriteLine("Max error: {0} ({1} times precision)", deltas.Max(), deltas.Max() / Fix64.Precision);
		Console.WriteLine("Average precision: {0} ({1} times precision)", deltas.Average(), deltas.Average() / Fix64.Precision);
	}
	
	[Test]
	public void Atan2() {
		var deltas = new List<decimal>();
		// Identities
		Assert.AreEqual(Fix64.Atan2(Fix64.Zero, -Fix64.One), Fix64.Pi);
		Assert.AreEqual(Fix64.Atan2(Fix64.Zero, Fix64.Zero), Fix64.Zero);
		Assert.AreEqual(Fix64.Atan2(Fix64.Zero, Fix64.One), Fix64.Zero);
		Assert.AreEqual(Fix64.Atan2(Fix64.One, Fix64.Zero), Fix64.PiOver2);
		Assert.AreEqual(Fix64.Atan2(-Fix64.One, Fix64.Zero), -Fix64.PiOver2);

		// Precision
		for (var y = -1.0; y < 1.0; y += 0.01) {
			for (var x = -1.0; x < 1.0; x += 0.01) {
				var yf = (Fix64) y;
				var xf = (Fix64) x;
				var actual = (double) Fix64.Atan2(yf, xf);
				var expected = (double) Math.Atan2((double) yf, (double) xf);
				var delta = Math.Abs(actual - expected);
				deltas.Add((decimal) delta);
				Assert.AreEqual(expected, actual, (double) 0.005, string.Format("Precision: Atan2({0}, {1}): expected {2} but got {3}", yf, xf, expected, actual));
			}
		}

		// Scalability and edge cases
		foreach (var y in TestCases) {
			foreach (var x in TestCases) {
				var yf = Fix64.FromRaw(y);
				var xf = Fix64.FromRaw(x);
				var actual = (double) Fix64.Atan2(yf, xf);
				var expected = (double) Math.Atan2((double) yf, (double) xf);
				var delta = Math.Abs(actual - expected);
				deltas.Add((decimal) delta);
				Assert.AreEqual(expected, actual, (double) 0.005, string.Format("Scalability: Atan2({0}, {1}): expected {2} but got {3}", yf, xf, expected, actual));
			}
		}
		Console.WriteLine("Max error: {0} ({1} times precision)", deltas.Max(), deltas.Max() / Fix64.Precision);
		Console.WriteLine("Average precision: {0} ({1} times precision)", deltas.Average(), deltas.Average() / Fix64.Precision);
	}

	/*
	[Test]
	public void Acos() {
		var deltas = new List<decimal>();

		Assert.AreEqual(Fix64.Zero, Fix64.Acos(Fix64.One));
		Assert.AreEqual(Fix64.PiOver2, Fix64.Acos(Fix64.Zero));
		Assert.AreEqual(Fix64.Pi, Fix64.Acos(-Fix64.One));

		// Precision
		for (var x = -1.0; x < 1.0; x += 0.001) {
			var xf = (Fix64) x;
			var actual = (double) Fix64.Acos(xf);
			var expected = Math.Acos((double) xf);
			var delta = Math.Abs(actual - expected);
			deltas.Add((decimal) delta);
			Assert.AreEqual(expected, actual, 12 * (double) Fix32Ext.Precision, string.Format("Precision: Acos({0}): expected {1} but got {2}", xf, expected, actual));
		}

		for (int i = 0; i < TestCases.Count; ++i) {
			var b = Fix64.FromRaw(TestCases[i]);

			if (b < -Fix64.One || b > Fix64.One) {
				Assert.Throws<ArgumentOutOfRangeException>(() => Fix64.Acos(b));
			}
			else {
				var expected = Math.Acos((double) b);
				var actual = (double) Fix64.Acos(b);
				var delta = Math.Abs(expected - actual);
				deltas.Add((decimal) delta);
				Assert.AreEqual(expected, actual, 16 * (double) Fix32Ext.Precision, string.Format("Acos({0}) = expected {1} but got {2}", b, expected, actual));
			}
		}
		Console.WriteLine("Max error: {0} ({1} times precision)", deltas.Max(), deltas.Max() / Fix64.Precision);
		Console.WriteLine("Average precision: {0} ({1} times precision)", deltas.Average(), deltas.Average() / Fix64.Precision);
	}

	*/

#if UNITY_EDITOR // Only override Console inside Unity
	public static class Console {
		public static void WriteLine(string text) {
			UnityEngine.Debug.Log(text);
		}
		public static void WriteLine(string text, params object[] args) {
			UnityEngine.Debug.LogFormat(text, args);
		}
	}
#endif

	static double Saturate(double v) {
		return Math.Max(Fix32.MinValue.ToDouble(), Math.Min(Fix32.MaxValue.ToDouble(), v));
	}

	static string Delta(double a, double b, List<double> deltas) {
		double delta = Math.Abs(a - b);
		if (deltas != null) deltas.Add(delta);
		return " [Delta: " + delta + "]";
	}

	static void PrepareStatistics(out List<double> deltas, out Stopwatch swF, out Stopwatch swD) {
		deltas = new List<double>();
		swF = new Stopwatch();
		swD = new Stopwatch();
	}

	static void PrintStatistics(List<double> deltas, Stopwatch swF, Stopwatch swD) {
		if (deltas.Count > 0) {
			Console.WriteLine("Delta statistics");
			Console.WriteLine("Max error: {0} ({1} times precision)", deltas.Max(), deltas.Max() / (double) Fix32Ext.Precision);
			Console.WriteLine("Min error: {0} ({1} times precision)", deltas.Min(), deltas.Min() / (double) Fix32Ext.Precision);
			var median = deltas.OrderBy(d => d).ElementAt(deltas.Count / 2);
			Console.WriteLine("Med error: {0} ({1} times precision)", median, median / (double) Fix32Ext.Precision);
			Console.WriteLine("Avg error: {0} ({1} times precision)", deltas.Average(), deltas.Average() / (double) Fix32Ext.Precision);
		}

		Console.WriteLine("Fixed:  {0} ms total", swF.Elapsed.TotalMilliseconds);
		Console.WriteLine("Double: {0} ms total", swD.Elapsed.TotalMilliseconds);
		Console.WriteLine("Fixed is {0} times faster than Double (Greater is better)", swD.Elapsed.TotalMilliseconds / swF.Elapsed.TotalMilliseconds);
	}
}
