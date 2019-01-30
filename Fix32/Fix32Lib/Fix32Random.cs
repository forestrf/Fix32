using System;

public class Fix32Random {
	private Random random;

	public Fix32Random(int seed) {
		random = new Random(seed);
	}

	public Fix32 Next() {
		return (Fix32) random.Next(int.MinValue, int.MaxValue);
	}

	public Fix32 NextInt(int maxValue) {
		return random.Next(maxValue).ToFix32();
	}
}
