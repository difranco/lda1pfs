package net.anthonydifranco.LDA;

public class Sampler {
	private static RandomNumberGenerator mt = new MersenneTwisterFast();
	
	public static final int sampleFrom(double[] p) {
		int index;
		final int L = p.length;
		// compute normalization
		double sum = 0.0;
		for (index = 1; index < L; index++) {
			sum += p[index];
		}
		// all zeros is uniform
		if (sum <= 0) return mt.nextInt(L);
		// draw a double in same range as cdf
		double u = mt.nextDouble() * sum;
		// find which component of cdf was sampled
		double probe = 0.0;
		for (index = 0; index < L - 1; index++) {
			probe += p[index];
			if (u < probe)
				break;
		}
		return index;
	}
	
	public static final int sampleFrom(int[] p) {
		int index;
		final int L = p.length;
		// compute normalization
		double sum = 0.0;
		for (index = 1; index < L; index++) {
			sum += p[index];
		}
		// all zeros is uniform
		if (sum <= 0) return mt.nextInt(L);
		// draw a double in same range as cdf
		double u = mt.nextDouble() * sum;
		// find which component of cdf was sampled
		double probe = 0.0;
		for (index = 0; index < L - 1; index++) {
			probe += p[index];
			if (u < probe)
				break;
		}
		return index;
	}
	
	public static final int sampleFrom(Integer[] p) {
		int index;
		final int L = p.length;
		// compute normalization
		double sum = 0.0;
		for (index = 1; index < L; index++) {
			sum += p[index];
		}
		// all zeros is uniform
		if (sum <= 0) return mt.nextInt(L);
		// draw a double in same range as cdf
		double u = mt.nextDouble() * sum;
		// find which component of cdf was sampled
		double probe = 0.0;
		for (index = 0; index < L - 1; index++) {
			probe += p[index];
			if (u < probe)
				break;
		}
		return index;
	}
}
