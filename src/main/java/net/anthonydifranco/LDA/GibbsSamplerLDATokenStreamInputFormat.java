// LDA by Gibbs Sampler
// Anthony Di Franco

package net.anthonydifranco.LDA;

public class GibbsSamplerLDATokenStreamInputFormat {
	final int[] R, F;
	int[] Z;
	final int Rdim;
	final int Fdim;
	final int C;
	int[][] Nrz;
	int[][] Nzf;
	int[] Nz;
	int[] Nr;

	private static RandomNumberGenerator mt = new MersenneTwisterFast();

	public GibbsSamplerLDATokenStreamInputFormat(int[] R, int Rdim, int[] F, int Fdim, int C) {
		this.R = R;
		this.F = F;
		this.Z = new int[R.length];
		this.Rdim = Rdim;
		this.Fdim = Fdim;
		this.C = C;
		Nzf = new int[C][Fdim];
		Nrz = new int[Rdim][C];
		Nz = new int[C];
		Nr = new int[Rdim];

		for (int i = 0; i < R.length; i++)
		{
			int row = R[i];
			int feature = F[i];
			int component = mt.nextInt(C);
			Z[i] = component;
			Nrz[row][component]++;
			Nzf[component][feature]++;
			Nr[row]++;
			Nz[component]++;
		}
	}

	public void gibbsIteration(double α, double β) {
		for (int i = 0; i < R.length; i++) {
			int row = R[i];
			int feature = F[i];
			int component = Z[i];
			Nzf[component][feature]--;
			Nrz[row][component]--;
			Nz[component]--;
			Nr[row]--;

			double[] p = new double[C];
			// compute unnormalized pdf
			for (int k = 0; k < C; k++) {
				p[k] = (Nzf[k][feature] + β) / (Nz[k] + Fdim * β)
				* (Nrz[row][k] + α) / (Nr[row] + C * α);
			}

			component = Sampler.sampleFrom(p);
			
			Z[i] = component;
			Nzf[component][feature]++;
			Nrz[row][component]++;
			Nz[component]++;
			Nr[row]++;
		}
	}

	public int[][] getNrz() {
		return Nrz;
	}

	public int[] getNr() {
		return Nr;
	}

	public int[][] getNzf() {
		return Nzf;
	}

	public int[] getNz() {
		return Nz;
	}

	public static void main(String[] args) {
		int[] R = {0,0,0,1,1,1,2,2,2,3,3,3};
		int Rdim = 4;
		int[] F = {0,0,0,1,1,1,2,2,2,3,3,3};
		int Fdim = 4;
		int C = 4;

		System.out.println("Latent Dirichlet Allocation by Gibbs Sampling.");
		System.out.println("Approximating 4x4 Identity scaled by 3 with 4 components.");
		System.out.println("Permutations of the identity should appear below with high probability.\n");

		/*
		 * "The Dirichlet parameter β was chosen to be constant 0.1 while α = 50/k
		 *  and during multi-corpus inference α was constant 50/(k(s) + k(n))
		 *  (these are the default values in GibbsLDA++)."
		 *  - István Bíró, Jácint Szabó, András A. Benczúr
		 */
		double α = 50 / C;
		double β = .1;

		GibbsSamplerLDATokenStreamInputFormat LDA = new GibbsSamplerLDATokenStreamInputFormat(R, Rdim, F, Fdim, C);
		for(int i = 0; i < 100; i++)
			LDA.gibbsIteration(α, β);

		int[][] Nrz = LDA.getNrz();
		int[] Nr = LDA.getNr();
		int[][] Nfz = LDA.getNzf();
		int[] Nz = LDA.getNz();
	}
}
