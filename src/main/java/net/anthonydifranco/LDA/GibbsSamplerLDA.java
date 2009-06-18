package net.anthonydifranco.LDA;

import java.io.Serializable;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import cern.colt.matrix.impl.SparseDoubleMatrix3D;

public class GibbsSamplerLDA implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 3637244151534705048L;
	// I can't index by tuple in Java easily!
	// This is the nearest reasonable approximation.
	SparseDoubleMatrix3D Nrfz;  // This is really an int matrix
	final int Rdim;
	final int Fdim;
	final int Z;
	final DoubleMatrix2D Nrf;  // This is really an int matrix
	int[][] Nrz;
	int[][] Nzf;
	int[] Nz;
	int[] Nr;
	final int N;

	private static RandomNumberGenerator mt = new MersenneTwister();

	public GibbsSamplerLDA(DoubleMatrix2D nrf2, int Z) {
		this.Rdim = nrf2.rows();
		this.Fdim = nrf2.columns();
		this.Nrfz = new SparseDoubleMatrix3D(Rdim, Fdim, Z);
		Nrfz.assign(0);
		this.Nrf = nrf2;
		this.Z = Z;
		Nzf = new int[Z][Fdim];
		Nrz = new int[Rdim][Z];
		Nz = new int[Z];
		Nr = new int[Rdim];
		int tempN = 0;

		
		
		for (int row = 0; row < Rdim; row++)
		{
			for (int f = 0; f < Fdim; f++) {
				for (int e = 0; e < nrf2.getQuick(row, f); e++) {
					int component = mt.nextInt(Z);
					Nrfz.setQuick(row, f, component, Nrfz.getQuick(row, f, component) + 1);
					Nrz[row][component]++;
					Nzf[component][f]++;
					Nr[row]++;
					Nz[component]++;
					tempN++;
				}
			}
		}
		N = tempN;
	}

	public void gibbsIteration(double α, double β) {
		for (int i = 0; i < N; i++) {
			int row = Sampler.sampleFrom(Nr);
			int feature = Sampler.sampleFrom(Nrf.viewRow(row).toArray());
			int component = Sampler.sampleFrom(Nrfz.viewSlice(row).viewRow(feature).toArray());
			while(Nrfz.get(row, feature, component) < 1) {
				row = Sampler.sampleFrom(Nr);
				feature = Sampler.sampleFrom(Nrf.viewRow(row).toArray());
				component = Sampler.sampleFrom(Nrfz.viewSlice(row).viewRow(feature).toArray());
			}
			Nrfz.set(row, feature, component, Nrfz.getQuick(row, feature, component) - 1);
			Nzf[component][feature]--;
			Nrz[row][component]--;
			Nz[component]--;
			Nr[row]--;

			double[] p = new double[Z];
			// compute pdf
			for (int z = 0; z < Z; z++) {
				p[z] = (Nrz[row][z] + α) / (Nr[row] + Z * α)
				* (Nzf[z][feature] + β) / (Nz[z] + Fdim * β);
			}
			
			component = Sampler.sampleFrom(p);

			Nrfz.set(row, feature, component, Nrfz.getQuick(row, feature, component) + 1);
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
		double[][] INrf = {{3,0,0,0},{0,3,0,0},{0,0,3,0},{0,0,0,3}};
		SparseDoubleMatrix2D Nrf = new SparseDoubleMatrix2D(INrf);
		int Z = 4;

		System.out.println("Latent Dirichlet Allocation by Gibbs Sampling.");
		System.out.println("Approximating 4x4 Identity scaled by 3 with 4 components.");
		System.out.println("Permutations of the identity should appear below with high probability.");
		System.out.println("Θ should roughly equal Φ'\n");

		/*
		 * The Dirichlet parameter β was chosen to be constant 0.1 while α = 50/k
		 *  and during multi-corpus inference α was constant 50/(k(s) + k(n))
		 *  (these are the default values in GibbsLDA++).
		 *  - István Bíró, Jácint Szabó, András A. Benczúr
		 */
		double α = 50 / Z;
		double β = .1;

		GibbsSamplerLDA LDA = new GibbsSamplerLDA(Nrf, Z);
		for(int i = 0; i < 100; i++)
			LDA.gibbsIteration(α, β);

		int[][] Nrz = LDA.getNrz();
		int[] Nr = LDA.getNr();
		int[][] Nzf = LDA.getNzf();
		int[] Nz = LDA.getNz();
	}
}
