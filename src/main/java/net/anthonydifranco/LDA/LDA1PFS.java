package net.anthonydifranco.LDA;

import java.util.ArrayList;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.lang.Math;

import cern.colt.list.DoubleArrayList;
import cern.colt.list.IntArrayList;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.SparseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import cern.colt.matrix.DoubleFactory2D;
import cern.jet.math.Functions;
import cern.jet.math.PlusMult;
import cern.colt.matrix.linalg.SeqBlas;
import cern.jet.random.Normal;

public class LDA1PFS implements Serializable {
	private static final long serialVersionUID = 4782077567351681783L;
	Encoder rowCodes;
	Encoder colCodes;

	ArrayList<TreeMap<Integer, Integer>> NrfList;
	DoubleMatrix2D Nrf;  // This is really an int matrix
	int[] Nf;

	DoubleMatrix2D[] Nzfsamples;  // This is really an int matrix
	double[] Nzfsampleweights;
	DoubleMatrix2D phihat;
	DoubleMatrix2D thetahat;
	final double α, β;
	final int Z;
	final int gibbsSize;
	final int burnin;
	final int samplesize;

	static RandomNumberGenerator mt;
	GibbsSamplerLDA gsm;

	public LDA1PFS() {
		Z = 20;
		α = 50 / Z;
		β = 0.1;
		gibbsSize = 200;
		burnin = 20;
		samplesize = 5;	
	}

	public LDA1PFS(int Z, int gibbsSize, int burnin, int samplesize) {
		rowCodes = new Encoder();
		colCodes = new Encoder();
		α = 50 / Z;
		β = 0.1;
		this.Z = Z;
		this.NrfList = new ArrayList<TreeMap<Integer,Integer>>();
		this.gibbsSize = gibbsSize;
		this.burnin = burnin;
		this.samplesize = samplesize;
		LDA1PFS.mt = new MersenneTwisterFast();
	}

	public void addRow(Object rowID, Map<Object, Integer> row) {
		if(rowCodes.encode.containsKey(rowID)) {
			final int rownum = rowCodes.encoded(rowID);
			if(rownum > gibbsSize) {
				removeRow(row);
				updatePhiEstimateByRow(row, 20, 5);
			} else {
				;
			}
		} else {
			// new row
			final int rownum = rowCodes.encoded(rowID);
			if (rownum <= gibbsSize) {
				// not enough rows to do Gibbs yet, just record it and wait
				TreeMap<Integer, Integer> codemap = new TreeMap<Integer,Integer>();
				for (Entry<Object, Integer> e: row.entrySet()) {
					final int colnum = colCodes.encoded(e.getKey());
					int frequency = e.getValue();
					codemap.put(colnum, frequency);
				}
				NrfList.add(codemap);
			}
			if (rownum == gibbsSize) {
				// time for Gibbs, convert stored rows to matrix and run Gibbs
				gibbsSampleSubset();
			}
			if (rownum > gibbsSize) {
				// Gibbs is over, just maintain phi estimate
				updatePhiEstimateByRow(row, 12, 4);
			}
		}
	}

	private void gibbsSampleSubset() {
		Nrf = new SparseDoubleMatrix2D(rowCodes.codesize, colCodes.codesize);
		Nzfsamples = new SparseDoubleMatrix2D[samplesize];
		Nzfsampleweights = new double[samplesize];
		Nf = new int[colCodes.codesize];
		for(int rowindex = NrfList.size() - 1; rowindex >= 0; rowindex--) {
			for(Entry<Integer, Integer> entry: NrfList.get(rowindex).entrySet()) {
				int colindex = entry.getKey();
				int frequency = entry.getValue();
				Nrf.setQuick(rowindex, colindex, frequency);
				Nf[colindex] += frequency;
			}
			NrfList.remove(rowindex);
		}
		gsm = new GibbsSamplerLDA(Nrf, Z);
		phihat = new DenseDoubleMatrix2D(Z, colCodes.codesize);
		phihat.assign(0);
		thetahat = new DenseDoubleMatrix2D(Nrf.rows(), Z);
		thetahat.assign(0);
		for (int i = 0; i < burnin; i++) gsm.gibbsIteration(α, β);
		for (int i = 0; i < samplesize; i++) {
			gsm.gibbsIteration(α, β);
			Nzfsamples[i] = new SparseDoubleMatrix2D(Z, colCodes.codesize);
			Nzfsamples[i].assign(Utils.arrayCastDouble(gsm.getNzf()));
			Nzfsampleweights[i] = 1.0 / samplesize;
			phihat.assign(Nzfsamples[i], PlusMult.plusDiv(samplesize));
			thetahat.assign(new DenseDoubleMatrix2D(Utils.arrayCastDouble(gsm.getNrz())), PlusMult.plusDiv(samplesize));
		}
		for(int z = 0; z < Z; z++) {
			for(int f = 0; f < colCodes.codesize; f++) {
				phihat.setQuick(z, f, (phihat.getQuick(z, f) + β) / (Nf[f] + colCodes.codesize * β));
			}
		}
		for(int r = 0; r <= gibbsSize; r++) {
			int Nr = 0;
			for(int z = 0; z < Z; z++) Nr += thetahat.getQuick(r,z);
			for(int z = 0; z < Z; z++) {
				thetahat.setQuick(r, z, (thetahat.getQuick(r, z) + α) / (Nr + Z * α));
			}
		}
		//System.out.println("α: " + α);
		//System.out.println("Gibbs Phihat\n" + ArrayPrettyPrinter.mat2str(phihat.toArray()));
		//System.out.println("Gibbs Thetahat\n" + ArrayPrettyPrinter.mat2str(thetahat.toArray()));
		gsm = null;
	}

	private DoubleMatrix2D[] sampleZFofRow(Map<Object, Integer> row, int miniburnin, int minisamplesize) {
		// this samples Z assignments for a row from the global phihat WITHOUT updating phihat accordingly
		DoubleMatrix2D Nzf; // This is really an int matrix
		int[] Nrz = new int[Z];
		int Nr = 0;
		DoubleMatrix2D[] miniNzfsamples = new DoubleMatrix2D[minisamplesize];
		miniNzfsamples[0] = new SparseDoubleMatrix2D(Z, colCodes.codesize);
		Nzf = miniNzfsamples[0];
		DoubleMatrix1D codedrow = new SparseDoubleMatrix1D(colCodes.codesize);  // this is really an int vector

		for(Entry<Object, Integer> e: row.entrySet()) {
			Object key = e.getKey();
			// only take features that are in the model (phi is defined)
			if(colCodes.encode.containsKey(key)) {
				int col = colCodes.encode.get(key);
				Integer val = e.getValue();
				codedrow.setQuick(col, val);
				for(int i = 0; i < val; i++) {
					int z = mt.nextInt(Z);  // currently I'm killing this assignment in first iteration below
					Nzf.setQuick(z, col, Nzf.getQuick(z, col) + 1.0);
					Nrz[z]++;
				}
				Nr += val;
			}
		}

		DoubleMatrix1D p = new DenseDoubleMatrix1D(Z);
		// Smoothing for the below will be a little off due to different N,
		// but I'm renormalizing when I resample anyway.
		// Do I want p inverse here or this, p?:  p is somewhat right stochastic; normalization doesn't matter here.
		SeqBlas.seqBlas.dgemv(false, 1, phihat, codedrow, 0, p);

		for (int i = 0; i < miniburnin + minisamplesize; i++) {
			DoubleMatrix1D pi = p.copy();
			if(i <= miniburnin) {
				Nzf = miniNzfsamples[0];
			} else {
				miniNzfsamples[i - miniburnin] = miniNzfsamples[i - miniburnin - 1].copy();
				Nzf = miniNzfsamples[i - miniburnin];
			}
			for (int z = 0; z < Z; z++) if (i > 0)
				pi.setQuick(z, pi.getQuick(z) * (Nrz[z] + α) / (Nr + Z * α));
			IntArrayList zs = new IntArrayList(Nr), fs = new IntArrayList(Nr);
			DoubleArrayList counts = new DoubleArrayList(Nr);  // this is really an IntArrayList
			Nzf.getNonZeros(zs, fs, counts);
			for(int observation = 0; observation < counts.size(); observation++) {
				int z = zs.getQuick(observation);
				int f = fs.getQuick(observation);
				double count = counts.getQuick(observation);
				for(int c = 0; c < count; c++) {
					Nzf.setQuick(z, f, Nzf.getQuick(z, f) - 1.0);
					Nrz[z]--;
					z = Sampler.sampleFrom(pi.toArray());
					Nzf.setQuick(z, f, Nzf.getQuick(z, f) + 1.0);
					Nrz[z]++;
				}
			}
		}
		return miniNzfsamples;
	}

	public void removeRow(Map<Object, Integer> row) {
		// decrements Nzf sample counts and updates phihat by sampling Z for this row
		// (sample a Z sample uniformly, then sample an element,
		// threshold rather than push any Nzf entry below zero) N times.
		// (thresholding should happen with very low probability)
		// Then, update phi with an empty row to run the particle filter code with the new counts.

		DoubleMatrix2D[] miniNzfsamples = sampleZFofRow(row, 20, 5);

		ArrayList<Integer> fset = new ArrayList<Integer>(row.size());

		for(Entry<Object, Integer> e: row.entrySet()) {
			// only take features that are in the model (phi is defined)
			Object key = e.getKey();
			if(colCodes.encode.containsKey(key)) {
				int f = colCodes.encode.get(key);
				fset.add(f);
				// Decrement Nf counts
				Nf[f] -= e.getValue();
			}
		}
		// decrement Nzf counts
		for (int sample = 0; sample < samplesize; sample++) {
			for (int f: fset) {
				int minisample = mt.nextInt(5);
				for (int z = 0; z < Z; z++) {
					int difference = (int) miniNzfsamples[minisample].getQuick(z, f);
					Nzfsamples[sample].setQuick(z, f, Math.max(Nzfsamples[sample].getQuick(z, f) - difference, 0));
				}
			}
		}

		updatePhiEstimateByRow(new TreeMap<Object, Integer>(), 20, 5);
	}

	public void updatePhiEstimateByRow(Map<Object, Integer> row, int miniburnin, int minisamplesize) {
		// updates the global phihat and each of the Nzf samples by sampling Z for a given row
		// updates importance weights and resamples from kernel estimator if needed
		// (sample a Z sample uniformly, then sample an element) N times.
		DoubleMatrix1D codedrow = new SparseDoubleMatrix1D(colCodes.codesize);
		ArrayList<Integer> fset = new ArrayList<Integer>(row.size());
		TreeMap<Integer, Integer> newfs = new TreeMap<Integer, Integer>();

		for(Entry<Object, Integer> e: row.entrySet()) {
			Object key = e.getKey();
			Integer val = e.getValue();
			int f;
			// only take features that are in the model (phi is defined)
			if(colCodes.encode.containsKey(key)) {
				f = colCodes.encode.get(key);
				codedrow.setQuick(f,val);
			} else {
				// extend column code
				f = colCodes.encoded(key);
				newfs.put(f, val);
			}
			fset.add(f);
		}
		if ( ! newfs.isEmpty() ) {
			// extend Nf
			int[] newNf = new int[Nf.length + newfs.size()];
			java.lang.System.arraycopy(Nf, 0, newNf, 0, Nf.length);
			for (Entry<Integer, Integer> nf: newfs.entrySet()) {
				newNf[nf.getKey()] = nf.getValue();
			}
			Nf = newNf;
			// extend each of Nzfsamples
			for(int sample = 0; sample < Nzfsamples.length; sample++) {
				DoubleMatrix2D newcols = new SparseDoubleMatrix2D(Z, newfs.size());
				Nzfsamples[sample] = DoubleFactory2D.sparse.appendColumns(Nzfsamples[sample], newcols);
			}
			// extend phihat
			DoubleMatrix2D newphihatcols = new DenseDoubleMatrix2D(Z, newfs.size());
			phihat = DoubleFactory2D.dense.appendColumns(phihat, newphihatcols);
		}

		DoubleMatrix2D[] miniNzfsamples = sampleZFofRow(row, miniburnin, minisamplesize);

		// Kind of a move-before-resample here (Gilks, Berzuini) updating Nzf counts
		for (int sample = 0; sample < samplesize; sample++) {
			for (int f: fset) {
				int minisample = mt.nextInt(minisamplesize);
				for (int z = 0; z < Z; z++) {
					int difference = (int) miniNzfsamples[minisample].getQuick(z, f);
					Nzfsamples[sample].setQuick(z, f, Nzfsamples[sample].getQuick(z, f) + difference);
				}
			}
		}

		// Update Nf counts
		for(Entry<Object, Integer> e: row.entrySet()) {
			Nf[colCodes.encode.get(e.getKey())] += e.getValue();
		}

		// Update importance weights
		DoubleMatrix1D thetarhat = new DenseDoubleMatrix1D(estimateThetaOfRow(row, miniburnin, minisamplesize));
		DoubleMatrix1D fprobs = new SparseDoubleMatrix1D(colCodes.codesize);
		for (int sample = 0; sample < samplesize; sample++) {
			DoubleMatrix2D phishat = new SparseDoubleMatrix2D(Z, colCodes.codesize);
			for (int z = 0; z < Z; z++) {
				for (int f: fset) {
					phishat.setQuick(z, f, (Nzfsamples[sample].getQuick(z, f) + β) / (Nf[f] + colCodes.codesize * β));
				}
			}
			SeqBlas.seqBlas.dgemv(true, 1, phishat, thetarhat, 0, fprobs);
			double likelihood = SeqBlas.seqBlas.ddot(fprobs, codedrow);
			Nzfsampleweights[sample] = Nzfsampleweights[sample] * (likelihood + 10e-8);
		}
		// renormalize importance weights
		double importanceweightsum = 0;
		for (int sample = 0; sample < samplesize; sample++) importanceweightsum += Nzfsampleweights[sample];
		for (int sample = 0; sample < samplesize; sample++) Nzfsampleweights[sample] /= importanceweightsum;
		//for (int sample = 0; sample < samplesize; sample++) System.out.println("postweight " + Nzfsampleweights[sample]);

		// Compute Effective Sample Size
		importanceweightsum = 0;
		double importanceweightsquaresum = 0;
		for (int i = 0; i < samplesize; i++) {
			double weight = Nzfsampleweights[i];
			importanceweightsum += weight;
			importanceweightsquaresum += weight * weight;
		}

		double ESS = importanceweightsum * importanceweightsum / importanceweightsquaresum;

		if (ESS < 0.25 * samplesize) {
			System.out.print("Replenishing sample set...");
			// resample with replacement from sample set, importance weighted
			DoubleMatrix2D[] newNzfsamples = new DoubleMatrix2D[samplesize];
			double[] newsampleweights = new double[samplesize];
			for (int i = 0; i < samplesize; i++) {
				int sample = Sampler.sampleFrom(Nzfsampleweights);
				newNzfsamples[i] = Nzfsamples[sample].copy();
				newsampleweights[i] = 1.0 / samplesize;
			}
			// follow Stavropoulos and Titterington 2001
			// importance-weighted variance of sample components from importance-weighted means of components
			DoubleMatrix2D meanNzf = new DenseDoubleMatrix2D(Z, colCodes.codesize);
			DoubleMatrix2D bNzf = new DenseDoubleMatrix2D(Z, colCodes.codesize); // Gaussian kernel bandwidth
			for (int z = 0; z < Z; z++) {
				for (int f = 0; f < colCodes.codesize; f++) {
					for (int i = 0; i < samplesize; i++) {
						meanNzf.setQuick(z, f, meanNzf.getQuick(z, f) + newNzfsamples[i].getQuick(z, f));
					}
					meanNzf.setQuick(z, f, meanNzf.getQuick(z, f) / samplesize);
				}
			}
			double bN = Math.pow(4.0 / ((colCodes.codesize * Z + 2.0) * samplesize), 1.0 / (colCodes.codesize * Z + 4.0));
			for (int z = 0; z < Z; z++) {
				for (int f = 0; f < colCodes.codesize; f++) {
					for (int i = 0; i < samplesize; i++) {
						bNzf.setQuick(z, f, bNzf.getQuick(z,f) + newNzfsamples[i].getQuick(z, f) - meanNzf.getQuick(z, f));
					}
					double s = bNzf.getQuick(z, f);
					bNzf.setQuick(z, f, s * s * bN * bN);
				}
			}
			double a = Math.sqrt(1 - bN * bN);
			DoubleMatrix2D[] shrunkNzfsamples = new DoubleMatrix2D[samplesize];
			for (int i = 0; i < samplesize; i++) {
				shrunkNzfsamples[i] = meanNzf.copy();
				SeqBlas.seqBlas.dscal(1 - a, shrunkNzfsamples[i]);
				SeqBlas.seqBlas.daxpy(a, Nzfsamples[i], shrunkNzfsamples[i]);
				for (int z = 0; z < Z; z++) {
					for (int f = 0; f < colCodes.codesize; f++) {
						double sigmasquare = bNzf.getQuick(z, f);
						double mu = shrunkNzfsamples[i].getQuick(z, f);
						shrunkNzfsamples[i].setQuick(z, f, (int) Math.max(Normal.staticNextDouble(mu, sigmasquare), 0));
					}
				}
			}
			Nzfsamples = shrunkNzfsamples;
			Nzfsampleweights = newsampleweights;
			// renormalize importance weights
			importanceweightsum = 0;
			for (int sample = 0; sample < samplesize; sample++) importanceweightsum += Nzfsampleweights[sample];
			for (int sample = 0; sample < samplesize; sample++) Nzfsampleweights[sample] /= importanceweightsum;
			//for (int sample = 0; sample < samplesize; sample++) System.out.println("ESSpostweight " + Nzfsampleweights[sample]);
			System.out.println("done.");
		}
		// Need to actually update phihat now
		for(int z = 0; z < Z; z++) {
			for(int f = 0; f < colCodes.codesize; f++) {
				double Nzfsum = 0;
				for (int i = 0; i < samplesize; i++) {
					Nzfsum += Nzfsamples[i].getQuick(z, f) * Nzfsampleweights[i]; 
				}
				phihat.setQuick(z, f, (Nzfsum + β) / (Nf[f] + colCodes.codesize * β));
			}
		}
	}

	public double[][] currentPhiEstimate() {
		return phihat.toArray();
	}

	public double[] estimateThetaOfRow(Map<Object, Integer> row, int miniburnin, int minisamplesize) {
		DoubleMatrix2D[] Nzfsamples = sampleZFofRow(row, miniburnin, minisamplesize);
		double[] thetaRhat = new double[Z];
		double Nr = 0;  // theoretically an int but we need to avoid discretization issues in averaging samples
		for (int minisample = 0; minisample < minisamplesize; minisample++)
			for (int z = 0; z < Z; z++) {
				Nr += Nzfsamples[minisample].viewRow(z).aggregate(Functions.plus, Functions.identity);
			}
		Nr = Nr / (double) minisamplesize;
		for (int z = 0; z < Z; z++) {
			for (int minisample = 0; minisample < minisamplesize; minisample++) {
				thetaRhat[z] += Nzfsamples[minisample].viewRow(z).aggregate(Functions.plus, Functions.identity);
			}
			thetaRhat[z] = (thetaRhat[z] / minisamplesize + α) / (Nr + Z * α);
		}
		return thetaRhat;
	}

	public double[] predictForRow(Map<Object, Integer> row, int miniburnin, int minisamplesize) {
		DoubleMatrix2D thetarhat = new DenseDoubleMatrix2D(1, Z);
		thetarhat.viewRow(0).assign(estimateThetaOfRow(row, miniburnin, minisamplesize));
		DoubleMatrix2D fprobs = new DenseDoubleMatrix2D(1, colCodes.codesize);
		SeqBlas.seqBlas.dgemm(false, false, 1, thetarhat, phihat, 0, fprobs);
		double sum = fprobs.aggregate(Functions.plus, Functions.identity);
		return fprobs.assign(Functions.mult(1 / sum)).toArray()[0];
	}

	public Map<Object, Integer> getColumnCodes() {
		return new TreeMap<Object, Integer>(colCodes.encode);
	}

	// Serialization methods
	public void saveState(String filename) throws IOException {
		FileOutputStream fos = new FileOutputStream(filename);
		ObjectOutputStream out = new ObjectOutputStream(fos);
		out.writeObject(this);
	}

	public static LDA1PFS loadState(String filename) throws IOException, ClassNotFoundException {
		FileInputStream fis = new FileInputStream(filename);
		ObjectInputStream in = new ObjectInputStream(fis);
		return (LDA1PFS) in.readObject();
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// test driver                                                                                                                                        
		System.out.println("small matrix scaled by 20");
		LDA1PFS lm = new LDA1PFS(2, 20, 100, 20);
		Map evenrow = new TreeMap();
		Map oddrow = new TreeMap();
		evenrow.put(0, 20); evenrow.put(1, 0); evenrow.put(2, 20); evenrow.put(3, 0);
		oddrow.put(0, 0); oddrow.put(1, 20); oddrow.put(2, 0); oddrow.put(3, 20);
		for(Integer i = 0; i < 100; i++) {
			if (i % 2 == 0) {
				lm.addRow(i, evenrow);
			}
			else {
				lm.addRow(i, oddrow);
			}
		}
		double[][] p = lm.currentPhiEstimate();
		double[][] thetaz = new double[2][];
		thetaz[0] = lm.estimateThetaOfRow(evenrow, 1000, 10);
		thetaz[1] = lm.estimateThetaOfRow(oddrow, 1000, 10);

		System.out.println("bigger test, repeated identity(6), usually replenishes");
		System.out.println("also should hopefully enforce sparsity by only finding 6 / 12 components");
		Map row = null;
		lm = new LDA1PFS(12, 24, 1000, 20);
		for(int i = 0; i < 1000; i++) {
			row = new TreeMap<Integer, Integer>();
			row.put(i % 6, 5);
			lm.addRow(i, row);
		}
		for(int i = 0; i < 1000; i++) {
			// This will delete rows
			row = new TreeMap<Integer, Integer>();
			row.put(i % 6, 5);
			lm.addRow(i, row);
		}

		p = lm.currentPhiEstimate();
		double[][] q = new double[1][];
		q[0] = lm.predictForRow(row, 60, 5);
		try { lm.saveState("testsave"); }
		catch (Exception e) { System.out.println("Didn't save state!"); };
		try { lm = LDA1PFS.loadState("testsave");}
		catch (Exception e) { System.out.println("Didn't load state!"); }
		q[0] = lm.predictForRow(row, 60, 5);
	}
}
