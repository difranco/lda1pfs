package net.anthonydifranco.LDA;

import java.io.IOException;
import java.util.Map;

public class LDA1PFSServer implements LDA1PFSServerInterface {

	private LDA1PFS lm;

	public void initialize(int Z, int gibbsSize, int burnin, int samplesize) {
		lm = new LDA1PFS(Z, gibbsSize, burnin, samplesize);
	}

	public void addRow(Object rowID, Map<Object, Integer> row) {
		lm.addRow(rowID, row);
	}

	public double[][] currentPhiEstimate() {
		return lm.currentPhiEstimate();
	}

	public double[] estimateThetaOfRow(Map<Object, Integer> row,
			int miniburnin, int minisamplesize) {
		return lm.estimateThetaOfRow(row, miniburnin, minisamplesize);
	}

	public Map<Object, Integer> getColumnCodes() {
		return lm.getColumnCodes();
	}

	public double[] predictForRow(Map<Object, Integer> row, int miniburnin,
			int minisamplesize) {
		return lm.predictForRow(row, miniburnin, minisamplesize);
	}

	public void removeRow(Map<Object, Integer> row) {
		lm.removeRow(row);
	}

	public void saveState(String filename) throws IOException {
		lm.saveState(filename);
	}
	
	public void loadState(String filename) throws IOException, ClassNotFoundException {
		lm = LDA1PFS.loadState(filename);
	}

	public void updatePhiEstimateByRow(Map<Object, Integer> row,
			int miniburnin, int minisamplesize) {
		lm.updatePhiEstimateByRow(row, miniburnin, minisamplesize);
	}

}
