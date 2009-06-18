package net.anthonydifranco.LDA;

import java.io.IOException;
import java.util.Map;

public interface LDA1PFSServerInterface {

	public abstract void initialize(int Z, int gibbsSize, int burnin, int samplesize);
	
	public abstract void addRow(Object rowID, Map<Object, Integer> row);

	public abstract void removeRow(Map<Object, Integer> row);

	public abstract void updatePhiEstimateByRow(Map<Object, Integer> row,
			int miniburnin, int minisamplesize);

	public abstract double[][] currentPhiEstimate();

	public abstract double[] estimateThetaOfRow(Map<Object, Integer> row,
			int miniburnin, int minisamplesize);

	public abstract double[] predictForRow(Map<Object, Integer> row,
			int miniburnin, int minisamplesize);

	public Map<Object, Integer> getColumnCodes();
	
	// Serialization methods
	public abstract void saveState(String filename) throws IOException;

}