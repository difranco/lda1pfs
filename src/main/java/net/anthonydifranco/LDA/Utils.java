package net.anthonydifranco.LDA;

public class Utils {
	public static double[][] arrayCastDouble(int[][] a) {
		double[][] t = new double[a.length][a[0].length];
		for(int i = 0; i < a.length; i++) {
			for(int j = 0; j < a[0].length; j++) {
				t[i][j] = (double)(a[i][j]);
			}
		}
		return t;
	}
}
