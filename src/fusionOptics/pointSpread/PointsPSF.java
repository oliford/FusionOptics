package fusionOptics.pointSpread;

import java.io.Serializable;

import fusionOptics.types.Pol;

import uk.co.oliford.jolu.OneLiners;
import uk.co.oliford.jolu.RandomManager;
import algorithmrepository.Algorithms;

/** Represents the PSF directly with the originally collected points */
public class PointsPSF extends PointSpreadFunction {
	private static final long serialVersionUID = 691093674356642058L;

	private double pos[][], E[][][], I[];
	private double cdf[];
	private double minX, maxX, minY, maxY;
		
	public PointsPSF() {		
	}
	
	@Override
	public double[] getCharacterisationData() {
		return new double[0];
	}

	@Override
	public void setCharacterisationData(double[] data) {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public boolean isEmpty() { return pos == null || pos.length == 0; }
	
	@Override
	public void setPoints(double[][] pos, double[][][] E) {
		this.pos = pos;
		this.E = E;
		int n = pos.length;
		I = new double[n];
		
		cdf = new double[n+1];
		for(int i=0; i < n; i++) {
			I[i] = 0;
			for(int j=0; j < E[i].length; j++)
				I[i] += Pol.intensity(E[i][j]);
			
			cdf[i+1] = cdf[i] + I[i];
		}
		
		for(int i=0; i < cdf.length; i++)
			cdf[i] /= cdf[n];
		
		if(pos.length == 0){
			minX = -1; maxX = 1;
			minY = -1; maxY = 1;
			return;
		}
		
		double ret[][] = OneLiners.columnRanges(pos);
		minX = ret[0][0]; maxX = ret[1][0];
		minY = ret[0][1]; maxY = ret[1][1];
	}

	@Override
	public double[][] generatePoints(int nPoints) {
		double ret[][] = new double[2][nPoints];
	
		//pick points randomly biased by intensity
		for(int i=0; i < nPoints; i++){
			double a = RandomManager.instance().nextUniform(0, 1);
			int idx = OneLiners.getNearestLowerIndex(cdf, a);
			ret[0][i] = pos[idx][0];
			ret[1][i] = pos[idx][1];
		}
		
		return ret;
	}
	
	@Override
	public void generatePolarisedPoints(int nPoints, double[][] pos, double[][][] E) {
		
		//pick points randomly biased by intensity
		for(int i=0; i < nPoints; i++){
			double a = RandomManager.instance().nextUniform(0, 1);
			int idx = OneLiners.getNearestLowerIndex(cdf, a);
			pos[i] = this.pos[idx];
			E[i] = this.E[idx];
		}	
	}

	@Override
	public double[][] addToGrid(double[] x, double[] y, double[][] G, double I0) {
		for(int i=0; i < pos.length; i++) {
			int iX = OneLiners.getNearestLowerIndex(x, pos[i][0]);
			int iY = OneLiners.getNearestLowerIndex(y, pos[i][1]);
			
			if(iX > 0 && iY > 0 && iX < x.length && iY < y.length){
				G[iY][iX] += I0 * I[i];
			}
		}
		
		return G;
	}

	@Override
	public double[][] addToGridMonteCarlo(double[] x, double[] y, double[][] G,
			double I, int nPoints) {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public double getMinX() { return minX;	}
	@Override
	public double getMaxX() { return maxX;	}
	@Override
	public double getMinY() { return minY;	}
	@Override
	public double getMaxY() { return maxY;	}

	@Override
	public void multiplyIntensity(double Imul) {
		for(int i=0; i < I.length; i++)
			I[i] *= Imul;
	}

	@Override
	public double getIntensity() {
		double sum = 0;
		for(int i=0; i < I.length; i++)
			sum += I[i];
		return sum;
	}
}
