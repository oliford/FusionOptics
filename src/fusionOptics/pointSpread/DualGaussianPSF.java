package fusionOptics.pointSpread;

import java.text.DecimalFormat;

import otherSupport.RandomManager;


import algorithmrepository.exceptions.NotImplementedException;

/** PSF Characterisation using two 2D gaussians.
 * 
 * This one usually works quite well for the normal blue sports caused
 * dominated by spherical abberation.
 * 
 * The inner gaussian is determined by repeatedly recomputing 
 * the mean and stdev while aggressively throwing out points outside of 
 * them (which are moved to the outer one). 
 * 
 * @author oliford
 */
public class DualGaussianPSF extends PointSpreadFunction{
	private static final long serialVersionUID = 187071881902599409L;
	private int minStatsCount = 20;
	
	GaussianPSF inner, outer;

	public DualGaussianPSF() {
	}
	
	public DualGaussianPSF(int minStatsCount) {
		this.minStatsCount = minStatsCount;
	}
	
	@Override
	public boolean isEmpty() { return (inner == null || inner.isEmpty()) && (outer == null || outer.isEmpty()); }
	//public boolean isEmpty() { return (inner == null || inner.isEmpty() || outer == null || outer.isEmpty()); }
	
	
	@Override
	public void setPoints(double[][] pos, double[][][] E) {
		inner = new GaussianPSF(minStatsCount);
		int n = pos.length;
		
		double posInner[][] = pos.clone();
		double posOuter[][] = new double[n][];
		int nOuter = 0;

		for(int i=0; i < 10; i++){
			inner.setPoints(posInner, E);
			
			for(int j=0; j < pos.length; j++){
				double mx = pos[j][0] - inner.meanX, my = pos[j][1] - inner.meanY;
				double u = (mx * inner.ux + my * inner.uy) / inner.ul; 
				double v = (mx * inner.vx + my * inner.vy) / inner.vl;
				double r2 = u*u + v*v;
				if( posInner[j] != null && r2 > 6 ){ //move points outside 4sigma
					posOuter[j] = posInner[j];
					posInner[j] = null;
					nOuter++;
				}
			}			
		}
		
		outer = new GaussianPSF(minStatsCount);
		outer.setPoints(posOuter, E);
		
	}

	@Override
	public double[][] generatePoints(int nPoints) {
		double pos[][] = new double[2][nPoints];
		double innerProb = inner.I0 / (inner.I0 + outer.I0);
		for(int i=0; i < nPoints; i++){
			double p = RandomManager.instance().nextNormal(0, 1);
			double a = RandomManager.instance().nextNormal(0, 1);
			double b = RandomManager.instance().nextNormal(0, 1);
			if(p < innerProb){
				pos[0][i] = inner.meanX + a*inner.ul*inner.ux + b*inner.vl*inner.vx;
				pos[1][i] = inner.meanY + a*inner.ul*inner.uy + b*inner.vl*inner.vy;
			}else{
				pos[0][i] = outer.meanX + a*outer.ul*outer.ux + b*outer.vl*outer.vx;
				pos[1][i] = outer.meanY + a*outer.ul*outer.uy + b*outer.vl*outer.vy;
			}
		}
		
		return pos;
	}

	@Override
	public void generatePolarisedPoints(int nPoints, double[][] pos,
			double[][][] E) {
		throw new NotImplementedException();
	}

	@Override
	public double[][] addToGrid(double[] x, double[] y, double[][] G, double Imul) {
		inner.addToGrid(x, y, G, Imul);
		outer.addToGrid(x, y, G, Imul);
		return G;
	}

	public double[][] addToGridMath(double[] x, double[] y, double[][] G, double Imul) {
		inner.addToGridMath(x, y, G, Imul);
		outer.addToGridMath(x, y, G, Imul);
		return G;
	}
	
	public double[][] addToGridMonteCarlo(double[] x, double[] y, double[][] G, double Imul, int nPoints) {
		inner.addToGridMonteCarlo(x, y, G, Imul, nPoints / 2);
		outer.addToGridMonteCarlo(x, y, G, Imul, nPoints / 2);
		return G;
	}
	
	@Override
	public double getMinX() { return Math.min(inner.getMinX(), outer.getMinX()); }
	public double getMaxX() { return Math.max(inner.getMaxX(), outer.getMaxX()); }
	public double getMinY() { return Math.min(inner.getMinY(), outer.getMinY()); }
	public double getMaxY() { return Math.max(inner.getMaxY(), outer.getMaxY()); }

	@Override
	public double[] getCharacterisationData() {
		double ret[] = new double[16 + 12];
		System.arraycopy(inner.getCharacterisationData(), 0, ret, 0, 6);
		System.arraycopy(outer.getCharacterisationData(), 0, ret, 6, 6);
		System.arraycopy(meanM, 0, ret, 12, 16);
		/*if(inner.meanM != null){
			System.arraycopy(inner.meanM, 0, ret, 12, 16);
			if(outer.meanM != null){
				for(int i=0; i < 16; i++)
					ret[16+i] += outer.meanM[i];
			}
		}else if(outer.meanM != null){
			System.arraycopy(outer.meanM, 0, ret, 12, 16);
		}*/
		return ret;
	}
	
	@Override
	public void setCharacterisationData(double[] data) {
		double innerData[] = new double[6];
		double outerData[] = new double[6];
		meanM = new double[16];
		System.arraycopy(data, 0, innerData, 0, 6);
		System.arraycopy(data, 6, outerData, 0, 6);
		System.arraycopy(data, 12, meanM, 0, 16);

		inner = new GaussianPSF();
		inner.setCharacterisationData(innerData);
		outer = new GaussianPSF();
		outer.setCharacterisationData(outerData);
	}
	
	@Override
	public void combine(PointSpreadFunction[] psfs, double[] coeffs) {
		GaussianPSF inners[] = new GaussianPSF[psfs.length];
		GaussianPSF outers[] = new GaussianPSF[psfs.length];

		for(int i=0; i < psfs.length; i++){
			if(coeffs[i] == 0)
				continue;
			inners[i] = ((DualGaussianPSF)psfs[i]).inner;
			outers[i] = ((DualGaussianPSF)psfs[i]).outer;			
		}
		inner = new GaussianPSF();
		inner.combine(inners, coeffs);
		outer = new GaussianPSF();
		outer.combine(outers, coeffs);
		
		setSourceFrom(psfs, coeffs);
		setMuellerFrom(psfs,coeffs);
	}

	@Override
	public void multiplyIntensity(double Imul) {
		inner.multiplyIntensity(Imul);
		outer.multiplyIntensity(Imul);		
	}

	@Override
	public double getIntensity() {
		return inner.getIntensity() + outer.getIntensity();
	}

	public GaussianPSF getInner() { return inner; }
	public GaussianPSF getOuter() { return outer; }
	
}
