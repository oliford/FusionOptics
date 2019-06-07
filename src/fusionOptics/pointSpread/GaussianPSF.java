package fusionOptics.pointSpread;

import java.io.Serializable;

import fusionOptics.types.Pol;

import otherSupport.RandomManager;

import binaryMatrixFile.BinaryMatrixFile;

import algorithmrepository.Algorithms;
import algorithmrepository.exceptions.NotImplementedException;

import oneLiners.OneLiners;

import net.jafama.FastMath;

/** Simple Gaussian intensity PSF.
 * 
 */
public class GaussianPSF extends PointSpreadFunction {
	private static final long serialVersionUID = -766093963343259L;
	private int minStatsCount = 20;
	
	double meanX = 0, meanY = 0;
	double cXX = 0, cXY = 0, cYY = 0;
	
	/** Intensity incoherent sum over all polarisation components and all points*/
	double I0;
	
	/** Principal vectors and values */
	double ux,uy,ul, vx,vy,vl;
	
	public GaussianPSF() {
		this.minStatsCount = 20;
	}
	
	public GaussianPSF(int minStatsCount) {
		this.minStatsCount = minStatsCount;
	}
	
	public GaussianPSF(double I0, double meanX, double meanY, double cXX, double cXY, double cYY) {
		this.I0 = I0;
		this.meanX = meanX;
		this.meanY = meanY;
		this.cXX = cXX;
		this.cXY = cXY;
		this.cYY = cYY;
		this.meanM = new double[]{ 1,0,0,0,  0,1,0,0,  0,0,1,0,  0,0,0,1 };
		calcBasisVecs();
	}
	
	@Override
	public double[] getCharacterisationData() {
		double ret[] = new double[22];
		ret[0] = I0; ret[1] = meanX; ret[2] = meanY; 
		ret[3] = cXX; ret[4] = cXY; ret[5] = cYY;
		meanM = new double[16];
		System.arraycopy(meanM, 0, ret, 6, 16);
		return ret;		
	}
	
	@Override
	public void setCharacterisationData(double[] data) {
		I0 = data[0];
		meanX = data[1];
		meanY = data[2];
		cXX = data[3];
		cXY = data[4];
		cYY = data[5];
		if(data.length > 6)
			copyPols(data); //this is the slow bit here. memory allocing?
		else
			meanM = null;
		calcBasisVecs();
	}
	
	public void copyPols(double[] data){
		meanM = new double[16];
		System.arraycopy(data, 6, meanM, 0, 16);
	}

	@Override
	public boolean isEmpty() { return I0 <= 0; }
	
	@Override
	public void setPoints(double pos[][], double E[][][]) {
		int nPts = pos.length;
		int nPols = nPts <= 0 ? 0 : E[0].length;
		
		//per point intensity over all stokes vecs
		double I[] = new double[nPts];
		
		meanX = 0; meanY = 0; I0 = 0; 
		int nOK = 0;
		for(int i=0; i < nPts; i++){
			if(pos[i] != null){
				I[i] = 0;
				//stokes vectors already have the I in them, so just add incoherently
				for(int j=0; j < nPols; j++){ 
					I[i] += Pol.intensity(E[i][j]);
				}
				
				meanX += I[i] * pos[i][0];
				meanY += I[i] * pos[i][1];
				I0 += I[i];
				if(I[i] > 0)
					nOK++;
			}
		}
		
		if(I0 <= 0 || nOK < minStatsCount){
			invalidate();
			return;
		}
		
		meanX /= I0;
		meanY /= I0;
		
		cXX = 0; cXY = 0; cYY = 0;
		for(int i=0; i < nPts; i++){
			if(pos[i] != null){
				cXX += I[i] * (pos[i][0]-meanX)*(pos[i][0]-meanX);
				cXY += I[i] * (pos[i][0]-meanX)*(pos[i][1]-meanY);
				cYY += I[i] * (pos[i][1]-meanY)*(pos[i][1]-meanY);
			}
		}
		cXX /= I0;
		cXY /= I0;
		cYY /= I0;
		
		calcBasisVecs();
	}

	@Override
	public double[][] generatePoints(int nPoints) {
		double pos[][] = new double[2][nPoints];
			
		for(int i=0; i < nPoints; i++){
			double a = RandomManager.instance().nextNormal(0, 1);
			double b = RandomManager.instance().nextNormal(0, 1);
			pos[0][i] = meanX + a*ul*ux + b*vl*vx;
			pos[1][i] = meanY + a*ul*uy + b*vl*vy;
		}
		
		return pos;
	}
	
	@Override
	public void generatePolarisedPoints(int nPoints, double[][] pos,
			double[][][] E) {
		throw new NotImplementedException();
	}

	
	public void calcBasisVecs() {
		
		double ret[][] = Algorithms.eigenVecsAndVals2x2(new double[][]{{cXX,cXY*1}, {cXY,cYY}});

		ux = ret[0][0];
		uy = ret[0][1];
		ul = FastMath.sqrt(ret[2][0]);

		vx = ret[1][0];
		vy = ret[1][1];
		vl = FastMath.sqrt(ret[2][1]);
		
	}
	
	@Override
	public double[][] addToGrid(double[] x, double[] y, double G[][], double Imul) {
		return addToGridMonteCarlo(x, y, G, Imul, 40);
		//return addToGridMath(x, y, G, Imul);
	}

	public double[][] addToGridMonteCarlo(double[] x, double[] y, double G[][], double Imul, int nPoints) {
		double pos[][] = generatePoints(nPoints);
		for(int i=0; i < nPoints; i++){
			//int iiX = OneLiners.getNearestLowerIndex(x, pos[0][i]);
			//int iiY = OneLiners.getNearestLowerIndex(y, pos[1][i]);
			int iiX = (int)((pos[0][i] - x[0]) / (x[1] - x[0]));
			int iiY = (int)((pos[1][i] - y[0]) / (y[1] - y[0]));
			if(iiX >= 0 && iiY >= 0 && iiX < x.length && iiY < y.length)
				G[iiY][iiX] += Imul * I0 / nPoints;
		}
		return G;
	}
	
	public double[][] addToGridMath(double[] x, double[] y, double G[][], double Imul) {
		
		double rtDetC =  FastMath.sqrt(cXX * cYY - cXY * cXY);
		
		int nX = x.length, nY = y.length;
		for(int iX=0; iX < nX; iX++) {
			double x0 = (iX == 0) ? (x[0] - (x[1] - x[0])/2) : (x[iX] - (x[iX] - x[iX-1])/2); 
			double x1 = (iX == (nX-1)) ? (x[nX-1] + (x[nX-1] - x[nX-2])/2) : (x[iX] + (x[iX+1] - x[iX])/2); 
			
			for(int iY=0; iY < y.length; iY++) {
				double y0 = (iY == 0) ? (y[0] - (y[1] - y[0])/2) : (y[iY] - (y[iY] - y[iY-1])/2); 
				double y1 = (iY == (nY-1)) ? (y[nY-1] + (y[nY-1] - y[nY-2])/2) : (y[iY] + (y[iY+1] - y[iY])/2); 
				
				G[iY][iX] += Imul * I0 * averageQuadGauss(x0 - meanX, x1 - meanX, 
									y0 - meanY, y1 - meanY, 
									0.1) * (x1-x0) * (y1-y0) / (rtDetC * 2 * Math.PI);
			}	
		}
		return G;
	}
	
	/** Returns the average of the amplitude 1 (not normalised to 1) gaussian inside the given rectangle  */
	private double averageQuadGauss(
			double x0, double x1,
			double y0, double y1,
			double maxNormDX) {
		
		double u00 = x0 * ux + y0 * uy; 
		double v00 = x0 * vx + y0 * vy; 
		double u01 = x0 * ux + y1 * uy; 
		double v01 = x0 * vx + y1 * vy; 
		double u10 = x1 * ux + y0 * uy; 
		double v10 = x1 * vx + y0 * vy; 
		double u11 = x1 * ux + y1 * uy; 
		double v11 = x1 * vx + y1 * vy;
		
		double minU = Math.min(Math.min(u00, u01), Math.min(u10, u11));
		double maxU = Math.max(Math.max(u00, u01), Math.max(u10, u11));
		double minV = Math.min(Math.min(v00, v01), Math.min(v10, v11));
		double maxV = Math.max(Math.max(v00, v01), Math.max(v10, v11));
		
		if(minU > 6 || maxU < -6 || minV > 6 || maxV < -6)
			return 0;
		
			
		int nSplit = 0;
		
		if( (maxU - minU) > maxNormDX ) {
			nSplit = (int)(((maxU - minU) / maxNormDX) + 1);
			
		} else if( (maxV - minV) > maxNormDX ) {
			nSplit = (int)(((maxV - minV) / maxNormDX) + 1);
		}
		
		if(nSplit > 0) {
			double dX = (x1 - x0) / nSplit;
			double dY = (y1 - y0) / nSplit;
			double sum = 0;
			for(int iX = 0; iX < nSplit; iX++){
				for(int iY = 0; iY < nSplit; iY++){
						sum += averageQuadGauss(
									x0 + iX*dX, x0 + (iX+1)*dX,  
									y0 + iY*dY, y0 + (iY+1)*dY,
									maxNormDX) / nSplit / nSplit;
				}
			}
			return sum;
			
		} else{
			double mx = (x0+x1)/2, my = (y0 +y1)/2;
			double u = mx * ux + my * uy; 
			double v = mx * vx + my * vy;
			return FastMath.exp( -FastMath.pow2(u/ul)/2 -FastMath.pow2(v/vl)/2 );
			
			
		}
	}


	@Override
	public double getMinX() { return meanX - 6*Math.sqrt(cXX); }
	@Override
	public double getMaxX() { return meanX + 6*Math.sqrt(cXX); }
	@Override
	public double getMinY() { return meanY - 6*Math.sqrt(cYY); }
	@Override
	public double getMaxY() { return meanY + 6*Math.sqrt(cYY); }

	@Override
	public String toString() {
		return "GaussianPSF: mX="+meanX+", mY="+meanY+", sX="+FastMath.sqrt(cXX)+", sY="+FastMath.sqrt(cYY);
	}

	@Override
	public void combine(PointSpreadFunction[] psfs, double[] coeffs) {
		
		double cSum = 0;
		I0 = 0;
		meanX = 0; meanY = 0;
		cXX = 0; cXY = 0; cYY = 0;
		meanM = null;
		meanStartDir = new double[3];
		meanStartUp = new double[3];
		
		for(int i=0; i < psfs.length; i++) {
			GaussianPSF gpsf = (GaussianPSF)psfs[i];
			
			//if(gpsf == null)
				//continue;
			if(gpsf == null || gpsf.I0 <= 0){ //don't include 0 intensity (invalid) ones
				// If any of the points are missing or invalid, we can't 
				// continue, because all points in this volume will congregate
				// at the meanX,meanY of the PSFs at the grid edge
				//invalidate();
				//return;
				
				//try to generate this point by extrapolating some others,
				//Look in each of the 6 directions away from this point for two neighbours.
				//If found, extrapolate from them.
				//Take average of all extrapolated values found;
				
				//throw new RuntimeException("Interpolation from invalid point");
				continue;
			}
			
			double c = coeffs[i];
			cSum += c;
			
			I0 += gpsf.I0 * c;
			meanX += gpsf.meanX * c;
			meanY += gpsf.meanY * c;
			cXX += gpsf.cXX * c;
			cXY += gpsf.cXY * c;
			cYY += gpsf.cYY *c;
			
			if(psfs[i].meanM != null){
				if(meanM == null)
					meanM = new double[16];
				for(int j=0; j < 16; j++) {
					meanM[j] += gpsf.meanM[j] * c;
				}
			}
			
			if(psfs[i].meanStartDir != null){
				meanStartDir[0] += psfs[i].meanStartDir[0] * c;
				meanStartDir[1] += psfs[i].meanStartDir[1] * c;
				meanStartDir[2] += psfs[i].meanStartDir[2] * c;
				meanStartUp[0] += psfs[i].meanStartUp[0] * c;
				meanStartUp[1] += psfs[i].meanStartUp[1] * c;
				meanStartUp[2] += psfs[i].meanStartUp[2] * c;
			}
		}
		if(I0 <= 0) { // nothing
			invalidate();
			return;
		}
		
		I0 /= cSum;
		meanX /= cSum;
		meanY /= cSum;
		cXX /= cSum;
		cXY /= cSum;
		cYY /= cSum;
		
		if(meanM != null)
			for(int j=0; j < 16; j++)
				meanM[j] /= cSum;

		meanStartDir[0] /= cSum;
		meanStartDir[1] /= cSum;
		meanStartDir[2] /= cSum;
		meanStartUp[0] /= cSum;
		meanStartUp[1] /= cSum;
		meanStartUp[2] /= cSum;

		if(meanM != null){
			for(int j=0; j < 16; j++)
				meanM[j] /= cSum;
		}
		
		calcBasisVecs();
		
//		setSourceFrom(psfs, coeffs);
//		setMuellerFrom(psfs, coeffs);		
	}
	
	private void invalidate() {
		I0 = 0;
		meanX = Double.NaN; meanY = Double.NaN; cXX = Double.NaN; cXY = Double.NaN; cYY = Double.NaN;
		ux = Double.NaN; uy = Double.NaN; ul = Double.NaN;
		vx = Double.NaN; vy = Double.NaN; vl = Double.NaN;
		meanM = null;
	}

	@Override
	public void multiplyIntensity(double Imul) {
		I0 *= Imul;
	}
	@Override
	public double getIntensity() {
		return I0;
	}

	public double getMeanX() { return meanX; }
	public double getMeanY() { return meanY; }
	public double getFWHM(){ return 2*FastMath.sqrt(2*FastMath.log(2)) * FastMath.sqrt(cXX + cYY); }
}
