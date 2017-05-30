package fusionOptics.pointSpread;

import java.io.Serializable;

import fusionOptics.types.Pol;

import otherSupport.RandomManager;

import binaryMatrixFile.BinaryMatrixFile;

import oneLiners.OneLiners;
import algorithmrepository.Algorithms;
import algorithmrepository.Interpolation1D;
import algorithmrepository.Interpolation2D;
import algorithmrepository.LinearInterpolation1D;
import algorithmrepository.LinearInterpolation2D;
import algorithmrepository.exceptions.NotImplementedException;

/** Represents the PSF directly with the originally collected points */
public class MiniImagePSF extends PointSpreadFunction {
	private static final long serialVersionUID = -6188443371317171763L;

	private int nX, nY;
	private double minX, maxX, minY, maxY;
	private double dX, dY;
	
	private double pdf[][];
	
	LinearInterpolation2D imageInterp;
	
	LinearInterpolation1D colSumCDF;
	LinearInterpolation1D colCDFs[];
	
	public MiniImagePSF(int nX, int nY) {
		this.nX = nX;
		this.nY = nY;
		this.minX = Double.NaN;
	}
	
	@Override
	public void setPoints(double[][] pos, double[][][] E) {
		if(pos.length <= 2){
			minX = Double.NaN;
			pdf = null;
			return;
		}

		double ret[][] = OneLiners.columnRanges(pos);
		minX = ret[0][0]; maxX = ret[1][0];
		minY = ret[0][1]; maxY = ret[1][1];
		
		dX = (maxX - minX) / (nX-2);
		dY = (maxY - minY) / (nY-2);
		
		//build the interpolatable PDF, which should end up with 0 all around the border
		double x[] = OneLiners.linSpace(minX-dX/2, maxX+dX/2, dX); 
		double y[] = OneLiners.linSpace(minY-dY/2, maxY+dY/2, dY); 
		pdf = new double[nX][nY];
		
		double sum = 0;
		for(int i=0; i < pos.length; i++) {
			double I = 0;
			for(int j=0; j < E[i].length; j++)
				I += Pol.intensity(E[i][j]); 
			
			int iX = (int)((pos[i][0] - minX - dX*1e-5) / dX) + 1;
			int iY = (int)((pos[i][1] - minY - dY*1e-5) / dY) + 1;
			
			//if(iX > 0 && iY > 0 && iX < nX && iY < nY){
				pdf[iX][iY] += I;
				sum += I;
			//}
		}
		
		colSumCDF = null;
		imageInterp = new LinearInterpolation2D(x, y, pdf, 0.0);
		
		
	}
	
	@Override
	public boolean isEmpty() { return pdf == null; }
	
	private void calcCDFs(){
		colCDFs = new LinearInterpolation1D[nX];
		
		//the CDFs are build around the grid cell concept of the P[][]
		//not the intepolatable PDF.
		double x[] = OneLiners.linSpace(minX, maxX, nX-1);
		double y[] = OneLiners.linSpace(minY, maxY, nY-1);
		
		double csCDF[] = new double[nX-1];
		
		for(int iX=0; iX < (nX-1); iX++){
			double cdf[] = new double[nY-1];
			
			for(int iY=0; iY < (nY-1); iY++){
				cdf[iY] = ((iY == 0) ? 0 : cdf[iY-1]) + pdf[iX][iY];
			}
			double colSum = cdf[nY-2]; //total for the column
				
			if(colSum > 0){
				for(int iY=0; iY < (nY-1); iY++){
					cdf[iY] /= colSum;
				}
				colCDFs[iX] = new LinearInterpolation1D(cdf, y);
			}else{
				colCDFs[iX] = null; //no probability at all
			}
			
			csCDF[iX] = ((iX == 0) ? 0 : csCDF[iX-1]) + colSum;
		}
		
		for(int iX=0; iX < (nX-1); iX++){
			csCDF[iX] /= csCDF[nX-2];
		}	
		colSumCDF = new LinearInterpolation1D(csCDF, x);
	}

	@Override
	public double[][] generatePoints(int nPoints) {
		if(Double.isNaN(minX))
			return null;
		if(colSumCDF == null)
			calcCDFs();
		
		double ret[][] = new double[2][nPoints];
		for(int i=0; i < nPoints; i++){
			double a = RandomManager.instance().nextUniform(0, 1);
			double b = RandomManager.instance().nextUniform(0, 1);
			ret[0][i] = colSumCDF.eval(a);
			int iX = (int)((ret[0][i] - minX) / dX) + 1;
			ret[1][i] = colCDFs[iX].eval(b);
		}
		
		return ret;
	}

	@Override
	public void generatePolarisedPoints(int nPoints, double[][] pos,
			double[][][] E) {
		throw new NotImplementedException();
	}
	
	public double[][] addToGrid(double[] x, double[] y, double[][] G, double I0) {
		return addToGridMonteCarlo(x, y, G, I0, 1000);
	}
	
	@Override
	public double[][] addToGridMonteCarlo(double[] x, double[] y, double[][] G, double I0, int nPoints) {
		if(Double.isNaN(minX))return G;
		
		/* This is a nasty image rescaling integration/interpolation problem
		 * 
		 * Instead, we're just going to draw points from the image PDF
		 * and put them into the main image
		 * 
		 */
		
		double pos[][] = generatePoints(nPoints);
		for(int i=0; i < nPoints; i++){
			int iiX = OneLiners.getNearestLowerIndexProbablyRegular(x, pos[0][i]);
			int iiY = OneLiners.getNearestLowerIndexProbablyRegular(y, pos[1][i]);
			G[iiY][iiX] += I0 / nPoints;
		}
		
		/*
		 * Method 1: Eval mini image at centre of each target image pixel  
		int iX0 = OneLiners.getNearestLowerIndex(x, minX);
		int iX1 = OneLiners.getNearestLowerIndex(x, maxX)+1;
		int iY0 = OneLiners.getNearestLowerIndex(y, minY);
		int iY1 = OneLiners.getNearestLowerIndex(y, maxY)+1;
		
		for(int iX = iX0; iX < iX1; iX++){
			for(int iY = iY0; iY < iY1; iY++){
				G[iX][iY] += I0 * imageInterp.eval(x[iX], y[iY]);
			}	
		}*/
		

		/* Method 2: Add value of each mini image pixel to the pixel of the target
		 * image which contains the that mini image pixel's centre.   
		for(int iX = 0; iX < nX; iX++) {
			for(int iY = 0; iY < nY; iY++) {
				int iiX = OneLiners.getNearestLowerIndex(x, minX + iX*dX + dX/2);
				int iiY = OneLiners.getNearestLowerIndex(y, minY + iY*dY + dY/2);
				G[iiX][iiY] += I0 * pdf[iX][iY];
			}
		}
		*/
				
				
		return G;
	}

	@Override
	public double[] getCharacterisationData() {
		double ret[] = new double[6+nX*nY];
		ret[0] = nX;
		ret[1] = nY;
		if(Double.isNaN(minX)){
			ret[2] = 0;
			ret[3] = 0;
			ret[4] = 0;
			ret[5] = 0;
		} else {
			ret[2] = minX;
			ret[3] = minY;
			ret[4] = maxX;
			ret[5] = maxY;
			for(int iX=0; iX < nX; iX++) {
				System.arraycopy(pdf[iX], 0, ret, 6+iX*nY, nY);				
			}
		}
		return ret;
	}

	@Override
	public void setCharacterisationData(double[] data) {
		throw new NotImplementedException();
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
		for(int iX = 0; iX < nX; iX++)
			for(int iY = 0; iY < nY; iY++)
				pdf[iX][iY] *= Imul;
	}

	@Override
	public double getIntensity() {
		double sum = 0;
		for(int iX = 0; iX < nX; iX++)
			for(int iY = 0; iY < nY; iY++)
				sum += pdf[iX][iY];
		return sum;
	}
}
