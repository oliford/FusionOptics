package fusionOptics.pointSpread;

import java.io.Serializable;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;

import fusionOptics.types.Pol;
import seed.matrix.Mat;
import uk.co.oliford.jolu.OneLiners;

/** General base for point spread functions.
 * 
 * Handles averaging of source 'up' sense and emission direction for the PSF.
 * Also provides the general formulation of Muller matrix, by incoherent addition
 * of all the collected rays from the range of standard set of initial polarisations.
 * 
 * @author oliford
 *
 */
public abstract class PointSpreadFunction implements Serializable {
	private static final long serialVersionUID = 670176178804905200L;
	final static double rt2 = Math.sqrt(2);
	
	protected double meanStartDir[];
	protected double meanStartUp[];
	
	/** Mean Mueller matrix for effect on polarisation, in terms of the mean 'Up' source 
	 * frame, to the 'up' sense of the final polarisation states passed to calcAverageMuellerMatrix 
	 * (probably from the PSF collector's polarisation plane).
	 * 
	 * This is a flat 4x4 matrix.
	 */
	protected double meanM[];
	
	/** The 4 polarisations states (EuRe, EuIm, ErRe, ErIm) to start each ray with, that will
	 *  allow calculation of Mullear matricies from the final polarisations states.
	 *  
	 *  This was a static final field but I keep forgetting to clone() it and messing it up
	 *  so it is now a method and is created each time.
	 */
	public final static double[][] getInputStatesForMuellerCalc(){
		return new double[][]{
			{ 1, 0, 0, 0 }, // linearly polarised parallel with 'up' direction, phi=0
			{ 0, 0, 1, 0 }, // linearly polarised parallel with 'right' direction, phi=0
			{rt2/2,0,rt2/2,0}, // linearly polarised at 45' between 'up' and 'right'
			{rt2/2,0,0,-rt2/2}  //circularly 'right' polarised -
						//clockwise looking forward, moving forward, starting 'up' and going towards 'right'		
		};
	}
	
	/** The 4 stokes vectors to match inputStatesForMuellerCalc */
	public final static double inputStokesForMuellerCalc[][] = new double[][]{
		{ 1,  1, 0, 0 }, // linearly polarised parallel with 'up' direction, (phases lost)
		{ 1, -1, 0, 0 }, // linearly polarised parallel with 'right' direction, (phases lost)
		{ 1,  0, 1, 0 }, // linearly polarised at 45' between 'up' and 'right'
		{ 1,  0, 0, 1 }  //circularly 'right' polarised (phases lost)	
	};
	
	/** The inverse of inputStokesForMuellerCalc[][], used to calculate M from the output stokes vectors */
	private final static double invInputStokes[][] = new double[][]{
	   {  0.5,  0.5, 0, 0 },
	   {  0.5, -0.5, 0, 0 },
	   { -0.5, -0.5, 1, 0 },
	   { -0.5, -0.5, 0, 1 }
	};
	
	/** Calculates the average Mueller matrix, assuming all light in the PSF is merged into one,
	 *   or has the same behaviour and that the inital polarisation states of the ray 
	 *   was inputStatesForMuellerCalc[][]
	 * 
	 *  The ray tracer itself can't describe 'unpolarised' light without actual randomisation of the phases
	 *  in the pol state arrays but we can still do this by looking at the response to each of the fully polarised
	 *  light signals and solving the Mueller matrix from that.
	 *  
	 *  We know that the output stokes vector will always be sOut = M . sIn
	 *  so M = sOut . inv(sIn)
	 */
	public void calcAverageMuellerMatrix(double E[][][]){
		double meanSOut[][] = new double[4][4];

		//stokes vectors already have the I in them, so just add incoherently
		int n=0;
		for(int i=0; i < E.length; i++){
			if(E[i] != null){
				for(int j=0; j < 4; j++){ 
					meanSOut[j][0] += Pol.s0(E[i][j]);
					meanSOut[j][1] += Pol.s1(E[i][j]);
					meanSOut[j][2] += Pol.s2(E[i][j]);
					meanSOut[j][3] += Pol.s3(E[i][j]);
				}
				n++;
			}
		}
		
		for(int j=0; j < 4; j++){
			meanSOut[j][0] /= n;
			meanSOut[j][1] /= n;
			meanSOut[j][2] /= n;
			meanSOut[j][3] /= n;
		}
				
		meanM = OneLiners.flatten(Mat.mul(invInputStokes, meanSOut));
	}
	
	public double[] getMeanMueller(){ return meanM; }
	
	public double[] getNormalisedMeanMueller(){
		double m[] = meanM.clone();
		for(int i=0; i < m.length; i++)
			m[i] /= meanM[0];
		return m;
	}
	
	/** Sets the source point information for the PSF
	 * 
	 * @param dirs Starting ray directions to take average of
	 * @param ups Starting ray 'up' definitions, to take average of
	 * @param I Intensity
	 */
	public void setSourceInfo(double dirs[][], double ups[][], double E[][][]){
		meanStartDir = new double[3];
		meanStartUp = new double[3];
		
		double Isum = 0;
		for(int i=0; i < dirs.length; i++){
			double I = Pol.intensity(E[i]);
			meanStartDir[0] += I * dirs[i][0];
			meanStartDir[1] += I * dirs[i][1];
			meanStartDir[2] += I * dirs[i][2];
			meanStartUp[0] += I * ups[i][0];
			meanStartUp[1] += I * ups[i][1];
			meanStartUp[2] += I * ups[i][2];
			Isum += I;
		}
		meanStartDir[0] /= Isum;
		meanStartDir[1] /= Isum;
		meanStartDir[2] /= Isum;
		meanStartUp[0] /= Isum;
		meanStartUp[1] /= Isum;
		meanStartUp[2] /= Isum;		
	}
	
	public double[] getMeanSourceDir(){ return meanStartDir; }
	public double[] getMeanSourceUp(){ return meanStartUp; }
	public void setMeanSourceDir(double meanStartDir[]){ this.meanStartDir = meanStartDir; }
	public void setMeanSourceUp(double meanStartUp[]){ this.meanStartUp = meanStartUp; }
	
	public abstract boolean isEmpty();
	
	/** Sets the PSF's internal info (stats etc) from the given points
	 * 
	 * E is the polarisation state {EuRe,EuIm,ErRe,ErIm} from each of the initial
	 * polarisations. Usually, these should have been: {1,0,0,0},{0,0,1,0} so
	 * that we get a sense of what happens to 'up' and what happens to 'right'
	 */
	public abstract void setPoints(double pos[][], double E[][][]);
	
	
	/** Generate random samples from the PSF
	 * 
	 * @return double[x/y][pointIdx]
	 */
	public abstract double[][] generatePoints(int nPoints);
	
	/** Fills the given array with positions and polarisations 
	 * generated randomly according to the PSF.
	 * 
	 * @return double[{x/y},{pols}][pointIdx]
	 */
	public abstract void generatePolarisedPoints(int nPoints, double pos[][], double E[][][]);
	

	public final double[][] getGridProbability(double[] x, double[] y) {
		double P[][] = new double[y.length][x.length];
		return addToGrid(x, y, P, 1);
	}
	
	/** Adds the PSF's probability in the given grid cells, to the given grid */ 
	public abstract double[][] addToGrid(double[] x, double[] y, double G[][], double I);
	
	public abstract double[][] addToGridMonteCarlo(double[] x, double[] y, double G[][], double I, int nPoints);
	
	public abstract double getMinX();
	public abstract double getMaxX();
	public abstract double getMinY();
	public abstract double getMaxY();
	
	/** Return characterisation data for the PSF (fixed length) */
	public abstract double[] getCharacterisationData();

	/** Init the PSF from its characterisation data (fixed length) */
	public abstract void setCharacterisationData(double data[]);
	
	/** Sets this PSF from weights average of the given PSFs.
	 * The base class implementation assumes all data can be averaged */
	public void combine(PointSpreadFunction psfs[], double coeffs[]){
		double avgData[] = null;
		
		for(int i=0; i < psfs.length; i++){
			double data[] = psfs[i].getCharacterisationData();
			if(avgData == null)
				avgData = new double[data.length];
			for(int j=0; j < data.length; j++)
				avgData[j] += data[j] * coeffs[i];
		}

		setCharacterisationData(avgData);
		setSourceFrom(psfs, coeffs);
		
	}

	protected void setSourceFrom(PointSpreadFunction psfs[], double coeffs[]){
		meanStartDir = new double[3];
		meanStartUp = new double[3];
		
		for(int i=0; i < psfs.length; i++){
			if(coeffs[i] == 0)
				continue;
			meanStartDir[0] += psfs[i].meanStartDir[0] * coeffs[i];
			meanStartDir[1] += psfs[i].meanStartDir[1] * coeffs[i];
			meanStartDir[2] += psfs[i].meanStartDir[2] * coeffs[i];
			meanStartUp[0] += psfs[i].meanStartUp[0] * coeffs[i];
			meanStartUp[1] += psfs[i].meanStartUp[1] * coeffs[i];
			meanStartUp[2] += psfs[i].meanStartUp[2] * coeffs[i];
		}
	}
	
	protected void setMuellerFrom(PointSpreadFunction psfs[], double coeffs[]){
		meanM = new double[16];
		
		for(int i=0; i < psfs.length; i++){
			if(coeffs[i] == 0)
				continue;
			for(int j=0; j < 16; j++)
				meanM[j] += psfs[i].meanM[j] * coeffs[i];
		}
	}
	
	public void dumpMueller() { dumpMueller(getMeanMueller());	}
	
	public void dumpNormalisedMueller() { dumpMueller(getNormalisedMeanMueller());	}
	
	public static void dumpMueller(double m[]) {
		
		DecimalFormat fmt = new DecimalFormat("0.###");
		DecimalFormatSymbols symbs = fmt.getDecimalFormatSymbols();
		symbs.setNaN("NaN");
		symbs.setInfinity("Inf");
		fmt.setDecimalFormatSymbols(symbs);
		
		for(int i=0; i < 4; i++){
			for(int j=0; j < 4; j++)
				System.out.print(fmt.format(m[i*4+j]) + ",\t");
			System.out.println();
		}
	}
	
	/** Multiply the PSF's intensity */
	public abstract void multiplyIntensity(double Imul);

	public abstract double getIntensity();
}
