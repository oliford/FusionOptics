package ipp.neutralBeams;

import java.util.LinkedList;

import oneLiners.OneLiners;
import algorithmrepository.Algorithms;
import fusionOptics.Util;
import fusionOptics.interfaces.NullInterface;
import fusionOptics.surfaces.Cylinder;
import fusionOptics.types.Element;
import fusionOptics.types.Optic;

/** Beam geometry calculated from the AUG-type beam geometry angles
 * (Used for AUG and W7X heating beams (but not RuDIX))
 * */
public abstract class SimpleBeamGeometryAUGType extends SimpleBeamGeometry {
	
	protected double Px[][];
	protected double u[][];
	protected double accelVoltage[];
	
	/** calcVectors() version taking per-box arrays (as in minerva) */
	protected void calcVectors(double R_P0[][], double theta[][], double psi[][], double alpha[][], double beta[][], double Dx[][], double Z0[][]) {
		int nBeams = R_P0.length * R_P0[0].length;
		
		//first sort out the geometry
		Px = new double[nBeams][3];
		u = new double[nBeams][3];
		
		for(int iBx=0; iBx < R_P0.length; iBx++){
			for(int iS=0; iS < R_P0[iBx].length; iS++){
				calcVectorsSingleBeam(iBx*4+iS, R_P0[iBx][iS], theta[iBx][iS], psi[iBx][iS], alpha[iBx][iS], beta[iBx][iS], Dx[iBx][iS], Z0[iBx][iS]);
			}	
		}
	}
	
	/** Calculates the point near beam-crossing and unit vector for a given beam */
	protected void calcVectorsSingleBeam(int iB, double R_P0, double theta, double psi, double alpha, double beta, double Dx, double Z0) {
		
		if(Dx == 0)
			Dx = 1e-3; //erm, Dx=0 breaks the uVec calc and I don't want to think about this AGAIN, so just make it 1mm... sorry
		
		double A[] = new double[3];
		
		//A is the horizontal beam crossing - the position where the beams with different angles in the horizontal plane cross each other 
		A[0] = R_P0 * Math.cos(theta);
		A[1] = R_P0 * Math.sin(theta);
		A[2] = Z0 + Dx * Math.tan(-beta); //sign of beta is set so that numbers from the AUG NBI website match the AUG NBI beams
		
		//Px is the vertical beam crossing - where beams with different vertical angles cross
		Px[iB][0] = A[0] + Dx * Math.cos(theta - psi - alpha);
		Px[iB][1] = A[1] + Dx * Math.sin(theta - psi - alpha);
		Px[iB][2] = Z0;
		
		double l = Dx / Math.cos(-beta); //length between A and Px, by def.
		
		if(l == 0)
			throw new RuntimeException("Can't calculate beam geometry properly. iB = " + iB + ", beta = " + beta + ", dx=" + Dx + ".");
		
		u[iB][0] = (A[0] - Px[iB][0]) / l;
		u[iB][1] = (A[1] - Px[iB][1]) / l;
		u[iB][2] = (A[2] - Px[iB][2]) / l;		
		
	}

	@Override
	public double[][] startAll() { return Px; }
	@Override
	public double[][] uVecAll() { return u; }
	@Override
	public double[] getVoltageAll() { return accelVoltage; }
	
	public void setVoltages(double accelVoltage[]){ this.accelVoltage = accelVoltage; }

	private int boxToBeams[][] = new int[][]{ { 0,1,2,3 }, { 4,5,6,7 }};	
	@Override
	public double[] startBox(int boxIdx) { return avgVec(Px, boxToBeams[boxIdx]);	}
	@Override
	public double[] uVecBox(int boxIdx)  { return OneLiners.reNorm(avgVec(u, boxToBeams[boxIdx]));	}
}
