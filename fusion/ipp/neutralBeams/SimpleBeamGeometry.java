package ipp.neutralBeams;

import java.util.LinkedList;

import algorithmrepository.Algorithms;
import fusionOptics.Util;
import fusionOptics.interfaces.NullInterface;
import fusionOptics.surfaces.Cylinder;
import fusionOptics.types.Element;
import fusionOptics.types.Optic;

/** Definition of neutral beam geometry, drawing etc 
 * For now just the fixed coords from minerva-land but later we'll
 * bring that geom stuff here
 * 
 * */
public abstract class SimpleBeamGeometry  {

	/** The beam indices actually match the AUG naming convention (except for the 1-based indexing)
	 * Q means 'Quelle' ('source' in German) */
	public static final int BEAM_BOXAVG = -1;
	public static final int BEAM_Q1 = 0;
	public static final int BEAM_Q2 = 1;
	public static final int BEAM_Q3 = 2;	
	public static final int BEAM_Q4 = 3;
	public static final int BEAM_Q5 = 4;
	public static final int BEAM_Q6 = 5;
	public static final int BEAM_Q7 = 6;
	public static final int BEAM_Q8 = 7;
	
		
	public final double[] start(int beamIdx){ return startAll()[beamIdx]; }
	public final double[] uVec(int beamIdx){ return uVecAll()[beamIdx];  }
	public final int nBeams(){ return startAll().length; }
	
	public abstract double[][] startAll();
	public abstract double[][] uVecAll();
	public abstract double[] getVoltageAll();
	
	public abstract double[] startBox(int boxIdx);
	public abstract double[] uVecBox(int boxIdx);
	
	protected double beamWidth;
	protected double plasmaR0;
	protected double plasmaR1;
	protected double sourceR;
	
	public double beamWidth(){ return beamWidth; }
	public double plasmaR0(){ return plasmaR0; }
	public double plasmaR1(){ return plasmaR1; }
	public double sourceR(){ return sourceR; }
		
	public double getLOfBeamAxisAtR(int beamIdx, double R) {
		double t[] = Algorithms.cylinderLineIntersection(start(beamIdx), uVec(beamIdx), new double[]{0,0,0}, new double[]{0,0,1}, R*R);
		return (t.length == 0) ? Double.NaN : Math.min(t[0],t[1]);
	}
	
	public double[] getPosOfBeamAxisAtR(int beamIdx, double R){		
		return Util.plus(start(beamIdx), Util.mul(uVec(beamIdx), getLOfBeamAxisAtR(beamIdx, R)));
	}
	
	public double[] getPosOfBoxAxisAtR(int boxIdx, double R){		
		double t[] = Algorithms.cylinderLineIntersection(startBox(boxIdx), uVecBox(boxIdx), new double[]{0,0,0}, new double[]{0,0,1}, R*R);
		return Util.plus(startBox(boxIdx), Util.mul(uVecBox(boxIdx), Math.min(t[0],t[1])));
	}
	
	public Optic makeAllBeamCylds(){
		return makeAllBeamCylds(0.05, 0.16);
	}
	
	public Optic makeAllBeamCylds(double dL, double fwhm){
		double l0 = 0.2;
		double l1 = 2.9;		
		return makeAllBeamCylds(dL, fwhm, l0, l1);
	}
		
	public Optic makeAllBeamCylds(double dL, double fwhm, double l0, double l1){
		LinkedList<Element> beams = new LinkedList<Element>();
		
		for(int i=0; i < nBeams(); i++){
			Optic beam = makeBeamCylds(i, dL, fwhm, l0, l1);
			beams.add(beam);
		}
		
		return new Optic("beamCylds", beams);
	}
	
	public Optic makeBeamCylds() {
		return makeBeamCylds(1, 0.10, 0.16);
	}
	
	public Optic makeBeamCylds(int iB, double dL, double fwhm) {
		double l0 = 0.2;
		double l1 = 2.9;		
		return makeBeamCylds(iB, dL, fwhm, l0, l1);
	}
		
	public Optic makeBeamCylds(int iB, double dL, double fwhm, double l0, double l1) {
		LinkedList<Element> cylds = new LinkedList<Element>();
		
		
		for(int i=0; i < (l1 - l0) / dL; i++){
			double l = l0 + i * dL;			
			Cylinder clyd = new Cylinder(
					"cyld"+i,
					new double[]{ 
							start(iB)[0] + l * uVec(iB)[0],
							start(iB)[1] + l * uVec(iB)[1],
							start(iB)[2] + l * uVec(iB)[2],							
					},
					uVec(iB),
					fwhm/2, //radius, HWHM of beam
					dL, //length 
					NullInterface.ideal());
			clyd.setDrawingDetails(8, 10);
			
			cylds.add(clyd);
		}
		
		return new Optic("beamCylds" + iB, cylds);
	}
	
	protected final double[] avgVec(double vecs[][], int idxs[]){
		double vec[] = new double[3];
		for(int i=0; i < idxs.length; i++)
			for(int k=0; k < 3; k++)
				vec[k] += vecs[idxs[i]][k] / idxs.length;
		return vec;
	}
	
}
