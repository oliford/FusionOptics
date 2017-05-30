package fusionOptics.optics;

import jafama.FastMath;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import fusionOptics.interfaces.IsoIsoStdFresnel;
import fusionOptics.surfaces.Aspheric;
import fusionOptics.surfaces.Disc;
import fusionOptics.surfaces.Dish;
import fusionOptics.types.Interface;
import fusionOptics.types.Material;
import fusionOptics.types.Medium;
import fusionOptics.types.Optic;
import fusionOptics.types.Surface;
import algorithmrepository.exceptions.NotImplementedException;

/** One spherical surface and one flat
 * 
 * Convex: +ve RadCurv or focalLength
 *        ----|
 *       /    |
 *      /     |
 *     /      |
 *    |       | ------> Normal
 *    |       |
 *     \      |
 *      \     |
 *       \    |
 *        ----|
 * 
 * 
 * Concave: -ve RadCurv or focalLength
 *     -------|
 *     \      |
 *      \     |
 *       \    |
 *        |   | ------> Normal
 *        |   |
 *       /    |
 *      /     |
 *     /      |
 *     -------|
 * 
 *  * 
 * @author A.Burckhart, oliford
 */
public class SimplePlanarAsphericLens extends Optic {

	private double centre[];
	private double radiusOfCurv;
	private Disc flatSurface;
	private Aspheric curvedSurface;
	private Medium lensMedium;
	private boolean concave; 
	
	public final static double focalLengthToRadCurv(double f, double n){
		return (n - 1.0) * f;
	}
		
	/** Creates a lens with the given focal length for the given wavelength at room temperature.
	 * Lens has the given edge thickness */
	public static SimplePlanarAsphericLens fromFocalLengthAndEdgeThickness(String name, double centreFlat[], double normal[], double radius, double focalLength, double edgeThickness, 
			double conicConst, double polyCoeffs[], Medium lensMedium, Interface iface, double wavelength) {
				
		double n = lensMedium.getMaterial().getRefractiveIndex(0, wavelength, 300);
		
		double R = focalLengthToRadCurv(focalLength, n);

		return fromRadiusOfCurvAndEdgeThickness(name, centreFlat,  normal, radius, R, edgeThickness, conicConst, polyCoeffs, lensMedium, iface);
	}
	
	/** Creates a lens with the given focal length for the given wavelength at room temperature.
	 * Lens has the given edge thickness */
	public static SimplePlanarAsphericLens fromFocalLengthAndCentreThickness(String name, double centreFlat[], double normal[], double radius, double focalLength, double centreThickness, 
							double conicConst, double polyCoeffs[], Medium lensMedium, Interface iface, double wavelength) {
		
		double n = lensMedium.getMaterial().getRefractiveIndex(0, wavelength, 300);
		
		double radiusOfCurvature = (n-1)*focalLength;
		
		return new SimplePlanarAsphericLens(name, centreFlat,  normal, radius, radiusOfCurvature, centreThickness, conicConst, polyCoeffs, lensMedium, iface);
	}
		
	/** Creates a lens with the given radius of curvature and the given central thickness */
	public static SimplePlanarAsphericLens fromRadiusOfCurvAndCentreThickness(String name, double centreFlat[], double normal[], double radius, double radiusOfCurvature, double centreThickness, 
			double conicConst, double polyCoeffs[], Medium lensMedium, Interface iface) {
		
		return new SimplePlanarAsphericLens(name, centreFlat,  normal, radius, radiusOfCurvature, centreThickness, conicConst, polyCoeffs, lensMedium, iface);
	}

	public static SimplePlanarAsphericLens fromRadiusOfCurvAndEdgeThickness(String name, double centreFlat[], double normal[], double radius, double radiusOfCurvature, double edgeThickness, 
			double conicConst, double polyCoeffs[], Medium lensMedium, Interface iface) {

		double centreThickness;
		if(radiusOfCurvature >= 0){
			centreThickness = edgeThickness + Dish.depth(radiusOfCurvature, radius);
		}else{
			centreThickness = edgeThickness - Dish.depth(-radiusOfCurvature, radius);
		}
				
		return new SimplePlanarAsphericLens(name, centreFlat, normal, radius, radiusOfCurvature, centreThickness, conicConst, polyCoeffs, lensMedium, iface);		
	}
	
	/** Private constructor, use the fromXXX() functions */
	protected SimplePlanarAsphericLens(String name, double centreFlat[], double normal[], double radius, double radiusOfCurvature, double centreThickness, 
								double conicConst, double polyCoeffs[], Medium lensMedium, Interface iface) {
		super(name);		
		this.lensMedium = lensMedium;
		
		double curvedCentre[] = new double[3];
		double reverseNormal[] = new double[3];
		centre = new double[3];
		
		//dishes have their normal facing toward the centre of their radius of curvature 
		for(int i=0;i<3;i++){
			curvedCentre[i] = centreFlat[i] - centreThickness * normal[i];
			reverseNormal[i] = -normal[i];
			centre[i] = centreFlat[i] - (centreThickness/2.0) * normal[i];
		}
		
		flatSurface = new Disc(getName() + "-flat", centreFlat, normal, radius, null, lensMedium, iface) ;
		
		if(radiusOfCurvature >= 0){ //convex
			this.concave = false;
			this.radiusOfCurv = radiusOfCurvature;
			curvedSurface = new Aspheric(getName() + "-curved", curvedCentre, normal, radiusOfCurv, radius, conicConst, polyCoeffs, lensMedium, null, iface);
			
		}else{ //concave
			this.concave = true;
			this.radiusOfCurv = -radiusOfCurvature;
			curvedSurface = new Aspheric(getName() + "-curved", curvedCentre, reverseNormal, radiusOfCurv, radius, conicConst, polyCoeffs, null, lensMedium, iface);
		}
		
		addElement(flatSurface);
		addElement(curvedSurface);
				
	}
		
	public static double thickLensFocalLength(double radius1, double radius2, double thickness, double n) {
		return 1.0 / ((n-1.0)*( 1/radius1 - 1/radius2 + (n-1)*thickness/(n*radius1*radius2) ));
	}

	public void setFocalLength(double f, double wavelength) {
		
		double n= lensMedium.getMaterial().getRefractiveIndex(0, wavelength, 300);
		double lensRadius = curvedSurface.getDishDiameter()/2;
		
		radiusOfCurv = focalLengthToRadCurv(f, n);
		
		//frontSurface.setRadiusOfCurv(radiusOfCurv);
		curvedSurface.setRadiusOfCurv(radiusOfCurv);
		
	}

	public Aspheric getAsphericSurface() { return curvedSurface;	}
	public Disc getPlanarSurface() { return flatSurface;	}
	
	public Dish getBackSurface() { return curvedSurface;	}
	public Disc getFrontSurface() { return flatSurface;	}

	public double getRadius() { return Math.max(flatSurface.getRadius(), curvedSurface.getDishDiameter()/2); }
	
}
