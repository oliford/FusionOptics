package fusionOptics.optics;

import net.jafama.FastMath;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import fusionOptics.Util;
import fusionOptics.interfaces.IsoIsoStdFresnel;
import fusionOptics.surfaces.Aspheric;
import fusionOptics.surfaces.Cylinder;
import fusionOptics.surfaces.Disc;
import fusionOptics.surfaces.Dish;
import fusionOptics.surfaces.TruncatedCylinder;
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
 *     \      |    ^
 *      \     |    | up
 *       \    |    |
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
 *       /    |    ^
 *      /     |    | up
 *     /      |    |
 *     -------|
 *     
 *     
 *     up is the axis along which the surface is curved
 *     the surface is straight along: right = normal x up
 * 
 *  */
public class CylindricalPlanarConvexLens extends Optic {

	private double centre[];
	private double radiusOfCurv;
	private Disc flatSurface;
	private TruncatedCylinder curvedSurface;
	private Medium lensMedium;
	private boolean concave;
	
	public final static double focalLengthToRadCurv(double f, double n){
		return (n - 1.0) * f;
	}
		
	/** Creates a lens with the given focal length for the given wavelength at room temperature.
	 * Lens has the given edge thickness */
	public static CylindricalPlanarConvexLens fromFocalLengthAndEdgeThickness(String name, double centreFlat[], double normal[], double up[], double radius, double focalLength, double edgeThickness, double length, Medium lensMedium, Interface iface, double wavelength) {
				
		double n = lensMedium.getMaterial().getRefractiveIndex(0, wavelength, 300);
		
		double R = focalLengthToRadCurv(focalLength, n);

		return fromRadiusOfCurvAndEdgeThickness(name, centreFlat,  normal, up, radius, R, edgeThickness, length, lensMedium, iface);
	}
	
	/** Creates a lens with the given focal length for the given wavelength at room temperature.
	 * Lens has the given edge thickness */
	public static CylindricalPlanarConvexLens fromFocalLengthAndCentreThickness(String name, double centreFlat[], double normal[], double up[], double radius, double focalLength, double centreThickness, double length, Medium lensMedium, Interface iface, double wavelength) {
		
		double n = lensMedium.getMaterial().getRefractiveIndex(0, wavelength, 300);
		
		double R = (n-1)*focalLength;
		
		return new CylindricalPlanarConvexLens(name, centreFlat,  normal, up, radius, R, centreThickness, length, lensMedium, iface);
	}
		
	/** Creates a lens with the given radius of curvature and the given central thickness */
	public static CylindricalPlanarConvexLens fromRadiusOfCurvAndCentreThickness(String name, double centreFlat[], double normal[], double up[], double radius, double radiusOfCurvature, double centreThickness, double length, Medium lensMedium, Interface iface) {
		
		return new CylindricalPlanarConvexLens(name, centreFlat,  normal, up, radius, radiusOfCurvature, centreThickness, length, lensMedium, iface);
	}

	public static CylindricalPlanarConvexLens fromRadiusOfCurvAndEdgeThickness(String name, double centreFlat[], double normal[], double up[], double radius, double radiusOfCurvature, double edgeThickness, double length, Medium lensMedium, Interface iface) {
		
		double centreThickness;
		if(radiusOfCurvature >= 0){
			centreThickness = edgeThickness + Dish.depth(radiusOfCurvature, radius);
		}else{
			centreThickness = edgeThickness - Dish.depth(-radiusOfCurvature, radius);
		}
		return new CylindricalPlanarConvexLens(name, centreFlat, normal, up, radius, radiusOfCurvature, centreThickness, length, lensMedium, iface);
		
	}
	
	/** Private constructor, use the fromXXX() functions */
	private CylindricalPlanarConvexLens(String name, double centreFlat[], double normal[], double up[], double radius, double radiusOfCurvature, double centreThickness, double length, Medium lensMedium, Interface iface) {
		super(name);
		this.lensMedium = lensMedium;

		
		double curvedCentre[] = new double[3];
		double cyldCentre[] = new double[3];
		double reverseNormal[] = new double[3];
		centre = new double[3];
		
		//dishes have their normal facing toward the centre of their radius of curvature 
		for(int i=0;i<3;i++){
			curvedCentre[i] = centreFlat[i] - centreThickness * normal[i];
			cyldCentre[i] = curvedCentre[i] + radiusOfCurvature * normal[i];
			reverseNormal[i] = -normal[i];
			centre[i] = centreFlat[i] - (centreThickness/2.0) * normal[i];
		}
		
		flatSurface = new Disc(getName() + "-flat", centreFlat, normal, radius, null, lensMedium, iface) ;
		
		if(radiusOfCurvature >= 0){ //convex
			this.concave = false;
			this.radiusOfCurv = radiusOfCurvature;
			double cyldAxis[] = Util.reNorm(Util.cross(normal, up));
			curvedSurface = new TruncatedCylinder(getName() + "-curved", cyldCentre, cyldAxis, up, radiusOfCurv, length, 2*radius, lensMedium, null, iface);
			
		}else{ //concave
			this.concave = true;
			this.radiusOfCurv = -radiusOfCurvature;
			double cyldAxis[] = Util.reNorm(Util.cross(reverseNormal, up));
			curvedSurface = new TruncatedCylinder(getName() + "-curved", cyldCentre, cyldAxis, up, radiusOfCurv, length, 2*radius, null, lensMedium, iface);
		}

		addElement(flatSurface);
		addElement(curvedSurface);
	}
		
	public static double thickLensFocalLength(double radius1, double radius2, double thickness, double n) {
		return 1.0 / ((n-1.0)*( 1/radius1 - 1/radius2 + (n-1)*thickness/(n*radius1*radius2) ));
	}

	public void setFocalLength(double f, double wavelength) {
		
		double n= lensMedium.getMaterial().getRefractiveIndex(0, wavelength, 300);
		// double lensRadius = curvedSurface.getDishDiameter()/2;
		
		radiusOfCurv = focalLengthToRadCurv(f, n);
		
		//frontSurface.setRadiusOfCurv(radiusOfCurv);
		// curvedSurface.setRadiusOfCurv(radiusOfCurv);
		throw new RuntimeException("Not implemented");
		
	}

	public TruncatedCylinder getConvexSurface() { return curvedSurface;	}
	public Disc getPlanarSurface() { return flatSurface;	}
	
	public TruncatedCylinder getBackSurface() { return curvedSurface;	}
	public Disc getFrontSurface() { return flatSurface;	}

	public double getRadius() { 
		//return Math.max(flatSurface.getRadius(), curvedSurface.getDishDiameter()/2);
		throw new RuntimeException("Not implemented");
	}
	
}
