package fusionOptics.optics;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import fusionOptics.interfaces.IsoIsoStdFresnel;
import fusionOptics.surfaces.Disc;
import fusionOptics.surfaces.Dish;
import fusionOptics.types.Interface;
import fusionOptics.types.Material;
import fusionOptics.types.Medium;
import fusionOptics.types.Optic;
import fusionOptics.types.Surface;

import algorithmrepository.exceptions.NotImplementedException;

/** Two equally curved spherical surfaces 
 * 
 * @author oliford
 */
public class SimpleDoubleConvexLens extends Optic {

	private double centre[];
	private double radiusOfCurv;
	private Dish frontSurface, backSurface;
	private Medium lensMedium;
	
	public final static double focalLengthToRadCurv(double f, double lensRadius, double n){
		double R = 2*(n-1)*f;
		for(int i=0 ; i < 10; i++)
			R = 1 / ( 1/(2*(n-1)*f) - (n-1)*lensRadius*lensRadius/(2*n*R*R*R));
		return R;
	}
		
	/** Creates a lens with the given focal length and edge thickness */
	public static SimpleDoubleConvexLens fromFocalLengthAndEdgeThickness(String name, double centre[], double normal[], double radius, double focalLength, double edgeThickness, Medium lensMedium, Interface iface, double wavelength) {
		
		double n = lensMedium.getMaterial().getRefractiveIndex(0, wavelength, 300);
		
		double R = focalLengthToRadCurv(focalLength, radius, n);

		return fromRadiusOfCurvAndEdgeThickness(name, centre,  normal, radius, R, edgeThickness, lensMedium, iface);
	}
	
	/** Creates a lens with the given focal length and centre thickness */
	public static SimpleDoubleConvexLens fromFocalLengthAndCentreThickness(String name, double centre[], double normal[], double radius, double focalLength, double centreThickness, Medium lensMedium, Interface iface, double wavelength) {
		
		double n = lensMedium.getMaterial().getRefractiveIndex(0, wavelength, 300);
		
		double R = focalLengthToRadCurv(focalLength, radius, n);

		return fromRadiusOfCurvAndCentreThickness(name, centre,  normal, radius, R, centreThickness, lensMedium, iface);
	}
	
	public static SimpleDoubleConvexLens fromRadiusOfCurvAndEdgeThickness(String name, double centre[], double normal[], double radius, double radiusOfCurvature, double edgeThickness, Medium lensMedium, Interface iface) {
		
		double centreThickness;
		if(radiusOfCurvature > 0){
			double dishDepth = radiusOfCurvature - Math.sqrt(radiusOfCurvature*radiusOfCurvature - radius*radius);
			centreThickness = edgeThickness + 2 * dishDepth;
		}else{
			radiusOfCurvature = - radiusOfCurvature;
			double dishDepth = radiusOfCurvature - Math.sqrt(radiusOfCurvature*radiusOfCurvature - radius*radius);
			centreThickness = -(edgeThickness - 2*dishDepth);
		}
		return new SimpleDoubleConvexLens(name, centre, normal, radius, radiusOfCurvature, centreThickness, lensMedium, iface);
	}
		
	public static SimpleDoubleConvexLens fromRadiusOfCurvAndCentreThickness(String name, double centre[], double normal[], double radius, double radiusOfCurvature, double centreThickness, Medium lensMedium, Interface iface) {
		
		return new SimpleDoubleConvexLens(name, centre,  normal, radius, radiusOfCurvature, centreThickness, lensMedium, iface);
	}
	
	private SimpleDoubleConvexLens(String name, double centre[], double normal[], double radius, double radiusOfCurvature, double centreThickness, Medium lensMedium, Interface iface) {
		super(name);
		this.centre = centre; 
		this.lensMedium = lensMedium;
		this.radiusOfCurv = radiusOfCurvature;
		
		double frontCentre[] = new double[3];
		double backCentre[] = new double[3];
		double frontNormal[] = new double[3];
		
		//dishes have their normal facing toward the centre of their radius of curvature 
		for(int i=0;i<3;i++){
			frontCentre[i] = centre[i] + (centreThickness/2) * normal[i];
			backCentre[i] = centre[i] - (centreThickness/2) * normal[i];
			frontNormal[i] = -normal[i];
		}
		
		if(centreThickness < 0){
			backSurface = new Dish(getName() + "-back", frontCentre, frontNormal, radiusOfCurv, radius, null, lensMedium, iface);
			frontSurface = new Dish(getName() + "-front", backCentre, normal, radiusOfCurv, radius, null, lensMedium, iface);
		}else{
			frontSurface = new Dish(getName() + "-front", frontCentre, frontNormal, radiusOfCurv, radius, lensMedium, null, iface);
			backSurface = new Dish(getName() + "-back", backCentre, normal, radiusOfCurv, radius, lensMedium, null, iface);
			
		}

		addElement(frontSurface);
		addElement(backSurface);
				
	}
		
	public static double thickLensFocalLength(double radius1, double radius2, double thickness, double n) {
		return 1.0 / ((n-1.0)*( 1/radius1 - 1/radius2 + (n-1)*thickness/(n*radius1*radius2) ));
	}

	public void setFocalLength(double f, double wavelength) {
		
		double n= lensMedium.getMaterial().getRefractiveIndex(0, wavelength, 300);
		double lensRadius = frontSurface.getDishDiameter()/2;
		
		radiusOfCurv = focalLengthToRadCurv(f, lensRadius, n);
		
		frontSurface.setRadiusOfCurv(radiusOfCurv);
		backSurface.setRadiusOfCurv(radiusOfCurv);
		
	}
	

	public Dish getBackSurface() { return backSurface;	}
	public Dish getFrontSurface() { return frontSurface;	}

	public double getRadius() { return Math.max(frontSurface.getDishDiameter()/2, backSurface.getDishDiameter()/2); }
	
}
