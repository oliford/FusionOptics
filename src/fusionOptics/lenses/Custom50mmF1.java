package fusionOptics.lenses;

import java.util.ArrayList;

import fusionOptics.interfaces.Absorber;
import fusionOptics.interfaces.IsoIsoStdFresnel;
import fusionOptics.materials.IsotropicFixedIndexGlass;
import fusionOptics.surfaces.Aspheric;
import fusionOptics.surfaces.Disc;
import fusionOptics.surfaces.Dish;
import fusionOptics.surfaces.Iris;
import fusionOptics.surfaces.Square;
import fusionOptics.types.Element;
import fusionOptics.types.Medium;
import fusionOptics.types.Optic;
import fusionOptics.types.Surface;

/** Custom 50mm 5 element lens */
public class Custom50mmF1 extends Optic {
	
	
	// Create Materials for arbitary fixed-index glasses, and create a medium for each 
	public static Medium SLAM7 = new Medium(new IsotropicFixedIndexGlass(1.74352613)); // OHARA S-LAM7
	public static Medium ZNSE= new Medium(new IsotropicFixedIndexGlass(2.57887416));   // Index 2.57887416 source:  Handbook of Optics Vol. II (need to recheck)
	public static Medium NSF57 = new Medium(new IsotropicFixedIndexGlass(1.83690172)); // SCHOTT N-SF57
	
	
	// Have to reverse all coefficients for negative surfaces?
	//public static double s1Dist = 0.5; //why is the lens 500mm away from 0?
	public static double s1Dist = 0.0;
	
	public static double s1Pos[] =  new double[]{s1Dist, 0.0, 0.0};
	public static double s1Curv = 3748.74e-3;
	public static double s1Norm[] = new double[]{-1.0, 0.0, 0.0};
	public static double s1Thick = 14.52579e-3;
	public static double s1Radius = 27.0e-3;
	public static Medium s1Medium1 = null;
	public static Medium s1Medium2 = NSF57;
	public static double s1ConicConst = 0.0; // As order 06/02/15
	public static double s1PolyCoeffsInMillimeters[] = new double[]{ 0.0, 3.4662607e-006, -9.5847605e-011, -4.6348623e-014 };
	public static double s1PolyCoeffs[] = Aspheric.rescaleCoeffs(s1PolyCoeffsInMillimeters, 1e-3);
	
	public static double s2Dist = s1Dist + s1Thick;
	public static double s2Pos[] =  new double[]{s2Dist, 0.0, 0.0};
	public static double s2Curv = 192.7492e-3;
	public static double s2Norm[] = new double[]{1.0, 0.0, 0.0};
	public static double s2Thick = 26.50819e-3;
	public static double s2Radius = 31.0e-3;
	public static Medium s2Medium1 = null;
	public static Medium s2Medium2 = NSF57;
	public static double s2ConicConst = 0.0; 
	public static double s2PolyCoeffsInMillimeters[] = new double[]{0, -2.8773777e-006, 7.5466219e-010, -5.3131379e-014};
	public static double s2PolyCoeffs[] = Aspheric.rescaleCoeffs(s2PolyCoeffsInMillimeters, 1e-3);

	public static double s3Dist = s2Dist + s2Thick;
	public static double s3Pos[] =  new double[]{s3Dist, 0.0, 0.0};
	public static double s3Curv = 170.8094e-3;
	public static double s3Norm[] = new double[]{-1.0, 0.0, 0.0};
	public static double s3Thick = 0.010;
	public static double s3Radius = 0.040;
	public static Medium s3Medium1 = null;
	public static Medium s3Medium2 = NSF57;
	
	public static double s4Dist = s3Dist + s3Thick;
	public static double s4Pos[] =  new double[]{s4Dist, 0.0, 0.0};
	public static double s4Curv = 91.92364e-3;
	public static double s4Norm[] = new double[]{-1.0, 0.0, 0.0};
	public static double s4Thick = 10.96722e-3;
	public static double s4Radius = 0.041;
	public static Medium s4Medium1 = NSF57;
	public static Medium s4Medium2 = null;
	
	public static double s5Dist = s4Dist + s4Thick;
	public static double s5Pos[] =  new double[]{s5Dist, 0.0, 0.0};
	public static double s5Curv = 138.1855e-3;
	public static double s5Norm[] = new double[]{1.0, 0.0, 0.0};
	public static double s5Thick = 14.92255e-3;
	public static double s5Radius = 46.5e-3;
	public static Medium s5Medium1 = ZNSE;
	public static Medium s5Medium2 = null;
	
	public static double s6Dist = s5Dist + s5Thick;
	public static double s6Pos[] =  new double[]{s6Dist, 0.0, 0.0};
	public static double s6Curv = 1046.628e-3;
	public static double s6Norm[] = new double[]{-1.0, 0.0, 0.0};
	public static double s6Thick = 2.0e-3;
	public static double s6Radius = 46.0e-3;
	public static Medium s6Medium1 = ZNSE;
	public static Medium s6Medium2 = null;
	
	public static double s7Dist = s6Dist + s6Thick;
	public static double s7Pos[] =  new double[]{s7Dist, 0.0, 0.0};
	public static double s7Curv = 41.44526e-3;
	public static double s7Norm[] = new double[]{1.0, 0.0, 0.0};
	public static double s7Thick = 14.9948e-3;
	public static double s7Radius = 35.0e-3;
	public static Medium s7Medium1 = ZNSE;
	public static Medium s7Medium2 = null;
	public static double s7ConicConst = 0.0; // As order 06/02/15
	public static double s7PolyCoeffsInMillimeters[] = new double[]{ 0, -1.4887119e-008, -1.5803008e-011, -1.4029144e-015 };
	public static double s7PolyCoeffs[] = Aspheric.rescaleCoeffs(s7PolyCoeffsInMillimeters, 1e-3);

	public static double s8Dist = s7Dist + s7Thick;
	public static double s8Pos[] =  new double[]{s8Dist, 0.0, 0.0};
	public static double s8Curv = 31.02156e-3;
	public static double s8Norm[] = new double[]{1.0, 0.0, 0.0};
	public static double s8Thick =  14.10915e-3;
	public static double s8Radius = 26e-3;
	public static Medium s8Medium1 = null;
	public static Medium s8Medium2 = ZNSE;
	
	public static double s9Dist = s8Dist + s8Thick;
	public static double s9Pos[] =  new double[]{s9Dist, 0.0, 0.0};
	public static double s9Curv = 146.0394e-3;
	public static double s9Norm[] = new double[]{1.0, 0.0, 0.0};
	public static double s9Thick = 9.782205e-3;
	public static double s9Radius = 24.0e-3;
	public static Medium s9Medium1 = SLAM7;
	public static Medium s9Medium2 = null;
	
	public static double s10Dist = s9Dist + s9Thick;
	public static double s10Pos[] =  new double[]{s10Dist, 0.0, 0.0};
	public static double s10Curv = 242.5449e-3;
	public static double s10Norm[] = new double[]{1.0, 0.0, 0.0};
	//static double s10Thick = 0.0;
	public static double s10Radius = 20.0e-3;
	public static Medium s10Medium1 = null;
	public static Medium s10Medium2 = SLAM7;
	public static double s10ConicConst = 0.0; 
	public static double s10PolyCoeffsInMillimeters[] = new double[]{ 0, 2.1353608e-006, 2.2377491e-010, 2.0783215e-012 };
	public static double s10PolyCoeffs[] = Aspheric.rescaleCoeffs(s10PolyCoeffsInMillimeters, 1e-3);

	public static double backFocalDistance = s10Dist + 17.5304e-3;
	
	public static double irisAperture = 50.0e-3;
	public static double irisDiameter = 100.0e-3;
	
	
	public Iris iris1 = new Iris("iris1", new double[] {s1Dist - 10.078e-3, 0.0, 0.0}, new double[]{ -1.0, 0.0, 0.0},  irisDiameter/2, irisAperture/2, Absorber.ideal());
	//public Iris iris2 = new Iris("iris2", s1Pos, new double[]{ -1.0, 0.0, 0.0}, s1Radius - 0.005, 0.001, Absorber.ideal());
	//public Iris iris3 = new Iris("iris3", new double[] {s8Dist + 0.012, 0.0, 0.0}, new double[]{ -1.0, 0.0, 0.0}, s8Radius + 0.05, 50.30929e-3/2.0, Absorber.ideal());
	
	//new Aspheric(name + "-curved", curvedSurfaceCentre, reverseAxis, -radiusOfCurvatureFront, 1.0001 * clearAperture / 2, conicConst, polyCoeffs, null, medium, IsoIsoInterface.ideal());
	public Aspheric lens1Front = new Aspheric("L1-front", s1Pos, s1Norm, s1Curv, s1Radius, s1ConicConst, s1PolyCoeffs, s1Medium1, s1Medium2, IsoIsoStdFresnel.ideal());
	//Dish lens1Front = new Dish("L1-front", s1Pos, s1Norm, s1Curv, s1Radius, s1Medium1, s1Medium2, IsoIsoStdFresnel.ideal());
	//Dish lens1Back = new Dish("L1-back", s2Pos, s2Norm, s2Curv, s2Radius, s2Medium1, s2Medium2, IsoIsoStdFresnel.ideal());
	public Aspheric lens1Back = new Aspheric("L1-back", s2Pos, s2Norm, s2Curv, s2Radius, s2ConicConst, s2PolyCoeffs, s2Medium1, s2Medium2, IsoIsoStdFresnel.ideal());
	public Dish lens2Front = new Dish("L2-front", s3Pos, s3Norm, s3Curv, s3Radius, s3Medium1, s3Medium2, IsoIsoStdFresnel.ideal());
	public Dish lens2Back = new Dish("L2-back", s4Pos, s4Norm, s4Curv, s4Radius, s4Medium1, s4Medium2, IsoIsoStdFresnel.ideal());
	public Dish lens3Front = new Dish("L3-front", s5Pos, s5Norm, s5Curv, s5Radius, s5Medium1, s5Medium2, IsoIsoStdFresnel.ideal());


	public Dish lens3Back = new Dish("L3-back", s6Pos, s6Norm, s6Curv, s6Radius, s6Medium1, s6Medium2, IsoIsoStdFresnel.ideal());
	//public Dish lens4Front = new Dish("L3-middle", s7Pos, s7Norm, s7Curv, s7Radius, s7Medium1, s7Medium2, IsoIsoStdFresnel.ideal());
	public Aspheric lens4Front = new Aspheric("L4-front", s7Pos, s7Norm, s7Curv, s7Radius, s7ConicConst, s7PolyCoeffs, s7Medium1, s7Medium2, IsoIsoStdFresnel.ideal());
	public Dish lens4Back = new Dish("L4-back", s8Pos, s8Norm, s8Curv, s8Radius, s8Medium1, s8Medium2, IsoIsoStdFresnel.ideal());
	public Dish lens5Front = new Dish("L5-front", s9Pos, s9Norm, s9Curv, s9Radius, s9Medium1, s9Medium2, IsoIsoStdFresnel.ideal());
	//public Dish lens5Back = new Dish("L4-back", s10Pos, s10Norm, s10Curv, s10Radius, s10Medium1, s10Medium2, IsoIsoStdFresnel.ideal());
	public Aspheric lens5Back = new Aspheric("L5-back", s10Pos, s10Norm, s10Curv, s10Radius, s10ConicConst, s10PolyCoeffs, s10Medium1, s10Medium2, IsoIsoStdFresnel.ideal());
	
	
	
	//public Square imagePlane = new Square("imagePlane", new double[]{ imageDistance, 0.0, 0.0}, new double[]{ 1.0, 0.0, 0.0}, new double[]{ 0.0, 0.0, 1.0 }, 0.50, 0.50, Absorber.ideal());
	
	public Custom50mmF1() {
		this(false);
	}
		// the 'all' Optic just contains all the other elements (optics and surfaces) 
	public Custom50mmF1(boolean addIrises) {
		super("Custom50mmF11");

		addElement(iris1);
		addElement(lens1Front);
		addElement(lens1Back);
		addElement(lens2Front);
		addElement(lens2Back);
		addElement(lens3Front); 
		addElement(lens3Back);
		addElement(lens4Front); 
		addElement(lens4Back);
		addElement(lens5Front); 
		addElement(lens5Back);
		
		double bodyDiamter = 0.150;

		//we require much higher than usual tolerance for this lens to work properly
		for(Surface s : new ArrayList<Surface>(getSurfaces())){
			if(s instanceof Aspheric){
				((Aspheric)s).setSurfaceFindingTolerance(1e-10, 2000);
			}
			
			if(addIrises){
				double c[] = null, r=Double.NaN;
				if(s instanceof Disc){
					c = ((Disc)s).getCentre();
					r = ((Disc)s).getRadius();
				}else if(s instanceof Dish){
					c = ((Dish)s).getRimCentre();
					r = ((Dish)s).getDishDiameter()/2;
				}
				
				if(c != null){
					addElement(new Iris("iris_" + s.getName(), 
							c, 
							new double[]{1,0,0}, 
							bodyDiamter/2, 
							r * 0.98,
							Absorber.ideal()));
				}
				
			}
		}

		//addElement(imagePlane);
	}

}
