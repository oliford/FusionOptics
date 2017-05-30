package fusionOptics.lenses;

import fusionOptics.interfaces.IsoIsoInterface;
import fusionOptics.materials.IsotropicFixedIndexGlass;
import fusionOptics.optics.SimplePlanarAsphericLens;
import fusionOptics.surfaces.Aspheric;
import fusionOptics.types.Interface;
import fusionOptics.types.Material;
import fusionOptics.types.Medium;

public class EdmundOptics50mmAspheric extends SimplePlanarAsphericLens {

	/* Edmund Optics 50mm Dia., 0.50 NA, Uncoated, Calibration Grade Aspheric Lens
	* [ http://www.edmundoptics.com/optics/optical-lenses/aspheric-lenses/calibration-grade-precision-aspheric-lenses/69144 
	*   http://www.edmundoptics.com/techsupport/resource_center/product_docs/prnt_69144.pdf ]
	*/ 
	
	public static double designWavelen = 587.6e-9;
	public final static Material mat = new IsotropicFixedIndexGlass(1.58913); //Ohara L-BAL35 at 587.6
	final static double rescale = 1e-3; //everything here is in mm, we want m
	//final static double focalLength = 50.00 / scale;
	public static double diameter = 46.00 * rescale;
	public static double thickness = 19.40 * rescale;
	public final static double backFocalDistance = 37.79 * rescale;	
	public static double curvatureRadius = 29.456 * rescale; // --> R=29.45mm
	public static double conicConstant = -1.431442E+00;
	public static double polyCoeffs_mm[] = new double[]{
			0,					// R^2
			4.289678E-06,		// R^4
			-1.553732E-10,		// R^6
			1.963218E-13,		// R^8
			-7.350410E-17,		// R^10
	};
		
	public EdmundOptics50mmAspheric(double[] pos, double normal[]) {
		super("edmundOptics50mmAspheric",
				pos, 
				normal, 
				diameter / 2, 
				curvatureRadius, 
				thickness, 
				conicConstant, 
				Aspheric.rescaleCoeffs(polyCoeffs_mm, rescale), 
				new Medium(mat), 
				IsoIsoInterface.ideal()); 
	}

	
}
