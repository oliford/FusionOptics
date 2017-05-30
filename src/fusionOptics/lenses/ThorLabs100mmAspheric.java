package fusionOptics.lenses;

import fusionOptics.interfaces.IsoIsoInterface;
import fusionOptics.materials.BK7;
import fusionOptics.materials.IsotropicFixedIndexGlass;
import fusionOptics.optics.SimplePlanarAsphericLens;
import fusionOptics.surfaces.Aspheric;
import fusionOptics.types.Interface;
import fusionOptics.types.Material;
import fusionOptics.types.Medium;

public class ThorLabs100mmAspheric extends SimplePlanarAsphericLens {

	/* Thor Labs AL50100-A N-BK7  350 - 700 nm, EFL=100.0mm, NA=0.23 Ø=50.0mm
	 *  
	 * https://www.thorlabs.de/newgrouppage9.cfm?objectgroup_id=7176#
	*/
	
	public static double designWavelen = 780.0e-9;
	public final static Material mat = new BK7(); //It's really BK7
	final static double rescale = 1e-3; //everything here is in mm, we want m
	//final static double focalLength = 50.00 / scale;
	public static double diameter = 50.00 * rescale;
	public static double thickness = 10.00 * rescale;
	public final static double backFocalDistance = 93.4 * rescale;	
	public static double curvatureRadius = 51.12 * rescale; // --> R=29.45mm
	public static double conicConstant = -0.575;
	public static double polyCoeffs_mm[] = new double[]{
			0,					// R^2
			-4.8366264e-11,		// R^4
			-8.5756915e-12,		// R^6
			-2.0138223e-15,		// R^8
			-4.5977971e-19,		// R^10
	};
		
	public ThorLabs100mmAspheric(double[] pos, double normal[]) {
		super("thorLabs100mmAspheric",
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
