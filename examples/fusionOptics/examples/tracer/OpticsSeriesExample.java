package fusionOptics.tracer;


import fusionOptics.MinervaOpticsSettings;
import fusionOptics.Util;
import fusionOptics.drawing.AsciiOutForWendel;
import fusionOptics.drawing.SVGRayDrawing;
import fusionOptics.drawing.VRMLDrawer;
import fusionOptics.interfaces.Absorber;
import fusionOptics.interfaces.IsoIsoInterface;
import fusionOptics.interfaces.IsoIsoStdFresnel;
import fusionOptics.interfaces.Reflector;
import fusionOptics.materials.IsotropicFixedIndexGlass;
import fusionOptics.optics.Box;
import fusionOptics.optics.SimpleDoubleConvexLens;
import fusionOptics.surfaces.Cylinder;
import fusionOptics.surfaces.Dish;
import fusionOptics.surfaces.Iris;
import fusionOptics.surfaces.Square;
import fusionOptics.surfaces.Triangle;
import fusionOptics.tracer.Tracer;
import fusionOptics.types.Element;
import fusionOptics.types.Material;
import fusionOptics.types.Medium;
import fusionOptics.types.Optic;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;
import fusionOptics.types.Surface;
import otherSupport.ColorMaps;

/** Test rays through a simple series of optics:
 * 1) An objective lens gets the rays roughly parallel.
 * 2) An iris cuts out some of the rays.
 * 3) Another lens and a spherical mirror act together as the imaging lens and turn the image through 90'
 * 4) A screen above is roughly at the imaging plane.
 * A back screen catches light that misses the mirror.
 *  
 *  */
public class OpticsSeriesExample {
	
	/** Output path for drawings and data */
	final static String outPath = MinervaOpticsSettings.getAppsOutputPath() + "/rayTracing";
	
	public static void main(String[] args) {
		/** Number of start positions (object plane) */
		int nStarts = 4;
		/** Number of rays from each start position */
		int nRays = 100;
		
		// Wavelength / m 
		double wavelen = 500e-9;
		
		// Create Materials for arbitary fixed-index glasses, and create a medium for each 
		Medium glass120 = new Medium(new IsotropicFixedIndexGlass(1.2));
		Medium glass130 = new Medium(new IsotropicFixedIndexGlass(1.3));
		Medium glass105 = new Medium(new IsotropicFixedIndexGlass(1.05));
				
		// Create SVG projections and a VRML output  
		SVGRayDrawing svg = new SVGRayDrawing(outPath + "/simpleOpticsSeries",  new double[]{ 0, -1, -1, 6, 1, 1 }, false);
		VRMLDrawer vrmlOut = new VRMLDrawer(outPath + "/simpleOpticSeries.vrml", 0.0005);
		vrmlOut.setDrawPolarisationFrames(true);
		
		//Box subject = new Box("subject", new double[]{ -0.2, -0.2, -0.2}, new double[]{ 0.2, 0.2, 0.2}, glass120, IsoIsoStdFresnel.ideal());
		
		// Create the optical elements 
		SimpleDoubleConvexLens lens1 = SimpleDoubleConvexLens.fromFocalLengthAndEdgeThickness("lens1", new double[]{ 2.0, 0.0, 0.0 }, new double[]{ -1.0, 0.0, 0.0}, 0.4, 2.0, 0.005, glass130, IsoIsoStdFresnel.ideal(), wavelen);
		Iris iris = new Iris("iris", new double[]{ 2.5, 0.0, 0.0 }, new double[]{ -1.0, 0.0, 0.0}, 0.3, 0.2, Absorber.ideal());
		SimpleDoubleConvexLens lens2 = SimpleDoubleConvexLens.fromFocalLengthAndEdgeThickness("lens2", new double[]{ 3.0, 0.0, 0.0 }, new double[]{ -1.0, 0.0, 0.0}, 0.5, 5.0, 0.005, glass120, IsoIsoStdFresnel.ideal(), wavelen);
		Dish mirror = new Dish("mirror", new double[]{ 4.0, 0.0, 0.0}, new double[]{ -0.707, 0.0, 0.707}, 4.0, 0.5, Reflector.ideal());
		Square screen = new Square("screen", new double[]{ 4.0, 0.0, 1.0}, new double[]{ 0.0, 0.0, -1.0}, new double[]{ 1.0, 0.0, 0.0}, 1, 1, Absorber.ideal());
		Square backstop = new Square("backstop", new double[]{ 6.0, 0.0, 0.0}, new double[]{ -1.0, 0.0, 0.0}, new double[]{ 0.0, 1.0, 0.0}, 2, 2, Absorber.ideal());
		
		Cylinder other = new Cylinder("other", new double[]{ 1.0, 0, 0}, new double[]{ 0, 1, 0}, 0.1, 1, glass105, null, IsoIsoStdFresnel.ideal());
		
		//Triangle other = new Triangle("other", new double[]{0.5,0,0}, new double[]{0.7,0.1,0}, new double[]{0.7,0,0.1}, Reflector.ideal());
		
		// the 'all' Optic just contains all the other elements (optics and surfaces) 
		Optic all = new Optic("all", new Element[]{ other, lens1, iris, lens2, mirror, screen, backstop });
		

		//Create colour table for lines in SVG output
		double colTab[][] = ColorMaps.jet(nStarts);
		svg.generateLineStyles(colTab, 0.001);

		for(int i=0; i < nStarts; i++){
			double startZ = (nStarts <= 1) ? 0 : ((double)i * 0.3)/(nStarts-1);
			
			svg.getSVG3D().startGroup("start"+i); //put all rays from each start-point in an SVG group
			vrmlOut.startGroup("start"+i);
			for(int j=0;j < nRays;j++){
				
				RaySegment ray = new RaySegment(); 
					
				//start the ray off heading randomly toward the first lens primary surface
				ray.startPos = new double[]{ 0, 0, startZ };				
				ray.dir = Tracer.generateRandomRayTowardSurface(ray.startPos, lens1.getSurfaces().get(1));
				
				//ray starts with infinte elength and no end-hit. These will get filled by the tracer.
				ray.length = Double.POSITIVE_INFINITY;
				ray.startHit = null;
				ray.E0 = new double[][]{ {1,0,0,1} }; // linear polarisation at 45deg - not important for this example.
				ray.wavelength = wavelen;
				ray.up = Util.createPerp(ray.dir); //sense of theta=0 for polarisation. Not important here.
				
				// Trace it through everything until it hits an absorber, or after 30 hits, or it gets below 1e-4 intensity
				Tracer.trace(all, ray, 30, 1e-4, true);
				
				//draw the line
				svg.drawRay(ray, i);
				vrmlOut.drawRay(ray, colTab[i]);
				
				// This must be called once the ray tree is finished with.
				Pol.recoverAll();
			}
			svg.getSVG3D().endGroup();
			vrmlOut.endGroup();
		}
		
		svg.drawElement(all);
		svg.destroy();

		vrmlOut.drawOptic(all);
		vrmlOut.destroy();
	}
}
