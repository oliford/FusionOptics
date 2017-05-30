package fusionOptics.load;

import fusionOptics.MinervaOpticsSettings;
import fusionOptics.Util;
import fusionOptics.drawing.SVGRayDrawing;
import fusionOptics.drawing.VRMLDrawer;
import fusionOptics.interfaces.Absorber;
import fusionOptics.interfaces.NullInterface;
import fusionOptics.lenses.Nikon135mmF28;
import fusionOptics.lenses.SchneiderXenon25mmF095;
import fusionOptics.load.ZemaxXMZFile;
import fusionOptics.optics.SequentialLensSeries;
import fusionOptics.optimisation.OptimiseIndices;
import fusionOptics.surfaces.Square;
import fusionOptics.tracer.Tracer;
import fusionOptics.types.Material;
import fusionOptics.types.Optic;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;
import otherSupport.ColorMaps;

/** Example for loading Zeemax lens descriptions into minerva-optics 
 * and using them to create an image. */
public class LoadZemaxFile {
	final static String outPath = MinervaOpticsSettings.getAppsOutputPath() + "/rayTracing/loadZMX";
	
	/** Number of start positions for rays */
	public final static int nStarts = 15;	
	/** Number of rays for each start position */
	public final static int nRays = 100;		
	/** Lens name */
	
	
	//public final static String lens = "Schneider75mm-f0.95";	
	/** Distance between object and lens */
	//public final static double objectDist = 3.000;	
	/** Distance between image and lens */
	//public final static double imageDist = 0.025;	
	/** Size of object to image */
	//public final static double objectSize = 1.3;
	//*/
	
	/* Alternate lens definitions: */
	
	// SequentialLensSeries 1 
	/*public final static String lens = "Nikon135mm-f2.8";
	public final static double objectDist = 3.000;
	public final static double imageDist = 0.135;
	public final static double objectSize = 1.0;
	
	// SequentialLensSeries 2 (50mm f/1.1)
	/*public final static String lens = "Nikon135mm-f2.8";
	public final static double objectDist = 3.000;
	public final static double imageDist = 0.050;
	public final static double objectSize = 1.5;
	//*/
	
	/*public final static String lens = "Nikon28-200mm";
	public final static double objectDist = 3.000;
	public final static double imageDist = 0.180;
	public final static double objectSize = 0.5;
	//*/
	
	public final static String lens = "Nikonf105mmN2";
	public final static double objectDist = 3.000;
	public final static double imageDist = 0.11464;
	public final static double objectSize =  2;
	//*/
	
	//String lens = "Canon200mmf2.8";
	//*/
	
	public final static double wavelength = 660e-9;
	
	public static void main(String[] args) {
		
		// with thanks to Nick Konidaris for this: http://sites.google.com/site/nickkonidaris/prescriptions 
		ZemaxXMZFile zmxFile = new ZemaxXMZFile("/work/ipp/scrap/zemax-lenses/"+lens+".ZMX");
		zmxFile.setIrisOuterRadius(0.03);
		zmxFile.load();

		// The optical system we'll be tracing
		Optic sys = zmxFile.getMainOptic();
		
		// Alternatively, use a class defined lens:  
		// SchneiderXenon25mmF095 sys = new SchneiderXenon25mmF095();
		
		// SVG and VRML drawings
		VRMLDrawer vrmlOut = new VRMLDrawer(outPath + "/"+lens+".vrml", 0.001);
		SVGRayDrawing svgOut = new SVGRayDrawing(outPath + "/"+lens, new double[]{ -3,-0.2,-0.2, 1,0.2,0.2}, true);

		Square imagePlane = new Square("imagePlane", new double[]{ imageDist, 0, 0}, new double[]{ 1,0,0 }, new double[]{ 0,0,1 }, 0.1, 0.1, null, null, Absorber.ideal());
		
		sys.addElement(imagePlane);
		
		for(int iS=0; iS < nStarts; iS++){
			double ang = -(iS * objectSize / (nStarts - 1)) / objectDist;
			sys.addElement(new Square("pos"+iS, 
					new double[]{ imageDist + 0.005, 0, ang * imageDist },
					new double[]{ 0, 0, 1 },
					new double[]{ 1, 0, 0 },
					0.010, 0.010, null, null, NullInterface.ideal()));
		}
		
		svgOut.setPrecision(8);
		double cols[][] = ColorMaps.jet(nStarts);
		svgOut.generateLineStyles(cols, 0.0001);
		
		//  Code for autofocusing the image plane 
		/*AutoFocus af = new AutoFocus();
		OptimiseIndices af = new OptimiseIndices(0.025);
		af.setTracingElements(sys, imagePlane);
		//af.initRaysImaging(sys.getSurfaces().get(0), new double[]{ -objectDist,0,0}, new double[]{ 0,0,0.1 }, 10, wavelength, 500);
		af.initRaysParallel(sys.getSurfaces().get(0), new double[]{1,0,0}, 0.7*(objectSize/objectDist), 0.2, wavelength, 10, 200);
		//af.setMovement(imagePlane, new double[]{1,0,0}, -0.1, 0.1);
		af.setModifications(new Material[]{
				sys.getMedia().get(0).getMaterial(),
				sys.getMedia().get(1).getMaterial(),
				sys.getMedia().get(2).getMaterial(),
				sys.getMedia().get(3).getMaterial(),
				sys.getMedia().get(4).getMaterial(),
				sys.getMedia().get(5).getMaterial(),
				sys.getMedia().get(6).getMaterial(),
				sys.getMedia().get(7).getMaterial(),
				}, 
				new double[]{ 1.20, 1.20, 1.20, 1.20, 1.20, 1.20, 1.20, 1.20 }, 
				new double[]{ 1.70, 1.70, 1.70, 1.70, 1.70, 1.70, 1.70, 1.70 });
		
		//af.setSVGOut(outPath + "/autoFocus", 77);
		af.optimise(10000);//*/
	
		
		svgOut.startGroup("rays");
		vrmlOut.startGroup("rays");
		for(int iS=0; iS < nStarts; iS++){
			double imagePos[] = new double[]{ -objectDist, 0, iS * objectSize / (nStarts - 1) };
			vrmlOut.startGroup("start" + iS);
			svgOut.startGroup("start" + iS);
			
			double infinityDir[] = Util.minus(sys.getSurfaces().get(0).getBoundarySphereCentre(), imagePos);
			
			for(int iR=0; iR < nRays; iR++){
				
				//double dir[] = Tracer.generateRandomRayTowardSurface(startPos, sys.getSurfaces().get(0));
				//double rayStart = imagePos;
				
				double dir[] = infinityDir;
				double rayStart[] = Tracer.generateRandomInfinityRayStartPos(dir, sys.getSurfaces().get(0), 1);
				
				RaySegment ray = new RaySegment(rayStart, dir, wavelength);
												
				Tracer.trace(sys, ray, 1000, 0.001, false);
				
				svgOut.drawRay(ray, iS);
				vrmlOut.drawRay(ray, cols[iS]);
				
				Pol.recoverAll();
			}
			
			svgOut.endGroup();
			vrmlOut.endGroup();
		}
		svgOut.endGroup();
		vrmlOut.endGroup();
		
		vrmlOut.drawOptic(sys);
		svgOut.drawElement(sys);
		
		vrmlOut.destroy();
		svgOut.destroy();
		
	}
}
