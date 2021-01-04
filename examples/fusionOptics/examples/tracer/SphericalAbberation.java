package fusionOptics.examples.tracer;

import fusionOptics.MinervaOpticsSettings;
import fusionOptics.Util;
import fusionOptics.collection.ImageCollector;
import fusionOptics.drawing.VRMLDrawer;
import fusionOptics.interfaces.Absorber;
import fusionOptics.interfaces.IsoIsoStdFresnel;
import fusionOptics.materials.IsotropicFixedIndexGlass;
import fusionOptics.optics.SimpleDoubleConvexLens;
import fusionOptics.optimisation.OptimiseLensCurvature;
import fusionOptics.pointSpread.DualGaussianPSF;
import fusionOptics.pointSpread.GaussianPSF;
import fusionOptics.pointSpread.MiniImagePSF;
import fusionOptics.pointSpread.PSFGrid;
import fusionOptics.pointSpread.PointSpreadBuilder;
import fusionOptics.pointSpread.PointSpreadFunction;
import fusionOptics.pointSpread.PointsPSF;
import fusionOptics.surfaces.Dish;
import fusionOptics.surfaces.Iris;
import fusionOptics.surfaces.Square;
import fusionOptics.tracer.Tracer;
import fusionOptics.types.Element;
import fusionOptics.types.Material;
import fusionOptics.types.Medium;
import fusionOptics.types.Optic;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;
import otherSupport.ColorMaps;
import binaryMatrixFile.BinaryMatrixFile;
import binaryMatrixFile.BinaryMatrixWriter;


import oneLiners.OneLiners;
import seed.optimization.HookeAndJeeves;


/** Imaging with a single lens
 * Attempt to reduce spherical abberation by adjusting radii of curvature
 * 
 * 
 * From wiki, best if:
 * 	(r1+r2)/(r1-r2) = 2(n^2 - 1)/(n+2) * (i+o)/(i-o) 
 * 
 * 
 * @author oliford
 *
 */
public class SphericalAbberation {
	public final static String outPath = MinervaOpticsSettings.getAppsOutputPath() + "/rayTracing/sphrAbbr";
	
	public static final double fR0 = 0.95, fR1 = 1.1;
	public static final int nR = 50;
	public static final int nRays = 2000;
	public static final double rt2 = Math.sqrt(2);
	
	public final static double wavelen = 654e-9;
	
	/** Distance to object */
	public static final double dObj = 3.0;
	public static final double dImg = 2.0;
	public static final double f = 1.0 / (1.0/dObj + 1.0/dImg);
	
	public static void main(String[] args) {
		
		System.out.println("f = " + f);
		
		//BinaryMatrixFile.mustWrite(outPath + "/object.bin", obj, true);
		
		Material lensMat = new IsotropicFixedIndexGlass(1.5);
		Medium lensMed = new Medium(lensMat);
		
		Square backPlane = new Square("backPlane", new double[]{ 0, 0, 0 }, new double[]{ 1, 0, 0 }, new double[]{ 0, 1, 0 }, 1, 0.5, Absorber.ideal());
		
		Iris lensIris = new Iris("lensIris", new double[]{ dObj+1.0, 0, 0 }, new double[]{ -1, 0, 0 }, 0.3, 0.2, null, null, Absorber.ideal());
		
		SimpleDoubleConvexLens lens = SimpleDoubleConvexLens.fromFocalLengthAndEdgeThickness("lens", new double[]{  dObj+1.0, 0, 0 }, new double[]{ -1, 0, 0 }, 0.2, f, 0.005, lensMed, IsoIsoStdFresnel.ideal(), wavelen);
		
		Square imgPlane = new Square("imgPlane", new double[]{ dObj+1.0+dImg, 0, 0 }, new double[]{ 1, 0, 0 }, new double[]{ 0, 1, 0 }, 1, 0.5, Absorber.ideal());
		
		Optic all = new Optic("all", new Element[]{ backPlane, lens, lensIris, imgPlane });
		
		OptimiseLensCurvature opt = new OptimiseLensCurvature(new SimpleDoubleConvexLens[]{ lens }, imgPlane, all, new double[]{ 1,0,0 }, wavelen, 1000);
		System.out.println("R0 = " + ((Dish)lens.getSurfaces().get(0)).getRadiusOfCurv());
		opt.optimise(new HookeAndJeeves(null), 100);
		System.exit(0);
		
		PointSpreadBuilder psfCollect = new PointSpreadBuilder(imgPlane);
		
		ImageCollector imgCollect = new ImageCollector(
				-0.040, 0.040, 800, 
				-0.030, 0.030, 600);	
		
		BinaryMatrixWriter infoOut = new BinaryMatrixWriter(outPath + "/info.bin", 6);
		
		VRMLDrawer vrmlOut = new VRMLDrawer(outPath + "/sphrAbbr.vrml", 0.01);
		if(vrmlOut != null) {
			vrmlOut.setDrawPolarisationFrames(false);
			vrmlOut.setSkipRays(17);
		}
		
		double col[][] = ColorMaps.alternating2D2x2(nR, nR);
		
		/** shift in image pos for different R's */ 
		double dx = 0.005;
		double R0 = ((Dish)lens.getSurfaces().get(0)).getRadiusOfCurv();
		
		double fwhm[][] = new double[nR][nR];
			
		for(int iR1 = 0; iR1 < nR; iR1++){
			double R1 = R0 * (fR0 + iR1 * (fR1 - fR0)/(nR-1.0));
			for(int iR2 = 0; iR2 < nR; iR2++){
				double R2 = R0 * (fR0 + iR2 * (fR1 - fR0)/(nR-1.0));
				
				((Dish)lens.getSurfaces().get(0)).setRadiusOfCurv(R1);
				((Dish)lens.getSurfaces().get(1)).setRadiusOfCurv(R2);
		
				double startPos[] = new double[]{ 1.0, -0.5*nR*dx+iR1*dx, -0.5*nR*dx+iR2*dx };
				psfCollect.startNewPSF(startPos, new GaussianPSF());
				
				for(int i=0; i < nRays; i++) {
					
					RaySegment ray = new RaySegment();
					ray.startPos = startPos;
					ray.dir = Tracer.generateRandomRayTowardSurface(ray.startPos, lens.getSurfaces().get(0));
					
					ray.length = Double.POSITIVE_INFINITY;
					ray.up = Util.cross(Util.reNorm(Util.cross(ray.dir, new double[]{0,0,1})), ray.dir);
					
					
					ray.E0 = PointSpreadFunction.getInputStatesForMuellerCalc();
					ray.wavelength = wavelen;
					
					Tracer.trace(all, ray, 10, 0.1, false);
					
					if(vrmlOut != null)
						vrmlOut.drawRay(ray, col[iR2*nR+iR1]);
					
					int nHits = ray.processIntersections(imgPlane, imgCollect, psfCollect);
					psfCollect.nextCoherentSet();
					
					Pol.recoverAll();
				}
				
				GaussianPSF psf = (GaussianPSF)psfCollect.psfDone(false);
				
				System.out.println(R0 + "\t" + R1/R0 + "\t" + R2/R0 + "\t" + 
									psf.getMeanX() + "\t" + psf.getMeanY() + "\t" + psf.getFWHM());
						
				infoOut.writeRow(R0, R1/R0, R2/R0, psf.getMeanX(), psf.getMeanY(), psf.getFWHM());
				
				fwhm[iR2][iR1] = psf.getFWHM();
			}

		}
		
		infoOut.close();
		BinaryMatrixFile.mustWrite(outPath + "/fwhm.bin", 
				OneLiners.linSpace(fR0*R0, fR1*R0, nR), 
				OneLiners.linSpace(fR0*R0, fR1*R0, nR), 
				fwhm, false);
		
		if(vrmlOut != null) {
			vrmlOut.drawOptic(all);
			vrmlOut.destroy();
		}
		
		if(imgCollect != null){
			imgCollect.writeImage(outPath + "/image.bin");
		}
	}
	
	private static double[][] objectPointsFromGrid(double minX, double maxX, double minY, double maxY, int nX, int nY) {
		double obj[][] = new double[2][nX*nY];
		double dX = (maxX - minX) / (nX - 1.0);
		double dY = (maxY - minY) / (nY - 1.0);
		int k=0;
		for(int iX = 0; iX < nX; iX++){
			for(int iY = 0; iY < nY; iY++){
				obj[0][k] = minX + iX * dX;
				obj[1][k] = minY + iY * dY;
				k++;
			}
		}
		System.out.println(obj[0].length + ", " + k);
		return obj;
	}
	
	private static double[][] objectPointsFromLine(double lineX[], double lineY[], int nLinePts) {
		double rangeX[] = OneLiners.getRange(lineX);
		double rangeY[] = OneLiners.getRange(lineY);		
		
		double obj[][] = new double[2][lineX.length*nLinePts+lineY.length*nLinePts];
		for(int i=0; i < lineX.length; i++){
			for(int j=0; j < nLinePts; j++) {
				obj[0][i*nLinePts+j] = lineX[i];
				obj[1][i*nLinePts+j] = rangeX[0] + j * (rangeX[1] - rangeX[0]) / ((double)nLinePts - 1.0);
			}
		}
		for(int i=0;i < lineY.length; i++){
			for(int j=0; j < nLinePts; j++) {
				obj[0][(lineX.length+i)*nLinePts+j] = rangeY[0] + j * (rangeY[1] - rangeY[0]) / ((double)nLinePts - 1.0);
				obj[1][(lineX.length+i)*nLinePts+j] = lineY[i];
			}
		}
		
		return obj;
	}
	
}
