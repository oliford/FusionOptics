package fusionOptics.examples.tracer;

import net.jafama.FastMath;

import java.util.List;

import fusionOptics.MinervaOpticsSettings;
import fusionOptics.Util;
import fusionOptics.drawing.VRMLDrawer;
import fusionOptics.interfaces.Absorber;
import fusionOptics.interfaces.IsoIsoStdFresnel;
import fusionOptics.interfaces.IsoUniaxialInterface;
import fusionOptics.interfaces.SimplePolariser;
import fusionOptics.lenses.Nikon135mmF28;
import fusionOptics.lenses.Nikon50mmF11;
import fusionOptics.lenses.SchneiderXenon25mmF095;
import fusionOptics.materials.Calcite;
import fusionOptics.materials.CrystalQuartz;
import fusionOptics.materials.IsotropicFixedIndexGlass;
import fusionOptics.materials.LithiumNiobate;
import fusionOptics.materials.UniaxialFixedIndexGlass;
import fusionOptics.optics.Box;
import fusionOptics.optics.SimpleDoubleConvexLens;
import fusionOptics.pointSpread.DualGaussianPSF;
import fusionOptics.pointSpread.GaussianPSF;
import fusionOptics.pointSpread.MiniImagePSF;
import fusionOptics.pointSpread.PSFGrid;
import fusionOptics.pointSpread.PointSpreadBuilder;
import fusionOptics.pointSpread.PointSpreadFunction;
import fusionOptics.pointSpread.PointsPSF;
import fusionOptics.surfaces.Iris;
import fusionOptics.surfaces.Plane;
import fusionOptics.surfaces.Square;
import fusionOptics.tracer.Tracer;
import fusionOptics.types.Element;
import fusionOptics.types.Intersection;
import fusionOptics.types.Material;
import fusionOptics.types.Medium;
import fusionOptics.types.Optic;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;

import otherSupport.ColorMaps;

import binaryMatrixFile.BinaryMatrixFile;
import binaryMatrixFile.BinaryMatrixWriter;



import oneLiners.OneLiners;


/** Savart test with parallel rays.
 * Also will later test the PSF collection coherent interference.
 * 
 * This doesn't currently work, due to phase problems in IsoUniaxialInterface.
 * 
 * @author oliford
 */
public class SavartTest {
	final static String outPath = MinervaOpticsSettings.getAppsOutputPath() + "/rayTracing/savartTest";
	
	final static double maxTheta = 30 * Math.PI / 180;
	final static int nRaysPerSource = 1000;
	final static double rt2 = Math.sqrt(2);
	
	
	final static double gridDef[][] = {
		//{ 0.5, 1.5, 5 },
			{ 0.98, 1.02, 3 }, 
			{ -0.05, 0.05, 20 }, 
			{ -0.05, 0.05, 20 }
		};
	
	final static double imageX[] = OneLiners.linSpace(-0.3, 0.3, 500);
	final static double imageY[] = OneLiners.linSpace(-0.3, 0.3, 500);
	
	final static double wavelen = 593e-9;
	
	public Material lensMat = new IsotropicFixedIndexGlass(1.3);
	public Medium lensMed = new Medium(lensMat);
	
	public double lensRadA = 0.2, lensRadB = 0.3;
	public Square backPlane = new Square("backPlane", new double[]{ -0.010, 0, 0 }, new double[]{ 1, 0, 0 }, new double[]{ 0, 1, 0 }, 0.100, 0.050, Absorber.ideal());
	
	//public SimpleDoubleConvexLens objLens = new SimpleDoubleConvexLens("objLens", new double[]{ 4, 0, 0 }, new double[]{ -1, 0, 0 }, lensRadA, lensMed, IsoIsoStdFresnel.ideal(), 3, wavelen);
	double objLensPos = 0.000 + 0.135; 
	public Nikon135mmF28 objLens = new Nikon135mmF28(new double[]{ objLensPos, 0, 0 });
	public Iris objLensIris = new Iris("objLensIris", new double[]{ objLensPos, 0, 0 }, new double[]{ -1, 0, 0 }, 0.050, objLens.getCaseRadius(), null, null, Absorber.ideal());
	
	public Square pol1 = new Square("pol1", new double[]{ 0.250, 0, 0}, new double[]{ -1, 0, 0}, new double[]{ 0, rt2/2, rt2/2}, 0.080, 0.080, 
								new SimplePolariser(new double[]{ 0, rt2/2, rt2/2 }, 0));
	
	public Box savart = new Box("savart", new double[]{ 0.300, 0, 0 }, 0.002, 0.080, 0.080, 
			new Medium(
					//new UniaxialFixedIndexGlass(1.5, 1.6, 1.0),
					new LithiumNiobate(),
					//new double[][]{{ rt2/2, 0, rt2/2 }},
					new double[][]{{ 0,0,1 }},
					300 ),
			IsoUniaxialInterface.ideal());
	

	public Square pol2 = new Square("pol1", new double[]{ 0.350, 0, 0}, new double[]{ -1, 0, 0}, new double[]{ 0, rt2/2, rt2/2}, 0.080, 0.080, 
								new SimplePolariser(new double[]{ 0, rt2/2, rt2/2 }, 0));
	
	//public Iris imgLensIris = new Iris("imgLensIris", new double[]{ 6, 0, 0 }, new double[]{ -1, 0, 0 }, 1.5*lensRadB, lensRadB, null, null, Absorber.ideal());
	//public SimpleDoubleConvexLens imgLens = new SimpleDoubleConvexLens("imgLens", new double[]{ 6, 0, 0 }, new double[]{ -1, 0, 0 }, lensRadB, lensMed, IsoIsoStdFresnel.ideal(), 2, wavelen);
	double imgPlanePos = 0.500; 
	double imgLensPos = imgPlanePos - 0.050;
	public Nikon50mmF11 imgLens = new Nikon50mmF11(new double[]{ imgLensPos, 0, 0 });
	public Iris imgLensIris = new Iris("imgLensIris", new double[]{ imgLensPos, 0, 0 }, new double[]{ -1, 0, 0 }, 0.1, imgLens.getCaseRadius(), null, null, Absorber.ideal());
	
	public Square imgPlane = new Square("imgPlane", new double[]{ imgPlanePos, 0, 0 }, new double[]{ -1, 0, 0 }, new double[]{ 0, 1, 0 }, 0.100, 0.050, Absorber.ideal());
	
	public Optic all = new Optic("all", new Element[]{ backPlane, objLens, objLensIris, pol1, savart, pol2, imgLens, imgLensIris, imgPlane });
	
	public SavartTest() {
		Util.rotateOnZ(objLens, new double[]{ objLensPos, 0, 0 }, Math.PI);
	}
	
	public static void main(String[] args) {
		//(new SavartTest()).buildPSFGrid();
		(new SavartTest()).fringeScan();
	}
	
	public void fringeScan() {		
		int nRaysPerSource = 100;
		int startsDrawSkip = 200;
		
		//Util.rotateOnZ(savart, savart.getBoundarySphereCentre(), 1 * Math.PI / 180);
		Util.rotateOnZ(savart.getSurfaces().get(0), 
				savart.getSurfaces().get(0).getBoundarySphereCentre(), 0.7 * Math.PI / 180);
		Util.rotateOnY(savart.getSurfaces().get(0), 
				savart.getSurfaces().get(0).getBoundarySphereCentre(), 0.4 * Math.PI / 180);
		
		double scanZ[] = OneLiners.linSpace(-0.010, 0.010, 1000);
		double ang[] = new double[scanZ.length];
		double img[] = new double[scanZ.length];
		double hitC[] = new double[scanZ.length];
		
		double a[][] = new double[][]{ 
				OneLiners.linSpace(-5*Math.PI/180, 5*Math.PI/180, 2000),
				new double[2000],
				new double[2000],
				new double[2000],
				new double[2000],
		};
		
		VRMLDrawer vrmlOut = new VRMLDrawer(outPath + "/fringeScan.vrml", 0.001);
		if(vrmlOut != null) {
			vrmlOut.setDrawPolarisationFrames(true);
			vrmlOut.setSkipRays(59);
		}
		double col[][] = ColorMaps.jet(20);
		
		
		//BinaryMatrixWriter dataOut = new BinaryMatrixWriter(outPath + "/scan.bin", 
		//		PointSpreadFunction.inputStatesForMuellerCalc.length+1);
		
		for(int i=0; i < scanZ.length; i++) {
			
			double startPos[] = new double[]{ 0, 0, scanZ[i] };
			for(int j=0; j < nRaysPerSource; j++) {
				
				RaySegment ray = new RaySegment();
				ray.startPos = startPos.clone();
				ray.dir = Util.reNorm(Util.minus(objLens.getBoundarySphereCentre(), ray.startPos));
				ang[i] = Math.asin(ray.dir[2]);
				if(nRaysPerSource > 1)
					ray.dir = Tracer.generateRandomRayTowardSurface(ray.startPos, objLens.getSurfaces().get(0));
				
				ray.length = Double.POSITIVE_INFINITY;
				ray.up = Util.cross(Util.reNorm(Util.cross(ray.dir, new double[]{0,0,1})), ray.dir);
				
				ray.E0 = PointSpreadFunction.getInputStatesForMuellerCalc();
				ray.wavelength = wavelen;
				
				Tracer.trace(all, ray, 100, 0.3, true);
				
				//if(vrmlOut != null)
				//	vrmlOut.drawRay(ray, col[i % col.length]);
				
				//ray.dumpPath();
				//ray.processIntersections(imgPlane, psfBuild);
				List<Intersection> hits = ray.getIntersections(imgPlane);
				double E[][] = new double[ray.E0.length][4];
				for(Intersection hit : hits) {
					hit.incidentRay.rotatePolRefFrame(imgPlane.getUp());
					for(int k=0; k < ray.E0.length; k++){
						E[k][0] += hit.incidentRay.E1[k][0];
						E[k][1] += hit.incidentRay.E1[k][1];
						E[k][2] += hit.incidentRay.E1[k][2];
						E[k][3] += hit.incidentRay.E1[k][3];
					}
					hitC[i]++;
				}				
				img[i] += Pol.intensity(E[0]);
				if(hits.size() > 0 && (i % startsDrawSkip) == 0 && vrmlOut != null){
					vrmlOut.drawRay(ray, col[(i / startsDrawSkip) % col.length]);	
				}
				
				hits = ray.getIntersections(savart.getSurfaces().get(1));
				if(hits.size() > 0){
					Intersection hit = hits.get(0);
					double plateAng = FastMath.atan2(hit.incidentRay.dir[2], hit.incidentRay.dir[0]);
					int idx = OneLiners.getNearestLowerIndexProbablyRegular(a[0], plateAng);
					double Ec[] = new double[4];
					hit.transmittedOrdinary.rotatePolRefFrame(new double[]{0, 0, 1 });
					hit.transmittedExtraordinary.rotatePolRefFrame(new double[]{0, 0, 1 });
					Ec[0] = hit.transmittedOrdinary.E0[0][0] + hit.transmittedExtraordinary.E0[0][0];
					Ec[1] = hit.transmittedOrdinary.E0[0][1] + hit.transmittedExtraordinary.E0[0][1];
					Ec[2] = hit.transmittedOrdinary.E0[0][2] + hit.transmittedExtraordinary.E0[0][2];
					Ec[3] = hit.transmittedOrdinary.E0[0][3] + hit.transmittedExtraordinary.E0[0][3];					
					a[1][idx] += Pol.intensity(E[0]);
					a[2][idx] = FastMath.atan2(Ec[1], Ec[0]);
					a[3][idx] = FastMath.atan2(Ec[3], Ec[2]);
					a[4][idx] = (hit.transmittedExtraordinary.length - hit.transmittedOrdinary.length) / wavelen;
				}	
				Pol.recoverAll();			
			}
			System.out.print(".");
			if((i % 100) == 0)
				System.out.println();
		}
		
		BinaryMatrixFile.mustWrite(outPath + "/imgScan.bin", new double[][]{ scanZ, ang, img, hitC}, true);
		BinaryMatrixFile.mustWrite(outPath + "/angScan.bin", a, true);

		if(vrmlOut != null) {
			vrmlOut.drawOptic(all);
			vrmlOut.destroy();
		}
	}
	
	public void buildPSFGrid() {
		
		PointSpreadBuilder psfBuild = new PointSpreadBuilder(imgPlane, outPath + "/psfData.bin");
		PSFGrid psfGrid = new PSFGrid("imgTest", gridDef, DualGaussianPSF.class);
		psfBuild.setMaxCoherentIntegrationRadius(0.005);
		
		VRMLDrawer vrmlOut = new VRMLDrawer(outPath + "/savartTest.vrml", 0.01);
		if(vrmlOut != null) {
			vrmlOut.setDrawPolarisationFrames(true);
			vrmlOut.setSkipRays(777);
		}
		
		//double col[][] = ColorMaps.jet(nRaysPerSource);
		double col[][] = ColorMaps.jet(psfGrid.getNY()*psfGrid.getNZ());
		
		for(int iX=0; iX < psfGrid.getNX(); iX++) {
			for(int iY=0; iY < psfGrid.getNY(); iY++) {
				for(int iZ=0; iZ < psfGrid.getNZ(); iZ++) {
					double startPos[] = psfGrid.gridPos(iX, iY, iZ);
					
					psfBuild.startNewPSF(startPos,
							//new GaussianPSF()
							new DualGaussianPSF(20)
							//new PointsPSF()
							//new MiniImagePSF(20, 20)
							);
					
					for(int i=0; i < nRaysPerSource; i++) {
						
						RaySegment ray = new RaySegment();
						ray.startPos = startPos;
						ray.dir = Tracer.generateRandomRayTowardSurface(ray.startPos, objLens.getSurfaces().get(0));
						
						ray.length = Double.POSITIVE_INFINITY;
						ray.up = Util.cross(Util.reNorm(Util.cross(ray.dir, new double[]{0,0,1})), ray.dir);
						
						ray.E0 = PointSpreadFunction.getInputStatesForMuellerCalc();
						ray.wavelength = wavelen;
						
						Tracer.trace(all, ray, 100, 0.05, true);
						
						int cIdx = iY*psfGrid.getNZ()+iZ;
						if(vrmlOut != null && (cIdx % 4) == 0)
							vrmlOut.drawRay(ray, col[cIdx]);
						
						psfBuild.nextCoherentSet();
						ray.processIntersections(imgPlane, psfBuild);
						
						Pol.recoverAll();
					}
					
					PointSpreadFunction psf = psfBuild.psfDone(true);
					psfGrid.put(iX, iY, iZ, psf);
				
					System.out.print(".");
				}
			}
			System.out.println("\n" + iX + ", ");
		}
		
		if(vrmlOut != null) {
			vrmlOut.drawOptic(all);
			vrmlOut.destroy();
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
