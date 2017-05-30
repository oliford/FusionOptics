package fusionOptics.polarisation;

import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.List;

import javax.net.ssl.SSLEngineResult.Status;

import fusionOptics.MinervaOpticsSettings;
import fusionOptics.Util;
import fusionOptics.drawing.VRMLDrawer;
import fusionOptics.interfaces.Absorber;
import fusionOptics.interfaces.IsoIsoInterface;
import fusionOptics.interfaces.IsoIsoStdFresnel;
import fusionOptics.interfaces.IsoUniaxialInterface;
import fusionOptics.interfaces.NullInterface;
import fusionOptics.interfaces.Reflector;
import fusionOptics.interfaces.SimplePolariser;
import fusionOptics.materials.BK7;
import fusionOptics.materials.Calcite;
import fusionOptics.materials.CrystalQuartz;
import fusionOptics.materials.IsotropicFixedIndexGlass;
import fusionOptics.materials.LithiumNiobate;
import fusionOptics.materials.SchottSFL6;
import fusionOptics.materials.UniaxialFixedIndexGlass;
import fusionOptics.optics.Box;
import fusionOptics.optics.SimpleDoubleConvexLens;
import fusionOptics.pointSpread.DualGaussianPSF;
import fusionOptics.pointSpread.GaussianPSF;
import fusionOptics.pointSpread.PointSpreadBuilder;
import fusionOptics.pointSpread.PointSpreadFunction;
import fusionOptics.surfaces.Cylinder;
import fusionOptics.surfaces.Disc;
import fusionOptics.surfaces.Dish;
import fusionOptics.surfaces.Iris;
import fusionOptics.surfaces.Square;
import fusionOptics.tracer.Tracer;
import fusionOptics.types.Element;
import fusionOptics.types.Interface;
import fusionOptics.types.Intersection;
import fusionOptics.types.Material;
import fusionOptics.types.Medium;
import fusionOptics.types.Optic;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;
import fusionOptics.types.Surface;

import binaryMatrixFile.BinaryMatrixFile;

import otherSupport.ColorMaps;
import otherSupport.RandomManager;
import otherSupport.StatusOutput;

import oneLiners.OneLiners;
import jafama.FastMath;
import junit.framework.TestCase;
/**
 * Test that the PSF Mueller matrix determination
 * obtains the correct matricies for polarisers
 * and 1/2 and 1/4 wave plates.  
 * 
 * @author oliford
 */
public class MullerTest extends TestCase {
	final static String outPath = MinervaOpticsSettings.getAppsOutputPath() + "/rayTracing/polTests";
	
	final static int maxHits = 1000;
	final static int maxRays = 10000;
	final static double maxTheta = 30 * Math.PI / 180;
	final static double rt2 = Math.sqrt(2);
	
	final static double wavelen = 593e-9;
	
	public void testMuellerSimpleBirefringence() {
		//we're testing waveplates and always firing rays perpendicular, so we can let the iso/iso interface function anyway.
		IsoUniaxialInterface.simpleBirefringence = true;
		doTestMueller();
	}
	
	public void testMuellerGeneralBirefringence() {
		//  THIS CURRENTLY FAILS!!
		// It should work, but I don't know how to fix it.
			
		//We should now be able todo plates properly and let the UniIsoInterface split the rays.
		//The PSFBuilder should do the appropriate coherent addition to correctly get the plate effects
		//However, this currently works for anything except exact normal incidence and/or
		//incidence plane at exactly 45' to optic axis. Odd.
		IsoUniaxialInterface.simpleBirefringence = false;
		doTestMueller();
	}	
	
	private void doTestMueller(){
		Material wavePlateMat = new UniaxialFixedIndexGlass(1.0, 1.001); // keep index difference low to avoid too much fresnel reflection
		
		//Quartz();
		Medium wavePlateMedium = new Medium(wavePlateMat, new double[][]{{ 0, 0, 1 }}, 300);
		
		double oneWaveDiffWidth = Util.calcWaveplateFullWaveThickness(wavePlateMat, wavelen);
		
		HashMap<String, double[][]> expectedMueller = new HashMap<String, double[][]>();
		
		
		
		// **************** Setup -1: Nothing *******************
		Square nothing = new Square("nothing", new double[]{ 2.5, 0, 0 }, new double[]{ -1, 0, 0 }, new double[]{ 0, 1, 0 }, 0.4, 0.4, 
				NullInterface.ideal());
		expectedMueller.put("nothing", new double[][]{{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}});

		// **************** Setup 0: Just some attenuating glass *******************
		Box attenuator = new Box("attenuator", new double[]{ 2.5, 0, 0 }, 2, 0.4, 0.4, new Medium(new IsotropicFixedIndexGlass(1.0, 0.5)), IsoIsoInterface.ideal());
		double a = 0.063; // no idea why 
		expectedMueller.put("attenuator", new double[][]{{a,0,0,0},{0,a,0,0},{0,0,a,0},{0,0,0,a}});
		
		// **************** Setup 1: Vertical Polariser *******************
		Square vertPol = new Square("vertPol", new double[]{ 2.5, 0, 0 }, new double[]{ -1, 0, 0 }, new double[]{ 0, 1, 0 }, 0.4, 0.4, 
							new SimplePolariser(new double[]{ 0, 0, 1 }, 0));
		expectedMueller.put("vertPol", new double[][]{{1,1,0,0},{1,1,0,0},{0,0,0,0},{0,0,0,0}});
		
		
		// **************** Setup 2: Vertical Polariser *******************
		Square horizPol = new Square("horizPol", new double[]{ 2.5, 0, 0 }, new double[]{ -1, 0, 0 }, new double[]{ 0, 1, 0 }, 0.4, 0.4, 
							new SimplePolariser(new double[]{ 0, 1, 0 }, 0));
		expectedMueller.put("horizPol", new double[][]{{1,-1,0,0},{-1,1,0,0},{0,0,0,0},{0,0,0,0}});
		
		// **************** Setup 3: 45' right Polariser *******************
		Square pos45Pol = new Square("pos45Pol", new double[]{ 2.5, 0, 0 }, new double[]{ -1, 0, 0 }, new double[]{ 0, 1, 0 }, 0.4, 0.4, 
							new SimplePolariser(new double[]{ 0, rt2/2, rt2/2 }, 0));
		Util.rotateOnX(pos45Pol, pos45Pol.getBoundarySphereCentre(),+45 * Math.PI / 180); //actually rotate it, this tests the rotation code too
		expectedMueller.put("pos45Pol", new double[][]{{1,0,-1,0},{0,0,0,0},{-1,0,1,0},{0,0,0,0}});	
		
		// **************** Setup 4: -45' right Polariser *******************
		Square neg45Pol = new Square("neg45Pol", new double[]{ 2.5, 0, 0 }, new double[]{ -1, 0, 0 }, new double[]{ 0, 1, 0 }, 0.4, 0.4, 
							new SimplePolariser(new double[]{ 0, -rt2/2, rt2/2 }, 0));
		Util.rotateOnX(neg45Pol, neg45Pol.getBoundarySphereCentre(), -45 * Math.PI / 180);
		expectedMueller.put("neg45Pol", new double[][]{{1,0,1,0},{0,0,0,0},{1,0,1,0},{0,0,0,0}});	
		

		// **************** Setup 5: 1/2 wave plate, vertical *******************		
		Box halfwaveVert = new Box("halfwaveVert", new double[]{ 2.5, 0.0, 0 }, oneWaveDiffWidth/2, 0.4, 0.4, wavePlateMedium.clone(), IsoUniaxialInterface.ideal());
		expectedMueller.put("halfwaveVert", new double[][]{{1,0,0,0},{0,1,0,0},{0,0,-1,0},{0,0,0,-1}});
		
		// **************** Setup 6: 1/2 wave plate, +45' *******************		
		Box halfwavePos45 = new Box("halfwavePos45", new double[]{ 2.5, 0.0, 0 }, oneWaveDiffWidth/2, 0.4, 0.4, wavePlateMedium.clone(), IsoUniaxialInterface.ideal());
		Util.rotateOnX(halfwavePos45, halfwavePos45.getBoundarySphereCentre(), -45 * Math.PI / 180);
		expectedMueller.put("halfwavePos45", new double[][]{{1,0,0,0},{0,-1,0,0},{0,0,1,0},{0,0,0,-1}});
		
		// **************** Setup 7: 1/2 wave plate, vertical *******************		
		Box quartwaveVert = new Box("quartwaveVert", new double[]{ 2.5, 0.0, 0 }, oneWaveDiffWidth/4, 0.4, 0.4, wavePlateMedium.clone(), IsoUniaxialInterface.ideal());
		expectedMueller.put("quartwaveVert", new double[][]{{1,0,0,0},{0,1,0,0},{0,0,0,1},{0,0,-1,0}});
		
		// **************** Setup 8: 1/2 wave plate, +45' *******************		
		Box quartwavePos45 = new Box("quartwavePos45", new double[]{ 2.5, 0.0, 0 }, oneWaveDiffWidth/4, 0.4, 0.4, wavePlateMedium.clone(), IsoUniaxialInterface.ideal());
		Util.rotateOnX(quartwavePos45, quartwavePos45.getBoundarySphereCentre(), +45 * Math.PI / 180);
		expectedMueller.put("quartwavePos45", new double[][]{{1,0,0,0},{0,0,0,-1},{0,0,1,0},{0,1,0,0}});
		
		// **************** Setup 9: Simple short focal length lens *******************
		Material lensMat = new  IsotropicFixedIndexGlass(1.5);
		SimpleDoubleConvexLens lens1 = SimpleDoubleConvexLens.fromFocalLengthAndCentreThickness("lens1", new double[]{ 2.0, 0, 0 }, new double[]{ -1, 0, 0 }, 0.5, 0.15, 0.005, new Medium(lensMat), IsoIsoStdFresnel.ideal(), wavelen);
		SimpleDoubleConvexLens lens2 = SimpleDoubleConvexLens.fromFocalLengthAndCentreThickness("lens2", new double[]{ 4.0, 0, 0 }, new double[]{ -1, 0, 0 }, 0.5, 0.3, 0.005, new Medium(lensMat), IsoIsoStdFresnel.ideal(), wavelen);
		Cylinder escapeShield = new Cylinder("escapeShield", new double[]{ 3.0, 0, 0 },  new double[]{ -1, 0, 0 }, 0.4, 4.0, Absorber.ideal());
		Optic shrtFocalLenses = new Optic("shrtFocalLenses", new Element[]{ lens1, lens2, escapeShield });
		//Util.rotateOnY(lens1, new double[]{ 2.5, 0, 0 }, 14 * Math.PI / 180);
		//Util.rotateOnY(lens2, new double[]{ 2.5, 0, 0 }, -14 * Math.PI / 180);
		
		// **************** Setup 10: The short focal length lens from the end of the AUG MSE system *******************
		Medium lensGlass = new Medium(new SchottSFL6());
		Interface lensIFace = IsoIsoStdFresnel.ideal();
		
		Iris lens3AIris = new Iris("lens3AIris", new double[]{ 1.5, 0,0 }, new double[]{ -1, 0, 0 }, 0.15, 0.040, Absorber.ideal());
		Dish lens3AFront = new Dish("lens3AFront", new double[]{ 1.52774, 0,0 }, new double[]{ 1, 0, 0 }, 0.128, 0.040, lensGlass, null, lensIFace);
		Disc lens3ABack = new Disc("lens3ABack", new double[]{ 1.535, 0,0 }, new double[]{ 1, 0, 0}, 0.040, null, lensGlass, lensIFace);
		
		Dish lens3BFront = new Dish("lens3BFront", new double[]{ 1.53596, 0,0 }, new double[]{ 1, 0, 0 }, 0.064, 0.037, lensGlass, null, lensIFace);
		Dish lens3BBack = new Dish("lens3BBack", new double[]{ 1.54318, 0,0 }, new double[]{ 1, 0, 0 }, 0.128, 0.037, null, lensGlass, lensIFace);

		Material reparaGlass = new IsotropicFixedIndexGlass(2.8);
		SimpleDoubleConvexLens reparaLens1 = SimpleDoubleConvexLens.fromFocalLengthAndCentreThickness("reparaLens", new double[]{ 1.612 + 0.18, 0, 0}, new double[]{ -1, 0, 0}, 0.4, 0.12, 0.005, new Medium(reparaGlass), IsoIsoInterface.ideal(), wavelen);
		SimpleDoubleConvexLens reparaLens2 = SimpleDoubleConvexLens.fromFocalLengthAndCentreThickness("reparaLens", new double[]{ 1.612 + 0.22, 0, 0}, new double[]{ -1, 0, 0}, 0.4, 0.12, 0.005, new Medium(reparaGlass), IsoIsoInterface.ideal(), wavelen);
		//Iris reparaIris = new Iris("reparaIris", new double[]{ 1.612 + 0.15 }, new double[]{ -1, 0, 0 }, 0.05, 0.04, iface)
		
		
		Optic augShrtLenses = new Optic("augShrtLenses", new Element[]{ lens3AFront, lens3ABack, lens3BFront, lens3BBack,reparaLens1,reparaLens2 });
		augShrtLenses.shift(new double[]{ -1.535 + 2.5, 0, 0 });
		
		
		
		
		// **************** Setup n-1: Depolarisation by lots of reflections down a long tube *******************
		//Square pol1 = new Square("pol1", new double[]{ 0.1, 0, 0 }, new double[]{ -1, 0, 0 }, new double[]{ 0, 1, 0 }, 0.4, 0.4, 
		//		new SimplePolariser(new double[]{ 0, 1/rt2, -1/rt2 }, 0));

		
		//iris to stop light skipping the whole device
		Iris depolIris = new Iris("depolIris", new double[]{ 0.1, 0, 0}, new double[]{ -1, 0, 0}, 0.2, 0.09, Absorber.ideal());
		
		//50% V pol and 50% H pol
		Square depolPol1 = new Square("depolPol1", new double[]{ 0.04, 0.15/2, 0 }, new double[]{ -1, 0, 0 }, new double[]{ 0, 0, 1 }, 0.3, 0.15, 
				new SimplePolariser(new double[]{ 0, 0, 1 }, 0));
		Square depolPol2 = new Square("depolPol2", new double[]{ 0.03, -0.15/2, 0 }, new double[]{ -1, 0, 0 }, new double[]{ 0, 0, 1 }, 0.3, 0.15, 
				new SimplePolariser(new double[]{ 0, 1, 0 }, 0));

		//long tube for light to bounce down
		Cylinder depolCyld = new Cylinder("depolCyld", new double[]{ 2.5, 0, 0 }, new double[]{ 1, 0, 0}, 0.1, 4.8, Reflector.ideal());
		
		//now some randomly slanted divergent lenses along the tube
		int nDepolLenes = 4; 
		
		Element depolElems[] = new Element[4 + nDepolLenes];
		depolElems[0] = depolPol1;
		depolElems[1] = depolPol2;
		depolElems[2] = depolIris;
		depolElems[3] = depolCyld;
		
		for(int i=0; i < nDepolLenes; i++) {
			double c[] = new double[]{ 0.15 + i * 4.8 / nDepolLenes, 0, 0 };
			depolElems[4+i] = SimpleDoubleConvexLens.fromRadiusOfCurvAndCentreThickness("depolLens"+i, c, new double[]{ -1, 0, 0 }, 
					0.12, //radius
					0.5, //radius of curvature
					-0.010, //thickness in centre
					new Medium(new IsotropicFixedIndexGlass(2.5)), IsoIsoInterface.ideal());
			Util.rotateOnY(depolElems[4+i], c, (RandomManager.instance().nextUniform(0, 1) * 80 - 40) * Math.PI / 180);
			Util.rotateOnZ(depolElems[4+i], c, (RandomManager.instance().nextUniform(0, 1) * 80 - 40) * Math.PI / 180);
		}
		
		Optic depol = new Optic("depol", depolElems);
		
		
		// ***********************************************************************************
		
		// The image plane normal needs to be in the ray travelling direction, so that the final polarisation has the same sense
		// as the initial.
		Square imgPlane = new Square("fwdPlane", new double[]{ 5, 6, 0 }, new double[]{ 1, 0, 0 }, new double[]{ 0, 0, 1 }, 0.5, 13, Absorber.ideal());		
				
		Element testOptics[] = new Element[]{ 
				imgPlane,
				nothing,
				attenuator, 
				vertPol, horizPol, 
				pos45Pol, neg45Pol, 
				halfwaveVert,
				halfwavePos45,
				quartwaveVert, 
				quartwavePos45,
				shrtFocalLenses,
				augShrtLenses,
				//depol,
			};
		
		int nTests = testOptics.length - 1;
		
		for(int i=0; i < nTests; i++){
			testOptics[i+1].shift(new double[]{ 0, i, 0 });			
		}
		
		Optic all = new Optic("all", testOptics);
		
		VRMLDrawer vrmlOut = new VRMLDrawer(outPath + "/muller.vrml", 0.01);
		vrmlOut.setDrawPolarisationFrames(true);
		vrmlOut.setSkipRays(maxHits / 27);
		
		double col[][] = ColorMaps.jet(maxRays);
		
		boolean allOK = true;
		
		PointSpreadBuilder psfBuild = new PointSpreadBuilder(imgPlane);
		psfBuild.setMaxCoherentIntegrationRadius(0.001);
		
		StatusOutput stat = new StatusOutput(MullerTest.class, nTests*maxHits);
		for(int iT=0; iT < nTests; iT++){
			psfBuild.startNewPSF(new double[]{ 0, 0, 0 }, new DualGaussianPSF());
			for(int i=0; i < maxRays; i++) {
					
				double c[] = testOptics[iT+1].getBoundarySphereCentre();
				double minR = 0.0, maxR = 0.040;
				
				RaySegment ray = new RaySegment();
				
				double r = FastMath.sqrt(minR*minR + (maxR * maxR - minR * minR) * RandomManager.instance().nextUniform(0,1));
				double phi = RandomManager.instance().nextUniform(0, 1) * 2 * Math.PI;
				ray.startPos = new double[]{ 0, 
							c[1] + r * FastMath.cos(phi),
							c[2] + r * FastMath.sin(phi)
						};
				
				//roughly perp, with a little bit of randomness
				ray.dir = Util.reNorm(new double[]{ 1, 
						(RandomManager.instance().nextUniform(0,1)*0.05)-0.025, 
						(RandomManager.instance().nextUniform(0,1)*0.05)-0.025 });
				
				//but...
				if(IsoUniaxialInterface.simpleBirefringence && testOptics[iT+1] instanceof Box){
					List<Medium> boxMedia = ((Box)testOptics[iT+1]).getUniqueMediaList();
					if(boxMedia.size() > 0 && boxMedia.get(0).getOpticAxes() != null && boxMedia.get(0).getOpticAxes().length > 0){
						//testing a uniaxial medium in simple mode, the ray must be exactly perpendicular
						ray.dir = new double[]{ 1, 0, 0 };
					}
				}
				
				ray.length = Double.POSITIVE_INFINITY;
				ray.up = Util.cross(Util.reNorm(Util.cross(ray.dir, new double[]{0,0,1})), ray.dir);
				
				ray.E0 = PointSpreadFunction.getInputStatesForMuellerCalc();
				//ray.E0 = new double[][]{ { 1, 0, 1, 0 }, {0,0,0,0},{0,0,0,0},{0,0,0,0} };
				
				ray.wavelength = wavelen;
				
				Tracer.trace(all, ray, 1000, 0.1, true);
				if(i==0 && false){
					ray.dumpPath();
					//RaySegment O = ray.endHit.transmittedOrdinary;
					//RaySegment E = ray.endHit.transmittedExtraordinary;
					RaySegment O = ray.endHit.transmittedOrdinary;
					RaySegment E = ray.endHit.transmittedExtraordinary;
					O.rotatePolRefFrame(new double[]{ 0, 0, 1 });
					E.rotatePolRefFrame(new double[]{ 0, 0, 1 });
					double OE[] = O.E0[2], EE[] = E.E0[2];
					OneLiners.dumpArray(OE);
					OneLiners.dumpArray(EE);
					System.out.println(FastMath.atan2(OE[Pol.uRe], OE[Pol.uIm])*180/Math.PI);
					System.out.println(FastMath.atan2(EE[Pol.rRe], EE[Pol.rIm])*180/Math.PI);
					System.out.println(Pol.psi(OE)*180/Math.PI);
					System.out.println(Pol.psi(EE)*180/Math.PI);
					
				}
				
				//ray.dumpPath();
				//System.exit(0);
				vrmlOut.drawRay(ray, col[i]);
				
				ray.processIntersections(imgPlane, psfBuild);
				psfBuild.nextCoherentSet();
				
				
				
				int n = psfBuild.getNPointsCollected();
				stat.doStatus(iT*maxHits + n);
				if(n >= maxHits)
					break;
			}

			DualGaussianPSF psf = (DualGaussianPSF) psfBuild.psfDone(false);
			System.out.println(testOptics[iT+1].getName() + ":");
			//if(testOptics[iT+1] == depol)
			//	psf.dumpNormalisedMueller();
			//else
			//	psf.dumpMueller();

			Pol.recoverAll();
			

			double M[] = psf.getMeanMueller();
			double Me[][] = expectedMueller.get(testOptics[iT+1].getName());
			
			DecimalFormat fmt = new DecimalFormat("0.00");
			for(int i=0; i < 4; i++){
				System.out.print("got[  ");
				for(int j=0; j < 4; j++)
					System.out.print(fmt.format(M[i*4+j])+"  ");
				
				System.out.print("  ]\texpect(");
				if(Me == null){
					System.out.print(" -- No expected Meuller --");
				}else{
					boolean lineOK = true;
					for(int j=0; j < 4; j++){
						System.out.print(fmt.format(Me[i][j])+"  ");
						if(FastMath.abs(Me[i][j] - M[i*4+j]) > 0.01)
							lineOK = false;
					}
					System.out.print(")");
					
					if(!lineOK){
						System.out.print("  <-- ** MISTMATCH **");
						allOK = false;
					}
				}
				
				System.out.print("\n");
				
			}
			
		}
		stat.done();
		

		vrmlOut.drawOptic(all);
		vrmlOut.destroy();

		assertTrue("Meuller calculation mistmatched expected. See stdout for details.", allOK);

	}
	
	
}
