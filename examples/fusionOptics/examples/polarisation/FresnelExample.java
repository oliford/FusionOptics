package fusionOptics.polarisation;

import java.util.Arrays;
import java.util.List;

import fusionOptics.MinervaOpticsSettings;
import fusionOptics.Util;
import fusionOptics.drawing.VRMLDrawer;
import fusionOptics.interfaces.Absorber;
import fusionOptics.interfaces.IsoIsoAntiReflective;
import fusionOptics.interfaces.IsoIsoInterface;
import fusionOptics.interfaces.IsoIsoStdFresnel;
import fusionOptics.materials.BK7;
import fusionOptics.materials.Calcite;
import fusionOptics.materials.CrystalQuartz;
import fusionOptics.materials.IsotropicFixedIndexGlass;
import fusionOptics.materials.SchottSFL6;
import fusionOptics.surfaces.Iris;
import fusionOptics.surfaces.Square;
import fusionOptics.tracer.Tracer;
import fusionOptics.types.Element;
import fusionOptics.types.Intersection;
import fusionOptics.types.Medium;
import fusionOptics.types.Optic;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;
import fusionOptics.types.Surface;

import otherSupport.ColorMaps;

import binaryMatrixFile.BinaryMatrixFile;

import oneLiners.OneLiners;
import net.jafama.FastMath;

/** Examine transmittance and reflectance vs angle and polarisation for a single interface */
public class FresnelExample {
	final static String outPath = MinervaOpticsSettings.getAppsOutputPath() + "/rayTracing/fresnelTest";
	
	final static double minAngle = 0 * Math.PI / 180;
	final static double maxAngle = 89.5 * Math.PI / 180;
	final static int nAngles = 100;
	final static double rt2 = Math.sqrt(2);
	final static double wavelength = 600e-9;
	
	public static void main(String[] args) {
		
		Medium glass = new Medium(new IsotropicFixedIndexGlass(1.8));
		//Medium glass = new Medium(new SchottSFL6());
		
		Square backPlane = new Square("backPlane", new double[]{ 0, 0, 0 }, new double[]{ rt2/2, -rt2/2, 0 }, new double[]{ rt2/2, rt2/2, 0 }, 2*rt2, 1.0, null, null, Absorber.ideal());
		Square fwdPlane = new Square("fwdPlane", new double[]{ 2, 0, 0 }, new double[]{ -rt2/2, -rt2/2, 0 }, new double[]{ rt2/2, -rt2/2, 0 }, 2*rt2, 1.0, null, null, Absorber.ideal());
		Square glassSurface = new Square("iface", new double[]{ 1, 0, 0 }, new double[]{ -1, 0, 0 }, new double[]{ 0, 1, 0 }, 1.0, 1.0, 
				null, glass, 
				IsoIsoStdFresnel.ideal());
				//new IsoIsoAntiReflective(1.38, wavelength/4/1.38, 0.0));				
			
		Optic all = new Optic("all", new Element[]{ backPlane, glassSurface, fwdPlane });
		
		VRMLDrawer vrmlOut = new VRMLDrawer(outPath + "/fresnel.vrml", 0.01);
		vrmlOut.setDrawPolarisationFrames(true);
		vrmlOut.setSkipRays(10);
		
		double frontIndex = (glassSurface.getFrontMedium() == null) ? 1.0 : glassSurface.getFrontMedium().getRefractiveIndex(0, wavelength);
		double backIndex = (glassSurface.getBackMedium() == null) ? 1.0 : glassSurface.getBackMedium().getRefractiveIndex(0, wavelength);
		
		double brewsterAng = Math.atan(backIndex / frontIndex);
		
		double ang[] = OneLiners.linSpace(minAngle, maxAngle, nAngles-1);
		ang = Arrays.copyOf(ang, nAngles);
		ang[nAngles-1] = brewsterAng;
		
		double targY[] = OneLiners.linSpace(-0.48, 0.48, nAngles);
		double targZ[] = OneLiners.linSpace(-0.43, 0.43, nAngles);
		double sLen = 0.5;
		
		double col[][] = ColorMaps.jet(nAngles);
		
		double out[][] = new double[nAngles][7];
			
		
		for(int i=0; i < nAngles; i++) {
			
			RaySegment ray = new RaySegment();
			ray.startPos = new double[]{
					1.0- sLen * Math.cos(ang[i]),
					targY[i] - sLen * Math.sin(ang[i]),
					targZ[i]
			};
			ray.dir = Util.reNorm(new double[]{
					sLen * Math.cos(ang[i]),
					sLen * Math.sin(ang[i]),
					0,
			});
			ray.length = Double.POSITIVE_INFINITY;
			ray.up = new double[]{ 0, 0, 1 };
			
			ray.E0 = new double[][]{ {1,0,1,0}, {1,0,-1,0} };
	
			
			ray.wavelength = wavelength;
			
			double initPsi = Pol.psi(ray.E0[0]);
			
			Tracer.trace(all, ray, 10, 1e-10, true);
			
			double tPsi=Double.NaN, rPsi=Double.NaN;
			out[i][0] = ang[i];
			
			List<Intersection> backHits = ray.getIntersections(backPlane);
			if(backHits.size() > 0){
				RaySegment rRay = backHits.get(0).incidentRay;
				rRay.rotatePolRefFrame(ray.up); //put polarisation back in initial frame
				out[i][1] = Pol.mag2Eu(rRay.E1[0]); //E0_up^2 (Rs)
				out[i][2] = Pol.mag2Er(rRay.E1[0]); //E0_right^2 (Rp)
				rPsi = Pol.psi(rRay.E1[0]);
				out[i][3] = rPsi;
				if(i < (nAngles-1)){
					if(ang[i] < brewsterAng) {
						if((rPsi * initPsi) > 0)				
							System.err.println("Expecting a change of polarisation angle sign, for reflected with theta < brewster");
					} else {
						if((rPsi * initPsi) < 0)				
							System.err.println("Expecting the same polarisation angle sign, for reflected with theta > brewster");
					}
				}
					
			}
			List<Intersection> frontHits = ray.getIntersections(fwdPlane);
			if(frontHits.size() > 0){
				RaySegment tRay = frontHits.get(0).incidentRay;
				tRay.rotatePolRefFrame(ray.up); //put polarisation back in initial frame
				out[i][4] = Pol.mag2Eu(tRay.E1[0]); //E0_up^2 (Ts)
				out[i][5] = Pol.mag2Er(tRay.E1[0]); //E0_right^2 (Tp)
				tPsi = Pol.psi(tRay.E1[0]);
				out[i][6] = tPsi;
				if(i < (nAngles-1) && ang[i] < brewsterAng && (rPsi * initPsi) > 0)
					System.err.println("Expecting the same polarisation angle sign, for transmissed rays ");
			}			
			
			System.out.println(ang[i]*180/Math.PI + "\t" + Pol.psi(ray.E0[0])*180/Math.PI + "\t" + rPsi*180/Math.PI + "\t" + tPsi*180/Math.PI);
					
			//What to expect:
			//Angles here are clockwise from 'up' looking in the ray travel direction
			// Starting with +45' ( s = [ 1 0 1 0 ] )
			// 
			// Transmitted should ALWAYS be the same - none of the transmission ampltiude coefficients are -ve 
			// 
			// For reflected:
			//   ni < nt (glass entry), ang < brewster: 			Rs=-ve, Rp=+ve ==> Flipped -45'
			//   ni < nt (glass entry), ang > brewster: 			Rs=-ve, Rp=-ve ==> Same +45'
			//   ni > nt (glass exit), ang < brewster: 				Rs=+ve, Rp=-ve ==> Flipped -45'
			//   ni > nt (glass exit), brewster < ang < critical: 	Rs=+ve, Rp=+ve ==> Same +45'
			//
			
			vrmlOut.drawRay(ray, col[i]);
		}
		vrmlOut.drawOptic(all);
		vrmlOut.destroy();
		
		BinaryMatrixFile.mustWrite(outPath + "/stokes.bin", out, false);
	}
	
	
}
