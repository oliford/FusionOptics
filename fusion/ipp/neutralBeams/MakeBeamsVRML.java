package ipp.neutralBeams;

import otherSupport.ColorMaps;
import fusionOptics.Util;
import fusionOptics.drawing.VRMLDrawer;
import fusionOptics.types.Optic;
import fusionOptics.types.Surface;
import jafama.FastMath;

public class MakeBeamsVRML {
	private static final double dL = 0.05;

	private static final double[][] cols = ColorMaps.jet(20);
	
	public static double vrmlTransformMatrix[][] = new double[][]{ {1000,0,0},{0,1000,0},{0,0,1000}};
	
	public static void makeBeamsAllGreen(String fileName, SimpleBeamGeometry beamGeom){
			
		VRMLDrawer vrmlOut = new VRMLDrawer(fileName, 0.005);
		vrmlOut.setTransformationMatrix(vrmlTransformMatrix);
		vrmlOut.setDrawPolarisationFrames(false);
		vrmlOut.drawOptic(beamGeom.makeAllBeamCylds());
		vrmlOut.destroy();
	}
	
	public static void makeBeamsPINIColoured(String fileName, SimpleBeamGeometry beamGeom){
		
		VRMLDrawer vrmlOut = new VRMLDrawer(fileName, 0.005);
		vrmlOut.setTransformationMatrix(vrmlTransformMatrix);
		vrmlOut.setDrawPolarisationFrames(false);
		
		double cols[][] = new double[][]{ 
				{1,1,0}, {0,1,0}, {1,0,1}, {1,0,0},
				{1,1,0}, {0,1,0}, {1,0,1}, {1,0,0},
			};
		
		for(int i=0; i < 8; i++){
			vrmlOut.addMat("beam_Q"+(i+1), cols[i][0] + " " + cols[i][1] + " " + cols[i][2], 0.7);
		}
				
		for(int i=0; i < beamGeom.nBeams(); i++){
			double l0 = beamGeom.getLOfBeamAxisAtR(i, beamGeom.plasmaR1);
			double l1 = beamGeom.getLOfBeamAxisAtR(i, beamGeom.plasmaR0);
			if(Double.isNaN(l1)){ l1 = 5; } //some beams dont hit the inner wall
			vrmlOut.startGroup("beamQ" +(i+1));
			vrmlOut.drawOptic(beamGeom.makeBeamCylds(i, dL, beamGeom.beamWidth(), l0, l1), "beam_Q"+(i+1), cols[i]);
			vrmlOut.endGroup();
		}
		
		vrmlOut.destroy();
	}
	
	public static void makeBeamsPINIColouredSeparate(String fileNamePrefix, SimpleBeamGeometry beamGeom){
		
		
		double cols[][] = new double[][]{ 
				{1,1,0}, {0,1,0}, {1,0,1}, {1,0,0},
				{1,1,0}, {0,1,0}, {1,0,1}, {1,0,0},
			};
		
				
		for(int i=0; i < beamGeom.nBeams(); i++){
			VRMLDrawer vrmlOut = new VRMLDrawer(fileNamePrefix + "-Q" + (i+1) + ".vrml", 0.005);
			vrmlOut.setTransformationMatrix(vrmlTransformMatrix);
			vrmlOut.setDrawPolarisationFrames(false);
			
			vrmlOut.addMat("beam_Q"+(i+1), cols[i][0] + " " + cols[i][1] + " " + cols[i][2], 0.7);
						
			double l0 = beamGeom.getLOfBeamAxisAtR(i, beamGeom.plasmaR1);
			double l1 = beamGeom.getLOfBeamAxisAtR(i, beamGeom.plasmaR0);
			if(Double.isNaN(l1)){ l1 = 5; } //some beams dont hit the inner wall
			vrmlOut.startGroup("beamQ" +(i+1));
			vrmlOut.drawOptic(beamGeom.makeBeamCylds(i, dL, beamGeom.beamWidth(), l0, l1), "beam_Q"+(i+1), cols[i]);
			vrmlOut.endGroup();
			
			vrmlOut.destroy();
		}
		
	}	
	
	public static void makeBeamsRadialColoured(String fileName, SimpleBeamGeometry beamGeom){
		VRMLDrawer vrmlOut = new VRMLDrawer(fileName, 0.005);
		for(int i=0; i < cols.length; i++)
			vrmlOut.addMat("beamR" + i, 
					cols[i][0] + " " + cols[i][1] + " " + cols[i][2],
					0.3);
		
		vrmlOut.addMat("beamOutOfRange", "0 0 0" , 0.8);
		vrmlOut.setTransformationMatrix(vrmlTransformMatrix);
		vrmlOut.setDrawPolarisationFrames(false);
		drawRadialColoured(vrmlOut, beamGeom.makeAllBeamCylds(), beamGeom);
		vrmlOut.destroy();
	}
	
	public static void drawRadialColoured(VRMLDrawer vrmlOut, Optic optic, SimpleBeamGeometry beamGeom){
		
		for(Optic subOptic : optic.getSubOptics()){
			drawRadialColoured(vrmlOut, subOptic, beamGeom);
		}
		
		for(Surface surf : optic.getSurfaces()){
			double pos[] = surf.getBoundarySphereCentre();
			double R = FastMath.sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
			
			int colIdx = (int)(cols.length * (R - beamGeom.plasmaR0) / (beamGeom.plasmaR1 - beamGeom.plasmaR0));
			
			if(colIdx >= 0 && colIdx < cols.length)
				vrmlOut.drawSurface(surf, "beamR" + colIdx, cols[colIdx]);
			//else
				//vrmlOut.drawSurface(surf, "beamOutOfRange");
			
		}
	}
}
