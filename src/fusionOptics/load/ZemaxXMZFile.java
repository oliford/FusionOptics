package fusionOptics.load;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;

import fusionOptics.interfaces.Absorber;
import fusionOptics.interfaces.IsoIsoInterface;
import fusionOptics.interfaces.IsoIsoStdFresnel;
import fusionOptics.materials.IsotropicLinearDispersiveGlass;
import fusionOptics.surfaces.Aspheric;
import fusionOptics.surfaces.Disc;
import fusionOptics.surfaces.Dish;
import fusionOptics.surfaces.Iris;
import fusionOptics.types.Material;
import fusionOptics.types.Medium;
import fusionOptics.types.Optic;
import fusionOptics.types.Surface;

import otherSupport.SimpleScriptParser;

import oneLiners.OneLiners;


/** Really basic loading of Zemax ZMX file 
 * Construct class instance with file name, then call load() 
 * 
 */
public class ZemaxXMZFile {
	
	private String fileName;
	
	private Optic mainOptic;
	
	private SurfaceInfo nextSurface;
	
	private Material currentMat;
	private Medium currentMedium;
	
	private double currentDist = 0;
	
	private double irisOuterRadius = Double.NEGATIVE_INFINITY;
	
	private class SurfaceInfo {
		public int number = -1 ;
		public String type;
		public double curv[];
		public String curvStr;
		public int hide[];
		public int mirr[];
		public int slab;
		public double diam[];
		public String diamStr;
		public double distToNext = Double.NaN;
		public Material materialToNext;
		public int pops[];
		public double conicConst;
		public HashMap<Integer,Double> params;
	}	

	/** @return An Optic object containing all the surfaces described by the ZMX file */ 
	public Optic getMainOptic() { return mainOptic; }
	
	public ZemaxXMZFile(String fileName) {
		this.fileName = fileName;
	}
	
	/** Loads the ZMX file encapsulated by this object. The
	 * created object can be retrieved by calling getMainOptic()  */
	public void load() {
		this.mainOptic = new Optic("LoadedZMX");
		currentDist = 0;
		currentMat = null;
		currentMedium = null;
		
		SimpleScriptParser.parse(this, fileName, "cmd", true);
	}
	
	/** Must be called before load() */
	public void setIrisOuterRadius(double irisOuterRadius){
		this.irisOuterRadius = irisOuterRadius;
	}
	
	public void cmdVERS(int a, int b, int c){
		System.out.println("ver"+a+","+b+","+c);
	}

	public void cmdNAME(String name){
		mainOptic.setName(name);
		System.out.println(name);
	}
	
	public void cmdSURF(int surfaceNum){
		if(nextSurface != null)
			surfaceDone();
		nextSurface = new SurfaceInfo();
		nextSurface.number = surfaceNum;
	}
	
	public void cmdTYPE(String surfType){ nextSurface.type = surfType;	}

	public void cmdCURV(double a, double b, double c, double d, double e, String str){
		nextSurface.curv = new double[]{ a, b, c, d, e };
		nextSurface.curvStr = str;
	}
	
	public void cmdHIDE(int hideData[]){ nextSurface.hide = hideData; }
	
	public void cmdMIRR(int mirrData[]){ nextSurface.mirr = mirrData; }
	
	public void cmdSLAB(int slab){ nextSurface.slab = slab; }
	
	public void cmdDISZ(double disz){ nextSurface.distToNext = disz; }
	
	public void cmdDIAM(double a, double b, double c, double d, double e, String str){
		nextSurface.diam = new double[]{ a, b, c, d, e };
		nextSurface.diamStr = str;
	}
	
	public void cmdPOPS(int popsData[]){ nextSurface.pops = popsData; }
	
	public void cmdGLAS(String name, int solveType, double unknown, double index_d, double abbe_d, int a, int b, int c, int d, int e, int f){
		nextSurface.materialToNext = new IsotropicLinearDispersiveGlass(index_d, abbe_d);		
	}
	
	public void cmdCONI(double conicConst){ nextSurface.conicConst = conicConst; }
	
	public void cmdPARM(int index, double val){
		if(nextSurface.params == null)
			nextSurface.params = new HashMap<Integer, Double>();
		nextSurface.params.put(index,val);		
	}
	
	private double scaleXY = 0.001;
	private double scaleZ = 0.001;
	private void surfaceDone(){
		double centre[] = new double[]{ currentDist, 0, 0 };

		Medium nextMedium = (nextSurface.materialToNext == null) 
					? null : new Medium(nextSurface.materialToNext);
		
		double irisRadius = -1;
		
		if("STANDARD".equalsIgnoreCase(nextSurface.type)){
			if(nextSurface.diam[0] <= 0){
				System.out.println("Invalid surface "+nextSurface.number+": diam < -0");
				return;
			}
			if(nextSurface.number > 1000)
				return;
			
			if(nextSurface.curv[0] != 0){
				double dishCentreNormal[];
				Medium frontMedium, backMedium;
				if(nextSurface.curv[0] > 0){
					dishCentreNormal = new double[]{ 1, 0, 0 };	
					frontMedium = nextMedium;
					backMedium = currentMedium;
				}else{
					nextSurface.curv[0] *= -1;
					dishCentreNormal = new double[]{ -1, 0, 0 };
					frontMedium = currentMedium;
					backMedium = nextMedium;				
				}
				
				Dish dish = new Dish("surf"+nextSurface.number+"-dish", 
						centre,
						dishCentreNormal, 
						(1/nextSurface.curv[0]) * scaleXY, 
						nextSurface.diam[0] * scaleXY,
						frontMedium, 
						backMedium,
						IsoIsoInterface.ideal());
				
				irisRadius = nextSurface.diam[0] * scaleXY;
				mainOptic.addElement(dish);
				
			}else{
				Disc disc = new Disc("surf"+nextSurface.number+"-disc", 
						centre, 
						new double[]{ 1, 0, 0 }, //normal 
						nextSurface.diam[0] * scaleXY,
						nextMedium,
						currentMedium,
						IsoIsoInterface.ideal());
				
				irisRadius = nextSurface.diam[0] * scaleXY;				
				mainOptic.addElement(disc);
			}
		}else if("EVENASPH".equalsIgnoreCase(nextSurface.type)){
			double dishCentreNormal[];
			Medium frontMedium, backMedium;
			if(nextSurface.curv[0] > 0){
				dishCentreNormal = new double[]{ 1, 0, 0 };	
				frontMedium = nextMedium;
				backMedium = currentMedium;
			}else{
				nextSurface.curv[0] *= -1;
				dishCentreNormal = new double[]{ -1, 0, 0 };
				frontMedium = currentMedium;
				backMedium = nextMedium;				
			}
			
			double polyCoeffs[];
			if(nextSurface.params == null){
				System.err.println("WARNING: No aspheric coefficients found on even aspheric surface");
				polyCoeffs = new double[0];
			}else{
					
				
				//find highest order
				int maxHalfOrder = 0;
				for(Entry<Integer, Double> entry : nextSurface.params.entrySet()){
					int halfOrder = entry.getKey();
					if(entry.getValue() != 0 && halfOrder > maxHalfOrder)
						maxHalfOrder = halfOrder;
				}
				
				polyCoeffs = new double[maxHalfOrder];
				for(Entry<Integer, Double> entry : nextSurface.params.entrySet()){
					double val = entry.getValue();
					if(val != 0)
						polyCoeffs[entry.getKey() - 1] = val;
				}
			}
				
			Aspheric asphere = new Aspheric("aspheric"+nextSurface.number+"-dish", 
					centre,
					dishCentreNormal, 
					(1/nextSurface.curv[0]) * scaleXY, 
					nextSurface.diam[0] * scaleXY,
					nextSurface.conicConst,
					polyCoeffs,
					frontMedium, 
					backMedium,
					IsoIsoInterface.ideal());
			
			irisRadius = nextSurface.diam[0] * scaleXY;
			mainOptic.addElement(asphere);
			
			
		}
		
		if(irisOuterRadius > 0 && irisRadius > 0){
			mainOptic.addElement(new Iris(
					"iris-surf" + nextSurface.number,
					centre.clone(),
					new double[]{ 1, 0, 0 },
					irisOuterRadius,
					irisRadius,
					Absorber.ideal()
					));
			
		}
		
		currentMedium = nextMedium;
		if(!Double.isInfinite(nextSurface.distToNext) && !Double.isNaN(nextSurface.distToNext))
			currentDist += nextSurface.distToNext * scaleZ;
	}
		

}
