package fusionOptics.optics;

import fusionOptics.interfaces.Absorber;
import fusionOptics.interfaces.IsoIsoInterface;
import fusionOptics.materials.IsotropicLinearDispersiveGlass;
import fusionOptics.surfaces.Disc;
import fusionOptics.surfaces.Dish;
import fusionOptics.surfaces.Iris;
import fusionOptics.types.Material;
import fusionOptics.types.Medium;
import fusionOptics.types.Optic;
import fusionOptics.types.Surface;
import jafama.FastMath;

/** Simple sequential lens series from data */
public class SequentialLensSeries extends Optic {
	
	public SequentialLensSeries(String name, double data[][], double scale, double irisOuterRadius) {
		this(data, scale, irisOuterRadius);
		this.name = name;
	}
		
	/** each entry of data[][] is:
	 * radius.curv,	distance to next,	index_d, abbe_d 
	 * 
	 * @param data
	 */
	public SequentialLensSeries(double data[][], double scale, double irisOuterRadius) {
		super("SequentialLensSeries");
		double currentDist = 0;
		Medium currentMedium = null;
		
		for(int i=0; i < data.length; i++){
			Surface surf;
			double radCurv = data[i][0] * scale;
			double rimRadius = data[i][2] * scale;
			double nd = data[i][3];
			double abbe_d = data[i][4];
			double centre[] = new double[]{ currentDist, 0, 0 };
					
			if(rimRadius > 0){
				
				Material nextMat = null; 
				Medium nextMedium = null;
				if(nd != 0){
					nextMat = new IsotropicLinearDispersiveGlass(nd, abbe_d);
					nextMedium = new Medium(nextMat);
				}
				
				if(radCurv != 0){
					surf = new Dish("dish"+i, 
							centre,
						new double[]{ (radCurv > 0) ? 1 : -1, 0, 0 },
						FastMath.abs(radCurv),
						rimRadius,
						(radCurv > 0) ? nextMedium : currentMedium,
						(radCurv > 0) ? currentMedium : nextMedium,
						IsoIsoInterface.ideal());
				}else{
					surf = new Disc("disc"+i, 
							centre,
							new double[]{ (radCurv > 0) ? 1 : -1, 0, 0 },
							rimRadius,
							(radCurv > 0) ? nextMedium : currentMedium,
							(radCurv > 0) ? currentMedium : nextMedium,
							IsoIsoInterface.ideal());
				}
				addElement(surf);
				
				if(irisOuterRadius > 0){
					addElement(new Iris(
							"iris-surf" + i,
							centre.clone(),
							new double[]{ 1, 0, 0 },
							irisOuterRadius,
							0.99*rimRadius,
							Absorber.ideal()
							));
					
				}
				currentMedium = nextMedium;
			}
			
			currentDist += data[i][1] * scale;
		}
	}

}
