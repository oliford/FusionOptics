package fusionOptics.interfaces;

import java.util.ArrayList;
import java.util.List;

import fusionOptics.Util;
import fusionOptics.types.Interface;
import fusionOptics.types.Intersection;
import fusionOptics.types.Medium;


/** Interfaces that can handle two media (although might not have either) */
public abstract class DualMediumInterface implements Interface {
	public boolean ignoreInterfaceCompatibility = false;
	
	public abstract void calcIntersection(Intersection hit, Medium incidentMedium, Medium transmissionMedium, double minIntensity);

	/** Calculates the ongoing rays etc at an intersection */
	public final void calcIntersection(Intersection hit, double minIntensity) {
		
		//first we need to see which side we are hitting the surface from
		double uDotN = Util.dot(hit.incidentRay.dir, hit.normal);
		
		Medium incidentMedium, transmissionMedium;
		if(uDotN <= 0){ //antiparallel, we are approaching from the surface 'front'
			incidentMedium = hit.surface.getFrontMedium();
			transmissionMedium = hit.surface.getBackMedium();
		}else{ //parallel, we are hitting from the back
			incidentMedium = hit.surface.getBackMedium();
			transmissionMedium = hit.surface.getFrontMedium();
			hit.normal[0] *= -1;
			hit.normal[1] *= -1;
			hit.normal[2] *= -1;
		}
		
		hit.incidentMedium = incidentMedium;
		hit.transmissionMedium = transmissionMedium;
		
		calcIntersection(hit, incidentMedium, transmissionMedium, minIntensity);
		
	}
}