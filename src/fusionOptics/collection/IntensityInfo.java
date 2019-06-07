package fusionOptics.collection;

import net.jafama.FastMath;

import java.security.KeyStore.Entry;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

import otherSupport.ScientificNumberFormat;

import fusionOptics.Util;
import fusionOptics.surfaces.Square;
import fusionOptics.types.Element;
import fusionOptics.types.Intersection;
import fusionOptics.types.Optic;
import fusionOptics.types.RaySegment;
import fusionOptics.types.Surface;


/** Records the total number and intensity of rays reaching
 * each component in a system 
 * 
 * Also records average total ray length
 * 
 * @author oliford
 *
 */
public class IntensityInfo implements IntersectionProcessor {
	ScientificNumberFormat fmt = new ScientificNumberFormat("#.###", "#.###E0", 3);
	
	private class Info { 
		String name; 
		int nRays = 0; 
		double intensity = 0;
		double length = 0;
	}
	
	private HashMap<Surface, Info> collect = new HashMap<Surface, IntensityInfo.Info>();
	
	public IntensityInfo(){
		
	}
	
	public IntensityInfo(Element elem) {
		fillNames(elem, "");
	}
	
	private void fillNames(Element elem, String nameStr) {
		nameStr += elem.getName();
		
		if(elem instanceof Optic){
			for(Optic optic : ((Optic)elem).getSubOptics())
				fillNames(optic, nameStr + ".");
			for(Surface surface : ((Optic)elem).getSurfaces())
				fillNames(surface, nameStr + ".");
		}else{
			Info info = new Info();
			info.name = nameStr;
			collect.put((Surface)elem, info);
		}
	}
	
	@Override
	public void nextIntersection(Intersection hit) {
		
		Info info = collect.get(hit.surface);
		if(info == null){
			info = new Info();
			collect.put(hit.surface, info);
			info.name = "[unknownSurface]." + hit.surface.getName();
		}
		
		info.nRays++;
		info.intensity += hit.incidentRay.endIntensity();
		
		double len = 0;
		RaySegment ray = hit.incidentRay;
		while(ray != null){
			len += ray.length;
			ray = ray.startHit == null ? null : ray.startHit.incidentRay;
		}
		info.length += len;
	}
	
	public void dump(){ dump(null, null, true, 0, 0); }
	
	public void dump(boolean sort){ dump(null, null, sort, 0, 0); };
	
	/**
	 * @param interested 	Collection of surfaces of interest (what to dump)
	 * @param relativeTo	If not null, also display info relative to rays hitting this surface
	 * @param sort			If true, sort the surfaces.
	 * @param firingTargetSolidAngle	Solid angle of the firing target used to generate rays. 
	 * 								If non-zero, the effective solid angle from the source of rays reaching the given surface is also displayed.
	 * 								 
	 * @param nAttemptedRays	The number of attempts at firing rays
	 */
	public void dump(Collection<Surface> interested, Surface relativeTo, boolean sort, double firingTargetSolidAngle, int nAttemptedRays){
				
		if(interested == null){
			interested = collect.keySet();
		}
		
		if(sort){
			List<Surface> sortedList = new ArrayList<Surface>(interested);
			Collections.sort(sortedList, new Comparator<Surface>() {
				@Override
				public int compare(Surface surf1, Surface surf2) {
					Info info1 = collect.get(surf1);
					Info info2 = collect.get(surf2);
					return Double.compare(
							(info1 != null) ? info1.intensity : 0,
							(info2 != null) ? info2.intensity : 0);
				}});
			interested = sortedList;
		}
		
		Info relativeInfo = null;
		if(relativeTo != null){
			relativeInfo = collect.get(relativeTo);
		}
		
		for(Surface surface : interested){
			Info info = collect.get(surface);
			
			if(info == null){
				System.out.println(surface.getName() + ": No Data");
				continue;
			}
			
			System.out.print(info.name + ":\tn=" + info.nRays + "\tI=" + fmt.format(info.intensity));
			if(firingTargetSolidAngle > 0){
				//effective solid angle = firing target solid angle * ray fraction making it to this surface
				System.out.print("\tSrcSldAng = " + fmt.format(firingTargetSolidAngle * info.nRays / nAttemptedRays) + "  sr" );
			}
			if(relativeInfo != null){
				System.out.print("\t" 
						+ "[ Rel: n=" + fmt.format(info.nRays * 100.0 / relativeInfo.nRays) + " %\t" 
						+ "I=" + fmt.format(info.intensity * 100.0 / relativeInfo.intensity) +" % ]");
			}			
			System.out.println("\tposOnRay=" + fmt.format(info.length / info.nRays));
			
		}
	}
	
	/** Clears the counts and intensities, not the surface list/names */
	public void reset(){ 
		for(java.util.Map.Entry<Surface, Info> entry : collect.entrySet()){
			Info info = entry.getValue();
			info.nRays = 0;
			info.intensity = 0;
			info.length = 0;
		}
	}

	/** Returns the effective solid angle from the source that is collected and ends up hitting the given surface */ 
	public double getSourceSolidAng(Surface surface, double firingTargetSolidAngle, int nAttemptedRays) {
		Info info = collect.get(surface);
		if(info == null)
			throw new RuntimeException("Surface " + surface.getName() + " not in list");
		
		return firingTargetSolidAngle * info.nRays / nAttemptedRays;
	}
}
