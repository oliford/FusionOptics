/**
 * Copyright 2011 Oliver Ford
 *
 * This file is part of the minerva-optics 'RayTracer'.
 *
 *   RayTracer is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   RayTracer is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with RayTracer.  If not, see <http://www.gnu.org/licenses/>.
 *   
 *   @author oliford <codes<at>oliford.co.uk>
 */
package fusionOptics.tracer;


import net.jafama.FastMath;

import java.util.List;

import fusionOptics.Util;
import fusionOptics.interfaces.Absorber;
import fusionOptics.interfaces.NullInterface;
import fusionOptics.types.Element;
import fusionOptics.types.Interface;
import fusionOptics.types.Intersection;
import fusionOptics.types.Medium;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;
import fusionOptics.types.Surface;
import otherSupport.RandomManager;


/**
 *  Top level tracer stuff. All static methods.
 *  
 *  The main process of the ray tracer is in Tracer.trace(). 
 * 
 * The guts of the calculations are in other places:
 * 		Surface - Ray direction and intersection calc.
 * 		Medium - Ray phase change and wave counting as it propagates.
 * 		Interface - Ray splitting, transmission/reflection and surface polarisation effects.   
 * 
 */
public abstract class Tracer {
	/** Distance after leaving an object before which we ignore new intersections with that object.
	 * Set globally to something significantly smaller than your problem size 
	 * {@link #generateRandomRayTowardSurface(double[], Element, boolean)} */ 
	public static double reHitTolerance = 1e-10;
	
	/** {@link generateRandomRayTowardSurface(RandomManager randomManager, double start[], Element target)} with default RandomManager */
	public static double[] generateRandomRayTowardSurface(double start[], Element target, boolean lengthAsSolidAngle){
		return generateRandomRayTowardSurface(RandomManager.instance(), start, target, lengthAsSolidAngle);
	}
	
	/** {@link generateRandomRayTowardSurface(RandomManager randomManager, double start[], Element target)} with default RandomManager and returns unit length vector*/
	public static double[] generateRandomRayTowardSurface(double start[], Element target){
		return generateRandomRayTowardSurface(RandomManager.instance(), start, target, false);
	}
	
	/** {@link generateRandomRayTowardSurface(RandomManager randomManager, double start[], Element target)} with unit length vector*/
	public static double[] generateRandomRayTowardSurface(RandomManager randomManager, double start[], Element target){
		return generateRandomRayTowardSurface(randomManager, start, target, false);
	}
		
	/** Returns vector randomly toward given element from given point.
	 * Vectors are uniformly distributed over solid angle. see:
	 * http://mathworld.wolfram.com/SpherePointPicking.html 
	 * 
	 * The length of the returned vector is equal to the solid angle (sterradians) subtended by the total target firing
	 * range (i.e. related to the fraction of rays that can be chosen compared to rays that would be fired in all directions)
	 * */
	public static double[] generateRandomRayTowardSurface(RandomManager randomManager, double start[], Element target, boolean lengthAsSolidAngle){
		
		double centre[] = target.getBoundarySphereCentre();
		double R = target.getBoundarySphereRadius();

		//vector from start point to centre of bounding sphere
		double vecToSphereCentre[] = new double[]{ 
				centre[0] - start[0], 
				centre[1] - start[1], 
				centre[2] - start[2],			
			};
		
		double distToSphereCentre = 
			FastMath.sqrt(vecToSphereCentre[0]*vecToSphereCentre[0] +
					vecToSphereCentre[1]*vecToSphereCentre[1] +
					vecToSphereCentre[2]*vecToSphereCentre[2]);
	
		vecToSphereCentre[0] /= distToSphereCentre;
		vecToSphereCentre[1] /= distToSphereCentre;
		vecToSphereCentre[2] /= distToSphereCentre;
		
		// calc max theta (away from axis to bounding sphere centre)
		//double maxTheta = FastMath.atan(R / distToSphereCentre);
		//double cosMaxTheta = FastMath.cos(maxTheta);
		
		//trig avoidance:
		double rtD2R2 = FastMath.sqrt(distToSphereCentre*distToSphereCentre + R*R);
		double cosMaxTheta = distToSphereCentre / rtD2R2;
		
		//we would generate uniform on cos(theta) in (-1, 1), but we only to go up to maxTheta
		double cosTheta = 1 - RandomManager.instance().nextUniform(0, 1) * (1 - cosMaxTheta);		
		double sinTheta = FastMath.sqrt(1 - cosTheta*cosTheta);
		
		double phi = randomManager.nextUniform(0, 1) * 2 * Math.PI;
		
		//generate in coord sys (a,b,c) with c as axis toward target 
		double a = sinTheta * FastMath.cos(phi);
		double b = sinTheta * FastMath.sin(phi);
		double c = cosTheta;
		
		//we now need a and b unit vector in any direction, but perp to each other and c
		double aU[] = Util.createPerp(vecToSphereCentre);
		double bU[] = Util.cross(vecToSphereCentre, aU);		
		Util.reNorm(aU);
		Util.reNorm(bU);
				
		//convert 
		double unit[] = new double[]{
			c * vecToSphereCentre[0] + a * aU[0] + b * bU[0],
			c * vecToSphereCentre[1] + a * aU[1] + b * bU[1],
			c * vecToSphereCentre[2] + a * aU[2] + b * bU[2]
		};
		
		Util.reNorm(unit);
		
		if(lengthAsSolidAngle){
			double targetSolidAngle = 2*Math.PI * (1.0 - cosMaxTheta); //exact, see https://en.wikipedia.org/wiki/Solid_angle#Cone,_spherical_cap,_hemisphere
			unit[0] *= targetSolidAngle;
			unit[1] *= targetSolidAngle;
			unit[2] *= targetSolidAngle;
		}
				
		return unit;
	}

	
	/** Build the RaySegments tree from rayStart by tracing it through the given Optic
	 * If multiPath is set, all ray paths are traced, recursively where necessary.
	 * Otherwise only the strongest intensity at each split is traced.
	 * Created rays under the given intensity are ignored.
	 * If the total generated ray segement count exceeds maxSegments, the tracing is aborted.
	 * The strongest ray path is traced first.
	 * 
	 * @param elem			Element or collection of elements to trace the ray through, in any order.
	 * @param rayStart		A initialised RaySegment describing the start of the ray.
	 * @param maxSegments	Maximum number of ray segment before aborting (strongest complete path will be returned first)
	 * @param minIntensity	Minimum intensity, below which ray segments will not be traced further.
	 * @param multiPath		If true, all split rays (part reflection, birefringent refraction etc) are followed recursively, strongest first.
	 * @return 			Total number of segments in tree.
	 */
	public static int trace(Element elem, RaySegment rayStart, int maxSegments, double minIntensity, boolean multiPath){
		return trace(elem, null, rayStart, maxSegments, minIntensity, multiPath);
	}
	
	/** Trace a single ray through the given predefined sequence of surfaces, building the ray tree.
	 * 
	 * This can be used to optimise tracing a single ray through a known optic sequence.
	 * (It won't work with birefingent objects, which only makes sense with multipath tracing)
	 *
	 * If initially unknown, the optic sequence can be obtained from ray.getPrimarySurfaceSequence() 
	 * after using the full trace(Element, ... ) to find a good ray (that reaches the final target).
	 * 
	 * However, be careful of e.g. lenses with irises built in. If the first ray goes through the 
	 * centre of all the lenses, never hitting an iris, further rays will not be tested against that iris
	 * and might pass through it.
	 * 
	 * @param seqeunce		List of surfaces to trace, in order.
	 * @param rayStart		A initialised RaySegment describing the start of the ray.
	 * @param maxSegments	Maximum number of ray segment before aborting (strongest complete path will be returned first)
	 * @param minIntensity	Minimum intensity, below which ray segments will not be traced further. 
	 * @return 			Total number of segments in tree.
	 */
	public static int trace(List<Surface> seqeunce, RaySegment rayStart, int maxSegments, double minIntensity){
		return trace(null, seqeunce, rayStart, maxSegments, minIntensity, false);
	}
	
	/** The guts of the tracer. See public trace() methods */
	private static int trace(Element elem, List<Surface> surfaceSequence, RaySegment rayStart, int maxSegments, double minIntensity, boolean multiPath){
		
		RaySegment ray = rayStart;
		int nSegments = 0;
				
		do{
			Intersection newHit = new Intersection();
			newHit.incidentRay = ray;
			
			//find if/what it hits. If we're tracing a given sequence, it MUST hit that surface, otherwise is can be anything.
			boolean anyHits;
			if(surfaceSequence == null)
				anyHits = elem.findEarlierIntersection(ray, newHit);
			else
				anyHits = (nSegments >= surfaceSequence.size()) ? false : surfaceSequence.get(nSegments).findEarlierIntersection(ray, newHit);
			
			if(!anyHits) //if nothing, give up
				return nSegments;
						
			ray.endHit = newHit;

			if(ray.length < 2*reHitTolerance && ray.startHit != null && ray.startHit.surface != null && ray.endHit.surface != null){
				Interface s1 = ray.startHit.surface.getInterface();
				Interface s2 = ray.endHit.surface.getInterface();
				if(!(s1 instanceof Absorber) && !(s2 instanceof Absorber) 
						//&& !(s2 instanceof NullInterface) && !(s1 instanceof NullInterface)  
						// !!- Null interfaces screw up as well, because the intersection finder can find an intersection
						// backwards along the ray
						){
					proximityWarning(ray);
				}
			}
			
			// apply changes to intensity/polarisation from the medium (if there is one)
			if(ray.medium != null)
				ray.medium.rayPropagation(ray);
			else
				Medium.rayPropagation(ray, 1.0, 1.0, 1.0, 1.0); //vacuum propagation
				
			//propagation in medium might have killed the ray's intesnity to the point we don't care
			if(ray.endIntensity() < minIntensity)
				return nSegments;
			
			// calculate what happens to the ray at the end hit (also creates new rays)
			newHit.surface.getInterface().calcIntersection(newHit, minIntensity);
			
			List<RaySegment> newRays = newHit.getSortedOutgoingRays();
			
			if(newRays.size() <= 0)
				return nSegments; //ray was absorbed at hit

			nSegments++; //count the segment we just did
			
			if(multiPath){ // (non-sequence only)
				//if multipath is active, then for all rays except the last one, 
				// we need to recurse this routine
				//do the strongest first, so we always do the stronger paths if we run out of maxSegments
				while(newRays.size() > 1){
					RaySegment newRay = newRays.remove(0);
					nSegments += trace(elem, newRay, maxSegments - nSegments, minIntensity, true);
				}
			}
			
			//we trace the last(weakest) one in here (or the strongest one, if multipath is disabled)
			ray = newRays.remove(0);
			newRays.clear();
			
		}while(nSegments < maxSegments);

		return nSegments;
	}

	/** Returns the ray perpendicular direction for the given ray, that projects into the 
	 * direction 'up' in the ray fan perpendicular plane. That plane being perpendicular to
	 * the line startPos -> target.
	 *
	 * Generating a ray fan using this will give each ray a slightly different sense of 'up'
	 * but when projected into a plane with Pol.projectToPlanesView(..., false), or collimated
	 * using a perfect lens, all the polarisations will be exactly parallel/equal.
	 * 
	 * I don't know yet if this matches any particular physical situation, but it does allow
	 * you to isolate the depolarisating effect of an imaging system, or components of it
	 * without the result being dominated by the problem of defining a consistent single initial
	 * polarisation for diverging rays. 
	 * 
	 * @param startPos	The start of all the diverging rays.
	 * @param rayDir	The direction of this ray.
	 * @param target	The centre target position at which the ray fan centre is directed.
	 * @param centralUp	The direction in the start-taget line perp plane in which the projected polarsation should be.
	 * @return	The ray perp 'up' which projects correctly.
	 */
	public static double[] generateRayFanConsistentPolarisationDefinition(
			double[] startPos, double[] rayDir, double[] target,
			double[] fanUp) {
		
		double n[] = Util.reNorm(Util.minus(startPos, target)); //the plane normal
		fanUp = Util.reNorm(Util.cross(Util.cross(n, fanUp), n)); //make sure its perp
		double s[] = Util.reNorm(Util.cross(rayDir, n));
		double p[] = Util.cross(s, rayDir);
		double As = Util.dot(fanUp, s);
		double Ap = Util.dot(fanUp, Util.cross(n, s));
		
		//System.out.println(Util.dot(s, s));
		//System.out.println(Util.dot(p, p));
		//System.out.println(Util.dot(Util.cross(n,s), Util.cross(n,s)));
		
		double rayUp[] = new double[]{
				As * s[0] + Ap * p[0],
				As * s[1] + Ap * p[1],
				As * s[2] + Ap * p[2],
			};
		
		//System.out.println(Util.dot(v, s) + "\t" + Util.dot(ray.up, s) + "\t" + As);
		//System.out.println(Util.dot(v, Util.cross(n, s)) + "\t" + Util.dot(ray.up, p) + "\t" + Ap);
		
		//System.out.println(Util.dot(ray.up, ray.up) + "\n");
		return rayUp;
	}

	/** Generates random parallel rays toward an elements bounding sphere */
	public static double[] generateRandomInfinityRayStartPos(double dir[], Element target, double rayLength) {
		
		double aU[] = Util.createPerp(dir);
		double bU[] = Util.cross(dir, aU);		
		Util.reNorm(aU);
		Util.reNorm(bU);
		
		//generate a random point on the circle made by the projection of the target bounding sphere
		// on to the perpendicular plane to rayDir0
		
		double r = target.getBoundarySphereRadius() * FastMath.sqrt(RandomManager.instance().nextUniform(0, 1));
		double ang = RandomManager.instance().nextUniform(0, 1) * 2*Math.PI;
		double lA = r * FastMath.sin(ang);
		double lB = r * FastMath.cos(ang);
		double cc[] = target.getBoundarySphereCentre();
		
		return new double[]{
				cc[0] + lA * aU[0] + lB * bU[0] - rayLength * dir[0],
				cc[1] + lA * aU[1] + lB * bU[1] - rayLength * dir[1],
				cc[2] + lA * aU[2] + lB * bU[2] - rayLength * dir[2],
		};
	}

	public static int proximityWarningCount = 0;
	public static final int maxProximityWarningCount = 100;
	
	/** This is warning for proximity of two surfaces that has casued a ray to hit both at almost exactly the same time.
	 * Surfaces are supposed to ignore hits if they themselves are hit again within the tolerance since,
	 * if the same surface is hit twice, we can probably safely ignore it. This is usually cases like curved mirrors or
	 * lenses where the numerics brings us back into, or back inside the surface. The first contact is usually a good
	 * approximation to what should really happen.
	 * 
	 * For two close (probably intersecting) surfaces, this is more of a problem, even if one or both are NullInterface surfaces.
	 * The reason is that the ray may continually hit the /other/ surface and never leave.
	 * 
	 * The alternative would be to entirely ignore any within tolerance hit, but then it is then possible to lose rays to the vacuum
	 * that think they are still inside a material, by skipping the second surface.
	 * 
	 * The proper solution is to separate any refracting/reflecting surfaces, and put an absorbing iris, cylinder or box to catch the edge rays.
	 * 
	 * We can safely ignore this if the new surface is a full Absorber, since the ray will soon stop regardless.
	 *  
	 */ 
	public static final void proximityWarning(RaySegment ray) {
		if(proximityWarningCount < maxProximityWarningCount){
			String surf1 = (ray.startHit == null || ray.startHit.surface == null) ? "NULL" : ray.startHit.surface.getName(); 
			String surf2 = (ray.endHit == null || ray.endHit.surface == null) ? "NULL" : ray.endHit.surface.getName(); 
					
			System.err.println("Tracer Surface Proximity Warning: Ray hit surface " + surf2 + " a distance " + ray.length + " after hitting " + surf1);
			
			if(proximityWarningCount >= maxProximityWarningCount){
				System.err.println("*** TOO MANY PROXIMITY WARNINGS, SURPRESSING FUTHER ERRORS *** ");
			}
		}		
		proximityWarningCount++;
	}
	
}