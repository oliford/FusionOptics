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
package fusionOptics.types;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;

import fusionOptics.Util;
import fusionOptics.collection.IntersectionProcessor;
import fusionOptics.pointSpread.PointSpreadFunction;
import fusionOptics.surfaces.Square;
import fusionOptics.tracer.Tracer;

import oneLiners.OneLiners;

import binaryMatrixFile.BinaryMatrixFile;
import binaryMatrixFile.BinaryMatrixWriter;


import jafama.FastMath;

/** Describes the striaght propagation of a ray through a medium*/
public class RaySegment {
	private static final DecimalFormat strFmt = new DecimalFormat("#.##");
	
	/** Starting position of the ray segement */
	public double startPos[];
	
	/** Unit vector direction */
	public double dir[];
	
	/** Vector perpendicular to dir, that defines 'up' for the polarisation described in the stokes vector */
	public double up[];
	
	/** polarisation states at start of ray.
	 * Precisely, these are electric field amplitudes integrated across the width of the ray and along the 
	 * distance the ray would cover in 1 unit of time.
	 * 
	 * A light ray of a given intensity I0 will have the same power/m^2 (energy flux per unit time) once it enters
	 * a  higher n medium, but will have more energy per unit volume and hence higher actual E field, because of the
	 * slower ray velocity. Here, this E doesn't include that, so remains proportional to sqrt(I)
	 */
	public double E0[][];
	
	/** polarisation states at end of ray segment */
	public double E1[][];
	
	/** Length, in real space meters (not optical dist) */
	public double length;
	
	/** The medium through which the ray segment travels */
	public Medium medium; 
	
	/** The wavelength of the ray in vacuum */
	public double wavelength;
	
	/** The refractive index of the ray segment in the medium it's in.
	 * 
	 * This will only exist for media where the ray segements must be treated
	 * with a single polarisation and the refractive index is a bit complicated to
	 * calculate. e.g. the extraordinary ray in birefringent medium
	 * 
	 * Otherwise it will be 0, or NaN
	 */
	public double raySpecificRefractiveIndex;
	
	/** Number of complete wavelengths the ray segment passes (in the medium) */
	public int nWaves;

	/** Intersection point that started the ray segment */
	public Intersection startHit;
	
	/** Intersection point at the end of the ray segment */
	public Intersection endHit;
	
	@Override
	/** Clone just this ray segement - linked segements have the same references as the original */
	protected RaySegment clone() {
		RaySegment newSeg = new RaySegment();
		
		newSeg.length = length;
		newSeg.wavelength = wavelength;
		newSeg.nWaves = nWaves;				
		newSeg.raySpecificRefractiveIndex = raySpecificRefractiveIndex;
		
		newSeg.startPos = startPos.clone();
		newSeg.dir = dir.clone();
		newSeg.up = up.clone();
		newSeg.E0 = Pol.copyAll(E0);
		newSeg.E1 = Pol.copyAll(E1);
		
		newSeg.startHit = startHit;
		newSeg.endHit = endHit;
		newSeg.medium = medium;
		
		return newSeg;
	}
	
	/* ------------------------ Constructors ------------------------*/
	
	/** Creates a blank ray segment to be filled in later */
	public RaySegment() {	
		//deliberatly invalidate everything that must be set
		this.wavelength = Double.NaN;
		this.nWaves = Integer.MIN_VALUE;
		//set the length though, because there really is only 1 sensible default
		this.length = Double.POSITIVE_INFINITY;
	}

	/** Creates a simple common ray segment with a given start position and direction.
	 * 
	 * Sets wavelength to 550nm, a nice greeny-yellow colour.
	 * Initial polarisation is set to linearly polarised parallel to the z direction 
	 * (or to y, if the ray dir is parallel to z)
	 *  
	 * @param start
	 * @param dir
	 */
	public RaySegment(double start[], double dir[]) {
		this.startPos = start;
		this.dir = dir;
		this.wavelength = 550e-9; //550nm, yellow-ey green
		double right[] = Util.cross(dir, new double[]{0,0,1});
		this.up = (Util.length(right) == 0) ? new double[]{0,1,0} :
					Util.reNorm(Util.cross(right, dir));
		this.E0 = new double[][]{{1,0,0,0}}; //lin.pol. in 'up' dir
		this.length = Double.POSITIVE_INFINITY;
		this.nWaves = 0;
		
	}
	
	/** Creates a simple common ray with a given start position, fired randomly at the given surface
	 * 
	 * Sets wavelength to 550nm, a nice greeny-yellow colour.
	 * Initial polarisation is set to linearly polarised parallel to the z direction 
	 * (or to y, if the ray dir is parallel to z)
	 *  
	 * @param start
	 * @param dir
	 */
	public RaySegment(double start[], Element target) {
		this.startPos = start;
		this.dir =  Tracer.generateRandomRayTowardSurface(start, target);
		this.wavelength = 550e-9; //550nm, yellow-ey green
		double right[] = Util.cross(dir, new double[]{0,0,1});
		this.up = (Util.length(right) == 0) ? new double[]{0,1,0} :
					Util.reNorm(Util.cross(right, dir));
		this.E0 = new double[][]{{1,0,0,0}}; //lin.pol. in 'up' dir
		this.length = Double.POSITIVE_INFINITY;
		this.nWaves = 0;
		
	}
	
	public RaySegment(double start[], double dir[], double wavelength) {
		this(start,dir);
		this.wavelength = wavelength;
	}
	
	
	/** Creates an infinite ray with a given start position, direction and 'up' definition.
	 * Sets the polarisation list to that used by the PSF system for Mueller characterisation.
	 * 	 *  
	 * @param start
	 * @param dir
	 * @param up	The sense of 'up' for the polarisation states
	 * @param wavelength
	 */
	public RaySegment(double start[], double dir[], double up[], double wavelength) {
		this.startPos = start;
		this.dir = dir;
		this.up = up;
		this.wavelength = wavelength;
		this.E0 = PointSpreadFunction.getInputStatesForMuellerCalc();
		this.length = Double.POSITIVE_INFINITY;
		this.nWaves = 0;
	}
	

	/* ------------------------ Helpers, Setters and Getters ------------------------*/
	
	/** Changes the the ray's sense of polarisation 'up', and
	 * appropriately adjusts the stokes vector so that the actual 
	 * polarisation doesn't change.
	 *  
	 * @param upIsh Any vector not parallel to the ray direction.
	 */
	public final void rotatePolRefFrame(double upIsh[]) {
		double newRight[] = Util.reNorm(Util.cross(dir, upIsh));
		double newUp[] = Util.cross(newRight, dir);
		
		double oUnU = Util.dot(up, newUp);
		double oUnR = Util.dot(up, newRight);
		
		if(E0 != null){
			for(int i=0; i < E0.length; i++)
				E0[i] = Pol.rotateFrame(E0[i], oUnU, oUnR);
		}
		if(E1 != null){
			for(int i=0; i < E1.length; i++)
				E1[i] = Pol.rotateFrame(E1[i], oUnU, oUnR);
		}
		
		this.up = newUp;
		
	}

	private class CollectHits implements IntersectionProcessor {
		List<Intersection> hitList;
		@Override
		public void nextIntersection(Intersection hit) {
			hitList.add(hit);
		}
	}
	
	
	/** Returns a list of intersections of the whole ray tree, starting at this one, with the given object */ 
	public final List<Intersection> getIntersections(Element element) {
		CollectHits hitCollector = new CollectHits();
		hitCollector.hitList = new LinkedList<Intersection>();
		if(element instanceof Surface){
			processIntersections((Surface)element, hitCollector, this);
		}else if(element instanceof Optic){
			for(Surface s : ((Optic)element).getSurfacesAll())
				processIntersections(s, hitCollector, this);				
		}
		return hitCollector.hitList;
	}
		
	/** Returns a list of intersections of the whole ray tree, starting at this one, with the given object */ 
	public final List<Intersection> getIntersections(Surface element, List<Intersection> hitList) {
		CollectHits hitCollector = new CollectHits();
		hitCollector.hitList = hitList;
		processIntersections(element, hitCollector, this);
		return hitCollector.hitList;
	}
	
	/** Calls the given IntersectionProcessors for every contact this ray makes with the given target surface */
	public final int processIntersections(Surface target, IntersectionProcessor... hitProcs) {
		return processIntersections(target, hitProcs, this);		
	}
	
	/** Calls the given IntersectionProcessors for every contact this ray makes with the given target surface */
	private static final int processIntersections(Surface target, IntersectionProcessor hitProc, RaySegment ray) {
		return processIntersections(target, new IntersectionProcessor[]{ hitProc }, ray);
	}
	
	/** Fills the given list with intesections from the whole ray tree, starting at this one, with the given object */ 
	private static final int processIntersections(Surface target, IntersectionProcessor[] hitProcs, RaySegment ray) {
		int nHits = 0;
		
		while(ray.endHit != null) {
			if(target == null || ray.endHit.surface == target){
				for(IntersectionProcessor hitProc : hitProcs)
					if(hitProc != null)
						hitProc.nextIntersection(ray.endHit);
				nHits++;
			}
			
			//keep going while there are outgoing rays...
			if(ray.endHit.isEnd())
				break;
		
			//get them, with the strongest first
			//List<RaySegment> nextSegs = ray.endHit.getSortedOutgoingRays();
			RaySegment[] nextSegs = ray.endHit.getSortedOutgoingRaysArray();

			//for each new branch... (every outgoing ray excluding the strongest) 
			for(int i=1; i < nextSegs.length; i++){
				if(nextSegs[i] != null)
					for(IntersectionProcessor hitProc : hitProcs)
						nHits += processIntersections(target, hitProc, nextSegs[i]);	
			}
			
			//continue with strongest
			ray = null;
			for(int i=0; i < nextSegs.length; i++){
				if(nextSegs[i] != null){
					ray = nextSegs[i];
					break;
				}
			}
		}
		
		return nHits;
	}
	
	/** Rough sense of total intensity for the ray segment
	 * @return Total intensity of all the polarisations. 
	 */
	public double startIntensity() {
		double I = 0;
		for(int i=0; i < E0.length; i++)
			I += Pol.intensity(E0[i]);
		return I;
	}
	
	/** Rough sense of total intensity for the ray segment
	 * @return Total intensity of all the polarisations. 
	 */
	public double endIntensity() {
		return Pol.intensity(E1);
	}
	
	/** Dumps the ray path/tree to stdout */
	public void dumpPath() { dumpPath(this, 0);	}
		
	/** Dumps the ray path/tree to stdout */
	public static void dumpPath(RaySegment ray, int depth){
		while(ray.endHit != null) {
			
			for(int i=0; i<depth;i++)
				System.out.print(" ");
			System.out.println(ray.toString());
			
			//keep going while there are outgoing rays...
			if(ray.endHit.isEnd())
				break;
		
			//get them, with the strongest first
			List<RaySegment> nextSegs = ray.endHit.getSortedOutgoingRays();

			//for each new branch... (every outgoing ray excluding the strongest) 
			for(int i=1; i < nextSegs.size(); i++){
				dumpPath(nextSegs.get(i), depth+1);	
			}
			
			//continue with strongest
			ray = nextSegs.get(0);			
		}
	}
	
	
	@Override
	public String toString() {
		if(dir[0] == 0 && dir[1] == 0 & dir[2] == 1)
			rotatePolRefFrame(new double[]{ 0, 1, 0 });
		else
			rotatePolRefFrame(new double[]{ 0, 0, 1 });

		StringBuilder str = new StringBuilder();
		str.append("Ray [");
		str.append((startHit == null) ? "null" : startHit.toString());
		if(E0 != null){
			str.append(",I=");
			str.append(strFmt.format(startIntensity()));
			for(int i=0; i < E0.length; i++)
				str.append("{"+Pol.toString(E0[i])+"}");
			
			str.append("] --> [");
		}else{
			str.append(",noPol");
		}
		
		str.append((endHit == null) ? "null" : endHit.toString());
		if(E1 != null){
			str.append(",I=");
			str.append(strFmt.format(endIntensity()));		
			for(int i=0; i < E1.length; i++)
				str.append("{"+Pol.toString(E1[i])+"}");
		}else{
			str.append(",noPol");
		}
		str.append("]");
		return str.toString();
	}

	/** Trace backwards through the ray path, from this segment until an intersection with the given surface is found */ 
	public Intersection findFirstEarlierIntersection(Surface target) {
		RaySegment ray = this;
		
		while(ray != null && ray.startHit != null){
			if(ray.startHit.surface == target)
				return ray.startHit;
			
			ray = ray.startHit.incidentRay;
		};
		return null;
	}

	
	/** Dumps the ray path/tree to stdout */
	public LinkedList<Surface> getPrimarySurfaceSequence() { 
		LinkedList<Surface> elemList = new LinkedList<Surface>();
		
		RaySegment ray = this;
		
		while(ray.endHit != null) {
			if(ray.endHit.surface != null)
				elemList.add(ray.endHit.surface);
			
			//get them, with the strongest first
			List<RaySegment> nextSegs = ray.endHit.getSortedOutgoingRays();
			if(nextSegs.size() <= 0)
				break;
			
			ray = nextSegs.get(0);
		}
		
		if(ray.endHit != null && ray.endHit.surface != null)
			elemList.add(ray.endHit.surface);
		
		return elemList;
	}

	/** Creates a ray which combines transmittedOrdinary and transmittedExtraordinary 
	 * from parallel ray paths. The result doesn't make any particular sense, however
	 * it's nice to look at. */
	public RaySegment createMergedRay() {
		LinkedList<RaySegment> raySegs = new LinkedList<RaySegment>();
		
		raySegs.add(this);
		RaySegment start = null;
		RaySegment current = null;
		
		do{
			double sumI = 0;
			double avgEndPos[] = new double[3];
			double avgNormal[] = new double[3];
			double sumE0[][] = null;
			double sumE1[][] = null;
			
			RaySegment ray0 = null;
			ListIterator<RaySegment> iter = raySegs.listIterator();
			while(iter.hasNext()){
				RaySegment ray = iter.next();
				iter.remove();
				
				if(ray.endHit == null) //ignore it
					continue;
				
				if(ray0 == null){ //compare everything with the first valid ray
					ray0 = ray;
				}else{
					//checks
					if(ray.endHit.surface != ray0.endHit.surface){ //they should at least all hit the same surface
						System.err.println("RaySegment.createMergedRay(): Ray paths diverge too much. Two 'parallel' segments hit different surfaces.");
						continue;
					}
					
					if(FastMath.acos(Util.dot(ray.dir, ray0.dir)) > 40*Math.PI/180){
						System.err.println("RaySegment.createMergedRay(): Ray paths diverge too much. Two 'parallel' segments diverge by more than 40°.");
						continue;
					}
					
					//and get the polarisations the same way up
					ray.rotatePolRefFrame(ray0.up);
				}
					
				//average the end position (intensity weighted)
				double I = Pol.intensity(ray.E1);
				avgEndPos[0] += I * ray.endHit.pos[0];
				avgEndPos[1] += I * ray.endHit.pos[1];
				avgEndPos[2] += I * ray.endHit.pos[2];
				avgNormal[0] += I * ray.endHit.normal[0];
				avgNormal[1] += I * ray.endHit.normal[1];
				avgNormal[2] += I * ray.endHit.normal[2];
				sumI += I;

				//add the polarisations
				if(sumE0 == null){
					sumE0 = Pol.copyAll(ray.E0);
				}else{
					for(int i=0; i < sumE0.length; i++)
						for(int j=0; j < sumE0[0].length; j++)
							sumE0[i][j] += ray.E0[i][j];
				}
					
				//add end polarisations for rays that have them
				if(ray.E1 != null){
					if(sumE1 == null){
						sumE1 = Pol.copyAll(ray.E1);
					}else{
						for(int i=0; i < sumE1.length; i++)
							for(int j=0; j < sumE1[0].length; j++)
								sumE1[i][j] += ray.E1[i][j];
					}
				}

				//replace the ray segment with all outgoing segs for the next round
				if(ray.endHit.transmittedOrdinary != null)
					iter.add(ray.endHit.transmittedOrdinary);
				if(ray.endHit.transmittedExtraordinary != null)
					iter.add(ray.endHit.transmittedExtraordinary);				
			}
			
			if(ray0 == null) //never found a complete valid segment
				break;
			
			avgEndPos[0] /= sumI;
			avgEndPos[1] /= sumI;
			avgEndPos[2] /= sumI;
			avgNormal[0] /= sumI;
			avgNormal[1] /= sumI;
			avgNormal[2] /= sumI;

			Intersection endHit = new Intersection();
			endHit.pos = avgEndPos;
			endHit.normal = avgNormal;			
			endHit.surface = ray0.endHit.surface;
			
			RaySegment newSeg = new RaySegment();
			newSeg.startPos = (current != null) ? current.endHit.pos.clone() : startPos.clone();
			newSeg.endHit = endHit;
			newSeg.E0 = sumE0;
			newSeg.E1 = sumE1;
			double rayVec[] = Util.minus(avgEndPos, newSeg.startPos);
			newSeg.length = Util.length(rayVec);
			newSeg.dir = Util.reNorm(rayVec);
			newSeg.medium = ray0.medium;
			newSeg.nWaves = 0; //this really doesn't make sense
			newSeg.wavelength = ray0.wavelength;
			newSeg.up = ray0.up.clone();
			endHit.incidentRay = newSeg;
			
			if(current == null){
				start = newSeg;				
			}else{
				current.endHit.transmittedOrdinary = newSeg;
				current.endHit.transmittedExtraordinary = null;
				current.endHit.reflectedOrdinary = null;
				current.endHit.reflectedExtraordinary = null;
			}
			current = newSeg;
			
		}while(true);
		
		
		return start;
	}
	
	
}
