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

/** Any 3D space element which can be intersected, 
 *  rotated, shifted and fit inside bounding sphere */
public interface Element {
	
	/** Returns the bounding sphere centre of this surface */  
	public abstract double[] getBoundarySphereCentre();
	
	/** Returns the bounding sphere radius of this surface */
	public abstract double getBoundarySphereRadius();
	
	/** Tests if the given ray intersects this surface and if so, fills some information 
	 * in to the given intersection object:  surface, pos[], normal[]
	 * This includes testing if it is within the current ray segment length, 
	 *   so calling this for multiple surfaces will leave it with the first intersection. */
	public abstract boolean findEarlierIntersection(RaySegment ray, Intersection hit);
		
	/** Applies a translation to this optic */
	public abstract void shift(double dX[]);
	
	/** Applies a rotation to this optic */
	public abstract void rotate(double point[], double matrix[][]);

	public abstract String getName();

	public abstract void setApproxDrawQuality(int approxDrawQuality);
		
	/** Some kind of 'centre' definition for the element */
	public abstract double[] getCentre();
}
