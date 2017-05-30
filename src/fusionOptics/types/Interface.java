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

import java.util.List;

import fusionOptics.Util;


/** Optical properties of a surface.
 * 
 *  Implementations of this class should handle all
 *  the various interface types and approximations.
 *  
 *  It does not need to be instantiated for each surface.
 *  
 */
public interface Interface {

	public void calcIntersection(Intersection newHit, double minIntensity);
	
	public void checkCompatibility(Surface surface);
	
}
