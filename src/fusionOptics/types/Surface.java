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

import java.util.LinkedList;
import java.util.List;


/** Base class of any physical surface.
 * Implementations of surface carry the 3D drawing and ray-intersection code.
 * The optical properties of the surface are carried in Interface.
 * 
 *  */
public abstract class Surface implements Element {
	private String name;
	
	protected Medium backMedium; //medium 'behind' the surface (i.e. backward/-ve w.r.t. to the normal)
	
	protected Medium frontMedium; //medium 'in front' of the surface (i.e. forward/+ve w.r.t. to the normal). The normal points into frontMedium
	
	protected List<Medium> media;

	/** optic of which this surfaces forms part */
	private Optic parentOptic;
	
	private  Interface iface;
	
	/** Approximate drawing quality, 0 - 100 where 0 = 1 triangle and 100 is spudily slow */
	protected int approxDrawQuality = 10;
		
	public Surface(String name, Medium frontMedium, Medium backMedium, Interface iface) {
		this.name = name;
		this.frontMedium = frontMedium;
		this.backMedium = backMedium;
		media = new LinkedList<Medium>(); 
		if(frontMedium != null)media.add(frontMedium);
		if(backMedium != null)media.add(backMedium);
		this.iface = iface;
		iface.checkCompatibility(this);
	}
	
	/** Returns the medium 'in front' of the surface (i.e. forward w.r.t. to the normal) */
	public Medium getFrontMedium(){ return frontMedium; }
	
	/** Returns the medium 'behind' the surface (i.e. backward w.r.t. to the normal) */
	public Medium getBackMedium(){ return backMedium; }
	
	/** Sets the medium 'in front' of the surface (i.e. forward w.r.t. to the normal) */
	public void setFrontMedium(Medium frontMedium){ this.frontMedium = frontMedium; }
	
	/** Sets the medium 'behind' the surface (i.e. backward w.r.t. to the normal) */
	public void setBackMedium(Medium backMedium){ this.backMedium = backMedium; }
	
	public final String getName(){ return name; }
	
	/** Returns the parent optic of this surface */
	public Optic getOptic(){ return parentOptic; }
	
	/** Returns lines (lists of double[x/y/z]s), describing the surface */ 
	public abstract List<double[][]> draw();

	public List<Medium> getMedia() { return media; }
	

	public void setInterface(Interface iface) {
		iface.checkCompatibility(this);	
		this.iface = iface;
	}
	
	public Interface getInterface(){ return iface; }

	@Override
	public final int hashCode() {
		int r = surfaceGeometryHashCode();
		r = 31 * r + ((backMedium == null) ? 0 : backMedium.hashCode());
		r = 31 * r + ((frontMedium == null) ? 0 : frontMedium.hashCode());
		r = 31 * r + ((iface == null) ? 0 : iface.hashCode());
		r = 31 * r + ((name == null) ? 0 : name.hashCode());
		return r;
	}
	
	@Override
	public final boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null || getClass() != obj.getClass())
			return false;
		Surface other = (Surface) obj;
		return ((backMedium == null) ? (other.backMedium == null) : backMedium.equals(other.backMedium))
				&& ((frontMedium == null) ? (other.frontMedium == null) : frontMedium.equals(other.frontMedium))
				&& ((iface == null) ? (other.iface == null) : iface.equals(other.iface))
				&& ((name == null) ? (other.name == null) : name.equals(other.name))
				&& surfaceGeometryEquals((Surface)obj);
	}

	/** HashCode of surface geometry implementation (excludes the media, name and interface)*/
	public abstract int surfaceGeometryHashCode();
	
	/** .equals() methods for surface geometry implementation  (excludes the media, name and interface)*/
	public abstract boolean surfaceGeometryEquals(Surface other);
	
	public void setApproxDrawQuality(int approxDrawQuality){
		this.approxDrawQuality = approxDrawQuality;
	}
	
}
