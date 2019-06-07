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
package fusionOptics.surfaces;

import net.jafama.FastMath;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import fusionOptics.Util;
import fusionOptics.tracer.Tracer;
import fusionOptics.types.Interface;
import fusionOptics.types.Intersection;
import fusionOptics.types.Medium;
import fusionOptics.types.RaySegment;
import fusionOptics.types.Surface;


/** Spherical Surface - All of the surface of sphere. */
public class Sphere extends Surface {
	
	public int nSectors = 7;
	public int nRingsPerHemisphere = 10;
	public int nPointsInSector = 3;
	
	protected double centre[];	
	protected double radius;	
	
	public Sphere(String name, double centre[], double radius, Interface iface) {
		this(name, centre, radius, null, null, iface);
	}
	/**
	 * @param name
	 * @param centre	The point on the surface, that is an equal distance from every point on the edge (rim). 
	 * @param dishCentreNormal	Unit vector that points towards the curvature centre, from the dish centre.
	 * @param radiusOfCurvature
	 * @param rimRadius
	 * @param frontMedium
	 * @param backMedium
	 * @param interfaceType
	 */
	public Sphere(String name, double centre[], double radius, Medium innerMedium, Medium outerMedium, Interface iface) {
		super(name, innerMedium, outerMedium, iface);
		this.centre = centre;
		this.radius = radius;
	}
	
	@Override
	public double[] getBoundarySphereCentre() { return centre; }

	@Override
	public double getBoundarySphereRadius() { return radius; }


	
	public void rotate(double point[], double matrix[][]){
		for(int i=0;i<3;i++){
			centre[i] -= point[i];
		}
			
		double newCentre[] = new double[3];
		
		for(int i=0;i<3;i++)
			for(int j=0;j<3;j++){
				newCentre[i] += matrix[i][j] * centre[j];
			}
				
		for(int i=0;i<3;i++){
			centre[i] = point[i] + newCentre[i];
		}
	}
	
	@Override
	public void shift(double[] dX) {
		for(int i=0;i<3;i++){
			centre[i] += dX[i];
		}
	}

	@Override
	public boolean findEarlierIntersection(RaySegment ray, Intersection hit) {
	
		double a=1,b=0,c=0;
		
		for(int i=0;i<3;i++){
			b += 2 * (ray.startPos[i] - centre[i]) * ray.dir[i];
			c += ((ray.startPos[i] - centre[i])*(ray.startPos[i] - centre[i]));
		}
		c -= radius * radius;

		if( (b*b) < (4*a*c) ) //no roots = no intersection of sphere		
			return false;
		
		double rtB24AC = FastMath.sqrt(b*b - 4*a*c);
		
		//two roots give two intersection points with sphere
		double s1 = ( -b - rtB24AC ) / (2*a);
		double s2 = ( -b + rtB24AC ) / (2*a);
		
		double X1[] = new double[3];
		double X2[] = new double[3];
		double N1[] = new double[3];
		double N2[] = new double[3];
		
		for(int i=0;i<3;i++){
			X1[i] = ray.startPos[i] + s1 * ray.dir[i];
			X2[i] = ray.startPos[i] + s2 * ray.dir[i];
			N1[i] = centre[i] - X1[i];
			N2[i] = centre[i] - X2[i];
		}
		
		double lenN1=0,lenN2=0; 
		for(int i=0;i<3;i++){
			lenN1 += N1[i]*N1[i];
			lenN2 += N2[i]*N2[i];			
		}
		lenN1 = FastMath.sqrt(lenN1);
		lenN2 = FastMath.sqrt(lenN2);
		
		for(int i=0;i<3;i++){
			N1[i] = N1[i] / lenN1; 
			N2[i] = N2[i] / lenN2; 
		}
			
		boolean lastSurfWasUs = (ray.startHit != null && ray.startHit.surface == this);
		double reHitTolerance =  lastSurfWasUs ? Tracer.reHitTolerance : 0;
	
		//This is copied from Dish, and I'm not sure if it makes sense here
		boolean hit1 = (s1 > reHitTolerance);
		boolean hit2 = (s2 > reHitTolerance);
		
		
		if(hit1){ 
			if(hit2){	//both hit, return the earliest
				if(s1 < s2){
					if(s1 < ray.length){
						hit.pos = X1; hit.normal = N1; ray.length = s1; hit.surface = this;
						return true;
					}
				}else{
					if(s2 < ray.length){						
						hit.pos = X2; hit.normal = N2; ray.length = s2; hit.surface = this;
						return true;
					}
				}
				
			}else{ //only #1 hit
				if(s1 < ray.length){						
					hit.pos = X1; hit.normal = N1; ray.length = s1; hit.surface = this;
					return true;
				}
			}
			
		}else{
			if(hit2){ //only #2 hit
				if(s2 < ray.length){						
					hit.pos = X2; hit.normal = N2; ray.length = s2; hit.surface = this;
					return true;
				}
			}else{ //none hit 
				
			}	
		}	 
		return false; //didn't hit, or hit after it's length (i.e. something somewhere else hit first)
	}

	/** Getters for defining properties */
	public double[] getCentre() { return centre.clone(); }
	public double getRadius() { return radius; }
	
	/** Setters for defining properties */
	public void setCentre(double centre[]){ this.centre = centre;  }
	public void setRadius(double radius) { this.radius = radius;  }
	
	protected double depthFromPlaneR2(double r2){
		return r2 / ( radius * (1.0 + FastMath.sqrt(1.0 - r2/radius/radius)) );		
	}
	
	public List<double[][]> draw() {
		
		ArrayList<double[][]> lines = new ArrayList<double[][]>();
		
		double x[] = new double[]{1, 0, 0};
		double y[] = new double[]{0, 1, 0};
		double z[] = new double[]{0, 0, 1};
		
		double rimRadius = radius*0.9999;
		
		for(int iR=0; iR < nRingsPerHemisphere; iR++){
			double r0, r1, xAS0, xAS1;
			if(iR < (nRingsPerHemisphere/2)){
				r0 = FastMath.sqrt(iR *2.0/ nRingsPerHemisphere) *rimRadius;
				r1 = FastMath.sqrt((iR+1) *2.0/ nRingsPerHemisphere) *rimRadius;				
				xAS0 = radius - depthFromPlaneR2(r0*r0);				
				xAS1 = radius - depthFromPlaneR2(r1*r1);
			}else{
				r0 = FastMath.sqrt((iR-nRingsPerHemisphere/2) *2.0/ nRingsPerHemisphere) *rimRadius;
				r1 = FastMath.sqrt(((iR-nRingsPerHemisphere/2)+1) *2.0/ nRingsPerHemisphere) *rimRadius;			
				xAS0 = -radius + depthFromPlaneR2(r0*r0);				
				xAS1 = -radius + depthFromPlaneR2(r1*r1);
			}
			
			for(int iS=0; iS < nSectors; iS++){
				double phi0 = iS * 2 * Math.PI / nSectors;			
				double phi1 = (iS+1) * 2 * Math.PI / nSectors;
				
				double line[][] = new double[3][2*nPointsInSector+1];
				
				for(int iP = 0; iP < nPointsInSector; iP++){
					double phi = phi0 + (phi1 - phi0) * iP / (nPointsInSector - 1.0);

					for(int k=0; k < 3; k++)
						line[k][iP] = centre[k] + xAS0 * x[k]  
						                       + r0 * FastMath.cos(phi) * y[k]
		            						   + r0 * FastMath.sin(phi) * z[k];
					
					phi = phi1 - (phi1 - phi0) * iP / (nPointsInSector - 1.0);
					for(int k=0; k < 3; k++)
						line[k][nPointsInSector+iP] = centre[k] + xAS1 * x[k]  
						                       + r1 * FastMath.cos(phi) * y[k]
		            						   + r1 * FastMath.sin(phi) * z[k];

					for(int k=0; k < 3; k++)
						line[k][2*nPointsInSector] = centre[k] + xAS0 * x[k]  
								                       + r0 * FastMath.cos(phi0) * y[k]
		 		            						   + r0 * FastMath.sin(phi0) * z[k];
					if(Double.isNaN(line[0][4]))
						throw new RuntimeException("NaN in Sphere");
					
					lines.add(line);
				}
				
			}
		}
		
		
		return lines;
		
	}

	@Override
	public int surfaceGeometryHashCode() {
		long t; int r = 1;
		r = 31 * r + Arrays.hashCode(centre);
		t = Double.doubleToLongBits(radius); 	r = 31 * r + (int) (t ^ (t >>> 32));
		return r;
	}
	
	@Override
	public boolean surfaceGeometryEquals(Surface obj) {
		Sphere other = ((Sphere)obj);
		return Arrays.equals(centre, other.centre)					
					&& Double.doubleToLongBits(radius) == Double.doubleToLongBits(other.radius);
	}
	
}
