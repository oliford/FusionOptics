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


/** Spherical Surface - Part of the surface of sphere.
 * (Also the base and first-order approximation class for Aspheric) */
public class Dish extends Surface {
	/*public int nSectors = 8;
	public int nRings = 6;
	public int nPointsInSector = 3;*/

	
	protected double centre[];
	protected double dishDiameter;
	protected double radiusOfCurv;
	protected double dishNormal[];
	protected double thetaCrit;
	protected double curvCentre[];
	protected double boundSphereRadius;
	
	/** @return The depth of a surface with given curvatureRadius at the given radius around it's symmetry axis */
	public static double depth(double radiusOfCurvature, double radius) {
		return radiusOfCurvature - Math.sqrt(radiusOfCurvature*radiusOfCurvature - radius*radius);
	}
	
	public Dish(String name, double centre[], double dishCentreNormal[], double radiusOfCurvature, double rimRadius, Interface iface) {
		this(name, centre, dishCentreNormal, radiusOfCurvature, rimRadius, null, null, iface);
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
	public Dish(String name, double centre[], double dishCentreNormal[], double radiusOfCurvature, double rimRadius,
			Medium frontMedium, Medium backMedium, Interface iface) {
		super(name, frontMedium, backMedium, iface);
		this.centre = centre.clone();
		this.dishNormal = dishCentreNormal;
		this.radiusOfCurv = radiusOfCurvature;
		this.dishDiameter = rimRadius*2;
		
		calcThetaCrit();
		calcCurvCentre();
		
	}
	
	protected void calcThetaCrit(){
		thetaCrit = FastMath.asin( dishDiameter / (2 * radiusOfCurv ) ); //angle away from axis of dish edge#
		boundSphereRadius = FastMath.sqrt( FastMath.pow(radiusOfCurv * (1 - FastMath.cos(thetaCrit)),2) + dishDiameter*dishDiameter/4 );
	}
	
	protected void calcCurvCentre(){
		curvCentre = new double[3];		
		for(int i=0;i<3;i++)
			curvCentre[i] = centre[i] + (radiusOfCurv * dishNormal[i]); //centre of cruvature
	}
	
	@Override
	public double[] getBoundarySphereCentre() { return centre; }

	@Override
	public double getBoundarySphereRadius() { return boundSphereRadius; }


	
	public void rotate(double point[], double matrix[][]){
		for(int i=0;i<3;i++){
			centre[i] -= point[i];
			curvCentre[i] -= point[i];
		}
			
		double newCentre[] = new double[3];
		double newCurvCentre[] = new double[3];
		double newNormal[] = new double[3];
		
		for(int i=0;i<3;i++)
			for(int j=0;j<3;j++){
				newCentre[i] += matrix[i][j] * centre[j];
				newCurvCentre[i] += matrix[i][j] * curvCentre[j];
				newNormal[i] += matrix[i][j] * dishNormal[j];
			}
				
		for(int i=0;i<3;i++){
			centre[i] = point[i] + newCentre[i];
			curvCentre[i] = point[i] + newCurvCentre[i];
		}
		
		dishNormal = newNormal;
	}
	
	@Override
	public void shift(double[] dX) {
		for(int i=0;i<3;i++){
			centre[i] += dX[i];
			curvCentre[i] += dX[i];
		}
	}

	@Override
	public boolean findEarlierIntersection(RaySegment ray, Intersection hit) {
	
		double a=1,b=0,c=0;
		
		for(int i=0;i<3;i++){
			b += 2 * (ray.startPos[i] - curvCentre[i]) * ray.dir[i];
			c += ((ray.startPos[i] - curvCentre[i])*(ray.startPos[i] - curvCentre[i]));
		}
		c -= radiusOfCurv * radiusOfCurv;

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

		//hit points and surface normals at those points
		for(int i=0;i<3;i++){
			X1[i] = ray.startPos[i] + s1 * ray.dir[i];
			X2[i] = ray.startPos[i] + s2 * ray.dir[i];
			N1[i] = curvCentre[i] - X1[i];
			N2[i] = curvCentre[i] - X2[i];
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
			
		double nnDotProd1=0, nnDotProd2=0;
		
		for(int i=0;i<3;i++){
			nnDotProd1 += dishNormal[i] * N1[i];
			nnDotProd2 += dishNormal[i] * N2[i];
		}
			
		//angle from dish normal axis 
		double theta1 = FastMath.acos( nnDotProd1 );
		double theta2 = FastMath.acos( nnDotProd2 );
		
		boolean lastSurfWasUs = (ray.startHit != null && ray.startHit.surface == this);
		double reHitTolerance =  lastSurfWasUs ? Tracer.reHitTolerance : 0;
	
		//which of those points are inside the dish's diameter? (& not backwards)
		boolean hit1 = (theta1 <= thetaCrit) & (s1 > reHitTolerance);
		boolean hit2 = (theta2 <= thetaCrit) & (s2 > reHitTolerance);
		
		
		if( hit1){ 
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
			}else{ //none hit the dish
				
			}	
		}	 
		return false; //didn't hit, or hit after it's length (i.e. something somewhere else hit first)
	}

	/** Getters for defining properties */
	public double[] getCentre() { return centre.clone(); }
	public double[] getDishNormal() { return dishNormal.clone(); }
	public double getRadiusOfCurv() { return radiusOfCurv; }
	public double getDishDiameter() { return dishDiameter; }
	
	/** Setters for defining properties */
	public void setCentre(double centre[]){ this.centre = centre; calcCurvCentre();	 }
	public void setDishNormal(double[] dishNormal) { this.dishNormal = dishNormal; calcCurvCentre();	}
	public void setRadiusOfCurv(double radiusOfCurv) { this.radiusOfCurv = radiusOfCurv; calcCurvCentre(); calcThetaCrit(); }
	public void setDishDiameter(double dishDiameter) { this.dishDiameter = dishDiameter; calcThetaCrit();	}

	/** Getters for calculated properties */
	public double[] getCurvCentre() { return curvCentre; }
	public double getThetaCrit() { return thetaCrit; }

	/** @return the depth of the surface, away from the x=0 plane at the given radius in (z,y) 
	 * in the 'dish' frame, which has x along the dish normal and the dish centre at (0,0,0) 
	 * @param r2	Squared radius away from dish central axis	 */
	protected double depthFromPlaneR2(double r2){
		return r2 / ( radiusOfCurv * (1.0 + FastMath.sqrt(1.0 - r2/radiusOfCurv/radiusOfCurv)) );		
	}

	/** @return the depth of the surface, away from the x=0 plane at the given radius in (z,y) 
	 * in the 'dish' frame, which has x along the dish normal and the dish centre at (0,0,0) 
	 * @param r		Radius away from dish central axis	 */
	public final double depthFromPlane(double r){
		return depthFromPlaneR2(r*r);
	}
	
	/** @return The position of the centre of the dish rim circle */ 
	public final double[] getRimCentre(){
		double rimRadius = dishDiameter / 2;
		double rimDepth = depthFromPlane(rimRadius);
		return Util.plus(centre, Util.mul(dishNormal, rimDepth));
	}
	
	public List<double[][]> draw() {
		int nSectors = 4 + (approxDrawQuality * 20 / 100); // min 4, max 24, 6 at def(10)
		int nRings = 2 + (approxDrawQuality * 20 / 100); // min 2, max 22, 4 at def(10)
		int nPointsInSector = 2 + (approxDrawQuality * 20 / 100); // min 2, max 22, 4 at def(10)
		
		ArrayList<double[][]> lines = new ArrayList<double[][]>();
		
		double u[] = Util.createPerp(dishNormal);
		double v[] = Util.cross(dishNormal, u);
		
		double rimRadius = dishDiameter / 2;
		
		for(int iR=0; iR < nRings; iR++){
			double r0 = iR*rimRadius / nRings;
			double r1 = (iR+1)*rimRadius / nRings;
			double xAS0 = depthFromPlaneR2(r0*r0);				
			double xAS1 = depthFromPlaneR2(r1*r1);				
			
			for(int iS=0; iS < nSectors; iS++){
				double phi0 = iS * 2 * Math.PI / nSectors;			
				double phi1 = (iS+1) * 2 * Math.PI / nSectors;
				
				double line[][] = new double[3][2*nPointsInSector+1];
				
				for(int iP = 0; iP < nPointsInSector; iP++){
					double phi = phi0 + (phi1 - phi0) * iP / (nPointsInSector - 1.0);

					for(int k=0; k < 3; k++)
						line[k][iP] = centre[k] + xAS0 * dishNormal[k]  
						                       + r0 * FastMath.cos(phi) * u[k]
		            						   + r0 * FastMath.sin(phi) * v[k];
					
					phi = phi1 - (phi1 - phi0) * iP / (nPointsInSector - 1.0);
					for(int k=0; k < 3; k++)
						line[k][nPointsInSector+iP] = centre[k] + xAS1 * dishNormal[k]  
						                       + r1 * FastMath.cos(phi) * u[k]
		            						   + r1 * FastMath.sin(phi) * v[k];

					for(int k=0; k < 3; k++)
						line[k][2*nPointsInSector] = centre[k] + xAS0 * dishNormal[k]  
								                       + r0 * FastMath.cos(phi0) * u[k]
		 		            						   + r0 * FastMath.sin(phi0) * v[k];
					
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
		r = 31 * r + Arrays.hashCode(dishNormal);
		t = Double.doubleToLongBits(dishDiameter);	r = 31 * r + (int) (t ^ (t >>> 32));
		t = Double.doubleToLongBits(radiusOfCurv); 	r = 31 * r + (int) (t ^ (t >>> 32));
		return r;
	}
	
	@Override
	public boolean surfaceGeometryEquals(Surface obj) {
		Dish other = ((Dish)obj);
		return Arrays.equals(centre, other.centre)			
					&& Arrays.equals(dishNormal, other.dishNormal)
					&& Double.doubleToLongBits(dishDiameter) == Double.doubleToLongBits(other.dishDiameter)			
					&& Double.doubleToLongBits(radiusOfCurv) == Double.doubleToLongBits(other.radiusOfCurv);
	}
	
}
