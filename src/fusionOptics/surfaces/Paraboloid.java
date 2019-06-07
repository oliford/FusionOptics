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


/** Surface made from some cut of a elliptic Paraboloid.
 *
 *  Internal (n,r,u) frame with n along paraboloid axis, r to the right and u as up
 *  n = r^2 / cr^2  + u^2 / cu^2 
 *  
 */
public class Paraboloid extends Surface {
	
	protected double tip[];
	protected double normal[];
	protected double up[];
	protected double right[];
	
	/** curvature in y and z directions */
	protected double curvatureUp;
	protected double curvatureRight;
	
	/** Boundary in (u,r) space */
	protected double u0, u1, r0, r1;
	
	protected double[] boundSphereCentre;
	protected double boundSphereRadius;

	/** Attempt to make a parabolic mirror for a specific imaging setup
	 * @param name	Text name of surface
	 * @param centre	Centre of mirror
	 * @param focus		Point on which to focus the rays
	 * @param normal	Direction of incoming rays	 
	 * @param radius 	Maximum mirror radius
	 * @param iface
	 */
	public Paraboloid(String name, double centre[], double focus[], double normal[], double radius,
			Medium frontMedium, Medium backMedium, Interface iface) {
		super(name, frontMedium, backMedium, iface);
		
		double mf[] = Util.minus(focus, centre); //vector from mirror to focus
		
		double up[];
		if(Util.length(Util.cross(mf, normal)) == 0) { //on axis mirror, so we need to make up an up
			up = Util.createPerp(normal);
		}else {			
			//off-axis mirror, up is perp to that and incoming rays so we're off-axis in the 'right' direction
			up = Util.reNorm(Util.cross(mf, normal));
		}
		double right[] = Util.reNorm(Util.cross(normal, up));
		
		 
		double Mr = Util.dot(mf, right); //coordinate of mirror from tip along 'right' direction
		double MFn = Util.dot(mf, normal); //distance along of mirror behind focus (along 'normal' direction)
		
		double f1 = MFn/2 + FastMath.sqrt(MFn*MFn + Mr*Mr)/2;
		double f2 = MFn/2 - FastMath.sqrt(MFn*MFn + Mr*Mr)/2;
		
		double f;
		if(f1 <= 0 && f2 > 0)
			f = f2;
		else if(f2 <= 0 && f1 > 0)
			f = f1;
		else
			throw new RuntimeException("Two mirror positions possible: f = (" + f1 + ", " + f2 + ")");
		
		double c = FastMath.sqrt(4*f);
		
		double tip[] = Util.plus(focus, Util.mul(normal, -f));
		
		this.tip = tip;
		this.normal = normal;
		this.up = up;
		this.curvatureUp = c;
		this.curvatureRight = c;
		this.u0 = -radius;
		this.r0 = Mr - radius;
		this.u1 = radius;
		this.r1 = Mr + radius;
		
		calc();
	}
	
	/**
	 * @param name	Text name of surface
	 * @param tip	Tip of paraboloid
	 * @param normal	Normal of surface at tip
	 * @param up	Vector of 'up' direction (should be normal to tip normal)
	 * @param curvatureUp Coefficient of curvature in 'up' direction
	 * @param curvatureRight Coefficient of curvature in 'right' direction
	 * @param limits 	Surface boundary in (u,r) plane { u0, r0, u1, r1 }
	 * @param iface
	 */
	public Paraboloid(String name, double tip[], double normal[], double up[], double curvatureUp, double curvatureRight, double limits[], Interface iface) {
		this(name, tip, normal, up, curvatureUp, curvatureRight, limits, null, null, iface);
	}
	
	/**
	 * @param name	Text name of surface
	 * @param tip	Tip of paraboloid
	 * @param normal	Normal of surface at tip
	 * @param up	Vector of 'up' direction (should be normal to tip normal)
	 * @param curvatureUp Coefficient of curvature in 'up' direction
	 * @param curvatureRight Coefficient of curvature in 'right' direction
	 * @param limits 	Surface boundary in (u,r) plane { u0, r0, u1, r1 }
	 * @param iface
	 */
	public Paraboloid(String name, double tip[], double normal[], double up[], double curvatureUp, double curvatureRight, double limits[],
			Medium frontMedium, Medium backMedium, Interface iface) {
		super(name, frontMedium, backMedium, iface);
		this.tip = tip;
		this.normal = normal;
		this.up = up;
		this.curvatureUp = curvatureUp;
		this.curvatureRight = curvatureRight;
		this.u0 = limits[0];
		this.r0 = limits[1];
		this.u1 = limits[2];
		this.r1 = limits[3];
		
		calc();
	}
		
	@Override
	public double[] getBoundarySphereCentre() { return boundSphereCentre; }

	@Override
	public double getBoundarySphereRadius() { return boundSphereRadius; }


	private void calc() {
		this.right = Util.reNorm(Util.cross(normal, up));
		
		//if either origin is also covered, then we should also include the tip since the 
		// surface can be strongly convex inside the boundary 
		boolean includeTip = (u0 < 0 & u1 > 0) || (r0 < 0 & r1 > 0);
		
		//calc points of all 4 corners
		double p00[] = new double[3];
		double p01[] = new double[3];
		double p10[] = new double[3];
		double p11[] = new double[3];
		boundSphereCentre = new double[3];
		
		for(int i=0; i < 3; i++) {
			p00[i] = tip[i] + r0 * right[i] + u0 * up[i] + height(r0, u0) * normal[i];
			p01[i] = tip[i] + r1 * right[i] + u0 * up[i] + height(r1, u0) * normal[i];
			p10[i] = tip[i] + r0 * right[i] + u1 * up[i] + height(r0, u1) * normal[i];
			p11[i] = tip[i] + r1 * right[i] + u1 * up[i] + height(r1, u1) * normal[i];
			
			boundSphereCentre[i] = p00[i] + p01[i] + p10[i] + p11[i];
			if(includeTip)
				boundSphereCentre[i] = (boundSphereCentre[i] + 2*tip[i]) / 6;
			else
				boundSphereCentre[i] /= 4;
		}
		
		double br00 = Util.length(Util.minus(p00, boundSphereCentre));
		double br01 = Util.length(Util.minus(p01, boundSphereCentre));
		double br10 = Util.length(Util.minus(p10, boundSphereCentre));
		double br11 = Util.length(Util.minus(p11, boundSphereCentre));
		
		boundSphereRadius = Math.max(Math.max(br00, br01), Math.max(br10, br11));
		
				
		if(includeTip)
			boundSphereRadius = Math.max(boundSphereRadius, Util.length(Util.minus(tip, boundSphereCentre)));
	}

	
	public void rotate(double point[], double matrix[][]){
		for(int i=0;i<3;i++){
			tip[i] -= point[i];
		}
			
		double newCentre[] = new double[3];
		double newNormal[] = new double[3];
		double newUp[] = new double[3];
		
		for(int i=0;i<3;i++)
			for(int j=0;j<3;j++){
				newCentre[i] += matrix[i][j] * tip[j];
				newNormal[i] += matrix[i][j] * normal[j];
				newUp[i] += matrix[i][j] * up[j];
			}
				
		for(int i=0;i<3;i++){
			tip[i] = point[i] + newCentre[i];
		}

		normal = newNormal;
		up = newUp;
		right = Util.reNorm(Util.cross(normal, up));
	}
	
	@Override
	public void shift(double[] dX) {
		for(int i=0;i<3;i++){
			tip[i] += dX[i];
		}
	}
	
	/** The paraboloid equation
	 *  @return distance along normal axis at a given (r,u) position */
	private double height(double r, double u) {
		return FastMath.pow2(r / curvatureRight) + FastMath.pow2(u / curvatureUp);
	}

	@Override
	public boolean findEarlierIntersection(RaySegment ray, Intersection hit) {
	
		//convert ray start to (n,u,r) frame, scaled by curvature
		double ST[] = Util.minus(ray.startPos, tip);
		double nS = Util.dot(ST, normal);
		double rS = Util.dot(ST, right) / curvatureRight;
		double uS = Util.dot(ST, up) / curvatureUp;

		//convert ray direction, scaled by curvature
		double nL = Util.dot(ray.dir, normal);
		double rL = Util.dot(ray.dir, right) / curvatureRight;
		double uL = Util.dot(ray.dir, up) / curvatureUp;
		
		double a = rL*rL + uL*uL;
		double b = 2*rS*rL + 2*uS*uL - nL;
		double c = rS*rS + uS*uS - nS;

		if( (b*b) < (4*a*c) ) //no roots = no intersection		
			return false;
		
		double s1, s2;
		if(a < Tracer.reHitTolerance) { 
			//special case of 1 root or that the other
			// root is very very far out, in which case the primary
			// one can be found this way to within the tracer's 'general precision'
			s1 = -c/b;
			s2 = Double.POSITIVE_INFINITY;
		}else {
			double rtB24AC = FastMath.sqrt(b*b - 4*a*c);
			
			//two roots give two intersection points with sphere
			s1 = ( -b - rtB24AC ) / (2*a);
			s2 = ( -b + rtB24AC ) / (2*a);
		}
		
		double X1[] = new double[3];
		double X2[] = new double[3];
		double N1[] = new double[3];
		double N2[] = new double[3];
		
		//hit points
		for(int i=0;i<3;i++){
			X1[i] = ray.startPos[i] + s1 * ray.dir[i];
			X2[i] = ray.startPos[i] + s2 * ray.dir[i];
		}
		
		//hit points in (u,r) space
		double ST1[] = Util.minus(X1, tip);
		double ST2[] = Util.minus(X2, tip);
		
		double h1u = Util.dot(ST1, up);
		double h2u = Util.dot(ST2, up);
		double h1r = Util.dot(ST1, right);
		double h2r = Util.dot(ST2, right);
		
		//normals, from surface gradient
		for(int i=0; i < 3; i++) {
			N1[i] = -2*h1r/curvatureRight/curvatureRight * right[i] - 2*h1u/curvatureUp/curvatureUp * up[i] + 1 * normal[i]; 
			N2[i] = -2*h2r/curvatureRight/curvatureRight * right[i] - 2*h2u/curvatureUp/curvatureUp * up[i] + 1 * normal[i];
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
	
		//which of those points are inside the dish's diameter? (& not backwards)
		boolean hit1 = h1r >= r0 && h1r < r1 && h1u >= u0 && h1u < u1 && (s1 > reHitTolerance);
		boolean hit2 = h2r >= r0 && h2r < r1 && h2u >= u0 && h2u < u1 && (s2 > reHitTolerance);
		
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
	public double[] getTip() { return tip.clone(); }
	public double[] getNormal() { return normal.clone(); }
	public double[] getUp() { return up.clone(); }
	public double getCurvatureUp() { return curvatureUp; }
	public double getCurvatureRight() { return curvatureRight; }
	
	@Override
	public double[] getCentre() { return tip; }
	
	/** @return The position of the Paraboloid's focus from parallel light directed at it on axis, from average of curvatures */
	public double[] getFocus() { 
		double curv = (curvatureRight + curvatureUp) / 2;
		return Util.plus(tip, Util.mul(normal, FastMath.pow2(curv) / 4)); 
	}
	
	
	/** Setters for defining properties */
	public void setTip(double tip[]){ this.tip = tip; calc(); }
	public void setNormal(double[] normal) { this.normal = normal; calc(); }
	public void setUp(double[] up) { this.up = up; calc(); }
	public void setCurvatureUp(double curvatureUp) { this.curvatureUp = curvatureUp; calc(); }
	public void setCurvatureRight(double curvatureRight) { this.curvatureRight = curvatureRight; calc(); }
	
	
	public List<double[][]> draw() {
		int n = 2 * approxDrawQuality;
		double dU = (u1 - u0) / n;
		double dR = (r1 - r0) / n;
		
		ArrayList<double[][]> lines = new ArrayList<double[][]>();
		
		for(int iU=0; iU < n; iU++){
			double uA = u0 + iU * dU; 
			double uB = u0 + (iU+1) * dU; 
			
			for(int iR=0; iR < n; iR++){
				double rA = r0 + iR * dR; 
				double rB = r0 + (iR+1) * dR;
				
				// rectangle AA - AB - BB - BA
				double r[] = { rA, rA, rB, rB, rA };
				double u[] = { uA, uB, uB, uA, uA };
				
				double line[][] = new double[3][r.length];
				for(int k=0; k < 3; k++) {
					for(int j=0; j < r.length; j++) {
						line[k][j] = tip[k] + r[j] * right[k] + u[j] * up[k] + height(r[j], u[j]) * normal[k];
					}
				}
				
				lines.add(line);
			}	
		}
		
		return lines;
		
	}

	@Override
	public int surfaceGeometryHashCode() {
		long t; int r = 1;
		r = 31 * r + Arrays.hashCode(tip);
		r = 31 * r + Arrays.hashCode(normal);
		r = 31 * r + Arrays.hashCode(up);
		t = Double.doubleToLongBits(curvatureUp);	r = 31 * r + (int) (t ^ (t >>> 32));
		t = Double.doubleToLongBits(curvatureRight);	r = 31 * r + (int) (t ^ (t >>> 32));
		t = Double.doubleToLongBits(u0);	r = 31 * r + (int) (t ^ (t >>> 32));
		t = Double.doubleToLongBits(u1);	r = 31 * r + (int) (t ^ (t >>> 32));
		t = Double.doubleToLongBits(r0);	r = 31 * r + (int) (t ^ (t >>> 32));
		t = Double.doubleToLongBits(r1);	r = 31 * r + (int) (t ^ (t >>> 32));
		return r;
	}
	
	@Override
	public boolean surfaceGeometryEquals(Surface obj) {
		Paraboloid other = ((Paraboloid)obj);
		return Arrays.equals(tip, other.tip)			
				&& Arrays.equals(normal, other.normal)
				&& Arrays.equals(up, other.up)
				&& Double.doubleToLongBits(curvatureUp) == Double.doubleToLongBits(other.curvatureUp)
				&& Double.doubleToLongBits(curvatureRight) == Double.doubleToLongBits(other.curvatureRight)
				&& Double.doubleToLongBits(u0) == Double.doubleToLongBits(other.u0)	
				&& Double.doubleToLongBits(u1) == Double.doubleToLongBits(other.u1)	
				&& Double.doubleToLongBits(r0) == Double.doubleToLongBits(other.r0)	
				&& Double.doubleToLongBits(r1) == Double.doubleToLongBits(other.r1);
	}

}
