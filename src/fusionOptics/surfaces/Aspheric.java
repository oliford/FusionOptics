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

import java.util.Arrays;

import fusionOptics.Util;
import fusionOptics.types.Interface;
import fusionOptics.types.Intersection;
import fusionOptics.types.Medium;
import fusionOptics.types.RaySegment;
import fusionOptics.types.Surface;

import net.jafama.FastMath;
import algorithmrepository.Algorithms;


/** Aspherical Surface - Small deformation of a part-sphere (Dish) 
 *   extends 'dish' class which provides the first order estimate of the ray intersection. 
 * 
 * Wikipedia [ http://en.wikipedia.org/wiki/Aspheric_lens ] has the definition:
 *  z(r) = r^2 / R*(1 + sqrt(1 - (1 + conicConst)*(r^2/R^2))) * sum_i[ alpha_i r^(2i) ]
 *  
 *  With the surface mostly in the (x,y) plane and intersected x=y=z=0
 * 
 * There is a very good chapter of a book here 
 * [ http://www.globalspec.com/reference/13858/121073/chapter-11-optical-systems-section-3-extract-ray-tracing ]
 * 
 * Which has the same definition, but with the conic constant fixed at zero.
 * 
 * The second link gives the main iterative proceedure used here to quickly find the surface to within the
 * specific numerical accurcy. It that fails (and it does quite often), a more brute force
 * section search is performed around the initial (sphere) guess. This is more robust but is slower. * 
 * 
 * It is unfortunately still possible for the ray intersection code to fail near the edge, since ray can hit the 
 * real surface without ever hitting the approximate surface. This is worse than just changing the effective size since 
 * only rays hitting the very edge at a glancing angle will fail while those directed straight on are OK.
 * 
 * If you put very strong polyCoeffs in, so that the edge of the disc is a long way away from the sphere, then it will fail.
 * Check the surfaceFailures counts! Warnings will also be printed to stderr.
 * 
 * There's also some interseting stuff here, which I've not read:
 * 
 * [ Gyeong-Il Kweon "Aspherical Lens Design by Using a Numerical Analysis"
 *   Journal of the Korean Physical Society, Vol. 51, No. 1, July 2007, pp. 93∼103
 *    http://www.google.co.uk/url?sa=t&rct=j&q=aspherical%20camera%20lens%20design%20radius&source=web&cd=6&sqi=2&ved=0CHIQFjAF&url=http%3A%2F%2Fwww.nanophotonics.kr%2Fsrc_doc%2FJ0707JKPS.pdf&ei=msFcT-2bDor0sgbq_Mz-Cw&usg=AFQjCNGUp54v5EsLJ4rSOawgUUDVNIS_Cw ]
 * 
 */

public class Aspheric extends Dish {
	private double surfaceFindingTolerance = 1e-8;	
	private int maxSurfaceFindingIterations = 500;
	
	/** Conic constant, from wikipedia's definition.
	 * It is implemented in depthFromPlane(), so the drawing definitely works with nonzero values.
	 * The iterative surface finding is tested using the correct formula so should work anyway, if a little inefficiently.
	 * The backup surface finding will work correctly.
	 * The biggest problem is that the calculation of the surface normal doesn't include it, fix this and everything will work.   
	 **/
	private double conicConstant = 0;
	/** Polynomial coefficients for:  z += a[0] r^2 + a[1] r^4 + a[2] r^6 + ... */
	private double polyCoeffs[];
	
	public int surfaceFailures = 0;
	public int surfaceBruteForces = 0;
	
	/**
	 * @param name				Element name.
	 * @param centre			The point on the surface, that is an equal distance from every point on the edge (rim). 
	 * @param dishCentreNormal	Unit vector that points towards the curvature centre, from the dish centre.
	 * @param radiusOfCurvature 
	 * @param rimRadius			Radius of dish rim.
	 * @param conicConst		Deformation of form of sphere.
	 * @param polyCoeffs		Coefficients for deformation from sphere away from centre:  a[0] r^2 + a[1] r^4 + a[2] r^6 + ...
	 * @param interfaceType
	 */
	public Aspheric(String name, double centre[], double dishCentreNormal[], double radiusOfCurvature, double rimRadius, 
			double conicConst, double polyCoeffs[], Interface iface) {
		super(name, centre, dishCentreNormal, radiusOfCurvature, rimRadius, iface);

		this.conicConstant = conicConst;
		this.polyCoeffs = polyCoeffs;
	}
	
	/**
	 * @param name				Element name.
	 * @param centre			The point on the surface, that is an equal distance from every point on the edge (rim). 
	 * @param dishCentreNormal	Unit vector that points towards the curvature centre, from the dish centre. i.e. points 'inside' the dish.
	 * @param radiusOfCurvature 
	 * @param rimRadius			Radius of dish rim.
	 * @param conicConst		Deformation of form of sphere. NOT SUPPORTED PROPERLY
	 * @param polyCoeffs		Coefficients for deformation from sphere away from centre:  a[0] r^2 + a[1] r^4 + a[2] r^6 + ...
	 * @param frontMedium		Medium on the 'inside' of the dish, that the normal points into.
	 * @param backMedium		Medium on the 'outside' of the dish, that the normal points away from.
	 * @param interfaceType
	 */
	public Aspheric(String name, double centre[], double dishCentreNormal[], double radiusOfCurvature, double rimRadius,
			double conicConst, double polyCoeffs[],
			Medium frontMedium, Medium backMedium, Interface iface) {
		super(name, centre, dishCentreNormal, radiusOfCurvature, rimRadius, frontMedium, backMedium, iface);
		
		this.conicConstant = conicConst;
		this.polyCoeffs = polyCoeffs;		
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
		
		dishNormal = Util.reNorm(newNormal);
	}
	
	@Override
	public void shift(double[] dX) {
		for(int i=0;i<3;i++){
			centre[i] += dX[i];
			curvCentre[i] += dX[i];
		}
	}
	
	@Override
	/** @return the depth of the surface, away from the x=0 plane at the given radius in (z,y) 
	 * in the 'dish' frame, which has x along the dish normal and the dish centre at (0,0,0) 
	 * @param r2	Squared radius away from dish central axis	 */
	protected double depthFromPlaneR2(double r2){
		double x = r2 / ( radiusOfCurv * (1.0 + FastMath.sqrt(1.0 - (1.0 + conicConstant)*r2/radiusOfCurv/radiusOfCurv)) );
		double rPow = r2;
		for(int i=0; i < polyCoeffs.length; i++){
			x += polyCoeffs[i] * rPow;
			rPow *= r2;
		}
		return x;
	}
	
	@Override
	public boolean findEarlierIntersection(RaySegment ray, Intersection hit) {

//		if(conicConstant == 0 && Algorithms.sum(polyCoeffs) == 0){ //check that zero behaviour matches Dish
//			return super.findEarlierIntersection(ray, hit);
//		}
		
		double oldRayLength = ray.length; //store that, because Dish.findEarlierIntersection() may overwrite it
		Intersection sphereHit = new Intersection();
		
		//Let Dish see if it hit the approximating spherical surface
		if(!super.findEarlierIntersection(ray, sphereHit))
			return false; //didn't hit the sphere
				/*  KNOWN BUG: It's possible that it does actually here the aspheric, bu
				 * didn't hit the sphere, if the apsheric is outside the sphere near the rim.
				 * It will also mess up if the ray length just catches the asphere but not the sphere. */
	
		/* create the dish frame, as in
		 * [http://beta.globalspec.com/reference/13858/121073/chapter-11-optical-systems-section-3-extract-ray-tracing]
		 * with x along the normal and (y,z) in the approximate plane of the dish and the dish centre at (0,0,0) */
		double yVec[] = Util.createPerp(dishNormal);
		double zVec[] = Util.cross(dishNormal, yVec);
		
		//direction of incident ray in the dish frame
		double X = Util.dot(ray.dir, dishNormal);
		double Y = Util.dot(ray.dir, yVec);
		double Z = Util.dot(ray.dir, zVec);

		//convert first order estimate of hit point (sphere) into that frame 
		double centreToSphereHit[] = Util.minus(sphereHit.pos, centre);
		double x0 = Util.dot(centreToSphereHit, dishNormal);
		double y0 = Util.dot(centreToSphereHit, yVec);
		double z0 = Util.dot(centreToSphereHit, zVec);
		
		//start of ray in dish coords
		double centreToRayStart[] = Util.minus(ray.startPos, centre);
		double rayStartDC[] = new double[]{
			Util.dot(centreToRayStart, dishNormal),
			Util.dot(centreToRayStart, yVec),
			Util.dot(centreToRayStart, zVec)
		};
			
		double x=x0, y=y0, z=z0;
		
		double C = 1 / radiusOfCurv;
		
		//radius in plane perp to dish axis
		double r2 = y*y + z*z;
		
		//calc aspheric depth at current r
		double xAS = depthFromPlaneR2(r2);
		
		double l0=0, m0=0, n0=0, G0=0;
		int i;
		for(i=0; i < maxSurfaceFindingIterations; i++){ //iterative improvement...
			
			//calc (l0,m0,n0) = the normal to surface at (xAS,y,z)
			//As far as I can tell, the conic constant gets added to the maths only in l0, and doesn't effect the rest of the differentials
			// Although I've not actually checked this
			l0 = FastMath.sqrt(1.0 - (1.0 + conicConstant)*C*C*r2);
			
			m0 = -y * C;
			n0 = -z * C;
			double rPow = 1.0;
			for(int j=0; j < polyCoeffs.length; j++){
				m0 -= y * l0 * (2*j+2) * polyCoeffs[j] * rPow;
				n0 -= z * l0 * (2*j+2) * polyCoeffs[j] * rPow;
				rPow *= r2;
			}
			
			//The distance along the ray to move
			G0 = l0 * (xAS - x) / (X*l0 + Y*m0 + Z*n0);
			
			x += G0 * X;
			y += G0 * Y;
			x += G0 * Z;
		
			//that calc gets us off the ray sometimes
			//so move to nearest point on the ray
			double rayVecDC[] = Util.minus(new double[]{x,y,z}, rayStartDC);
			double d= Util.dot(new double[]{X,Y,Z}, rayVecDC);
			x = rayStartDC[0] + d * X;
			y = rayStartDC[1] + d * Y;
			z = rayStartDC[2] + d * Z;
			
			r2 = y*y + z*z;
			xAS = depthFromPlaneR2(r2);
			
			if(FastMath.abs(xAS - x) < surfaceFindingTolerance) //we're now near enough,
				break;
		
		}
		
		if(i >= maxSurfaceFindingIterations){
			// If that method fails, we can do a more brute force attempt by sectioning the 
			//ray around the inital sphere estimate and just picking the closest approach to
			//the surface
			surfaceBruteForces++;
			
			double pos[] = findIntersectionBySplitting(X,Y,Z, x0,y0,z0);
			x = pos[0];
			y = pos[1];
			z = pos[2];
					
			r2 = y*y + z*z;
			
			///still need to calc the normal
			l0 = FastMath.sqrt(1.0 - (1.0 + conicConstant)*C*C*r2);
			m0 = -y * C;
			n0 = -z * C;
			double rPow = 1.0;
			for(int j=0; j < polyCoeffs.length; j++){
				m0 -= y * l0 * (2*j+2) * polyCoeffs[j] * rPow;
				n0 -= z * l0 * (2*j+2) * polyCoeffs[j] * rPow;
				rPow *= r2;
			}	
			
			xAS = depthFromPlaneR2(r2);
			
			if(FastMath.abs(xAS - x) > surfaceFindingTolerance){ //still didn't work
				System.err.println("WARNING: Aspheric.findEarlierIntersection() didn't find apsheric surface after "
										+ maxSurfaceFindingIterations + " iterations or by bisection.");
				surfaceFailures++;
			}
			
		}
				
		double distToSurace = FastMath.abs(x-xAS);
		
		if((y*y + z*z) > dishDiameter*dishDiameter/4){
			ray.length = oldRayLength;
			return false; //converged radius outside rim
		}
		
		//project position back into world coords
		double pos[] = new double[]{ 
				centre[0] + x*dishNormal[0] + y*yVec[0] + z*zVec[0],
				centre[1] + x*dishNormal[1] + y*yVec[1] + z*zVec[1],
				centre[2] + x*dishNormal[2] + y*yVec[2] + z*zVec[2],
			};
		
		double distToLine = Algorithms.distLinePoint(ray.startPos, Util.plus(ray.startPos, ray.dir), pos);
		
		//System.out.println(i+"\t"+distToSurace+"\t"+distToLine);
		
		ray.length = Util.dot(ray.dir,Util.minus(pos, ray.startPos)); //TODO: lazy, there's probably a more efficient way
		
		//if(ray.length < Tracer.reHitTolerance){
		//because of the numerics involved, we have to be a bit more rough with this one
		if(ray.startHit != null && ray.startHit.surface == this && ray.length < 10*surfaceFindingTolerance){
			ray.length = oldRayLength;
			return false; //ray too short - we're probably re-hitting this surface
		}
		
		//project normal into world coords
		hit.normal = new double[]{ 
				l0*dishNormal[0] + m0*yVec[0] + n0*zVec[0],
				l0*dishNormal[1] + m0*yVec[1] + n0*zVec[1],
				l0*dishNormal[2] + m0*yVec[2] + n0*zVec[2],
			};
		
		hit.normal = Util.reNorm(hit.normal); 
		
		hit.surface = this;
		hit.pos = pos;
		
		return true;
	
	}

	/** More brute force method of intersection finding, by sectioning of line to search for min in xAS 
	 * @param X,Y,Z 	Direction of ray in dish frame
	 * @param x0,y0,z0	Initial estimate of intersection (sphere) in dish frame
	 * @return Best estimate of intersection [x/y/z] in dish frame 
	 */
	private double[] findIntersectionBySplitting(double X, double Y, double Z, 
												double x0, double y0, double z0){
		
		int nIts = 50;
		int nSplit = 10;
		
		double d0 = -radiusOfCurv/5; // distances along line (x,y,z) + i(X,Y,Z)
		double d1 = +radiusOfCurv/5;
		
		double minXDiff = Double.POSITIVE_INFINITY;
		double dAtMin = Double.NaN;
		
		for(int i=0; i < nIts; i++){
			double dd = (d1 - d0) / (nSplit - 1);
				
			for(int j=0; j < nSplit; j++){
				double d = d0 + j * dd;
				
				double x = x0 + d * X;
				double y = y0 + d * Y;
				double z = z0 + d * Z;
				
				double r2 = y*y + z*z;
				
				double xAS = FastMath.abs(depthFromPlaneR2(r2) - x);
				
				if(xAS < minXDiff){
					minXDiff = xAS;
					dAtMin = d;
				}				
			}
			
			d0 = dAtMin - dd;
			d1 = dAtMin + dd;
			
			if(minXDiff < surfaceFindingTolerance)
				break; //close enough

		}
		
		return new double[]{
				x0 + dAtMin * X, 
				y0 + dAtMin * Y, 
				z0 + dAtMin * Z 
			};
	}
	
	public double[] getPolyCoeffs(){ return polyCoeffs; }
	public double getConicConstant(){ return conicConstant; }
	
	/** Sets the polynomial coefficients. { c(R^2), c(R^4), c(R^6) ... } */
	public void setPolyCoeffs(double[] polyCoeffs){ this.polyCoeffs = polyCoeffs; }
	
	public void setConicConstant(double conicConstant){ this.conicConstant = conicConstant; }
	
	@Override
	public int surfaceGeometryHashCode() {
		long temp = Double.doubleToLongBits(conicConstant);
		return (super.surfaceGeometryHashCode() * 31 + (int) (temp ^ (temp >>> 32)))
					* 31 + Arrays.hashCode(polyCoeffs);		
	}

	@Override
	public boolean surfaceGeometryEquals(Surface obj) {
		Aspheric other = (Aspheric) obj;
		return super.surfaceGeometryEquals(obj)
				&& Double.doubleToLongBits(conicConstant) == Double.doubleToLongBits(other.conicConstant)
				&& Arrays.equals(polyCoeffs, other.polyCoeffs);
	}

	/** Rescales the polyCoeffs to match a change in scale of the radius, curvatureRadius etc
	 * E.g if the coefficients were originally for dimensions in mm, but the lens object here is in m
	 * then call rescaleCoeffs(C, 1e-3) */
	public static double[] rescaleCoeffs(double[] polycoeffs, double scale) {
		double scaledCoeffs[] = new double[polycoeffs.length];
			
		double nextScale = 1.0 / scale;
		for(int i=0; i < polycoeffs.length; i++){
			scaledCoeffs[i] = polycoeffs[i] * nextScale;
			nextScale /= scale * scale;
		}
		return scaledCoeffs;
	}

	public void setSurfaceFindingTolerance(double surfaceFindingTolerance, int maxSurfaceFindingIterations) {
		this.surfaceFindingTolerance = surfaceFindingTolerance;
		this.maxSurfaceFindingIterations = maxSurfaceFindingIterations;
	}
}
