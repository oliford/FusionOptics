package fusionOptics.pointSpread;

import uk.co.oliford.jolu.OneLiners;

/** The idea here is to be able to generate a PSF for any 
 * point in a 3D volume, given a regular 3D grid of PSFs in
 * that volume.
 * 
 * It's done by interpolating the psfData() returned by the 
 * grid original PSFs.
 * 
 * @author oliford
 *
 */
public class PSFStatsSourceInterpolation {
	private boolean extrapolateMissingCorners = false; /** Experimental, works but is slow and doesn't help much */
	
	private double x0, x1, dx, y0, y1, dy, z0, z1, dz;
	private int nx, ny, nz;
	private PSFGrid psfGrid;
	
	public PSFStatsSourceInterpolation(PSFGrid psfGrid) {
		this.psfGrid = psfGrid;
		double gridDef[][] = (double[][])psfGrid.getGridDefinition();
		
		x0 = gridDef[0][0]; x1 = gridDef[0][1]; nx = (int)gridDef[0][2]; dx = (x1 - x0) / (nx-1.0);
		y0 = gridDef[1][0]; y1 = gridDef[1][1]; ny = (int)gridDef[1][2]; dy = (y1 - y0) / (ny-1.0);
		z0 = gridDef[2][0]; z1 = gridDef[2][1]; nz = (int)gridDef[2][2]; dz = (z1 - z0) / (nz-1.0);
	}
	
	public PointSpreadFunction getInterpolatedPSF(double x, double y, double z){
		if( x < x0 || x >= x1 || y < y0 || y >= y1 || z < z0 || z >= z1  ) 
			return null;
		
		int iX = (int)((x - x0) / dx);
		int iY = (int)((y - y0) / dy);
		int iZ = (int)((z - z0) / dz);
		
		PointSpreadFunction cornerPSFs[] = new PointSpreadFunction[]{
				getCornerPSF(iX, 	iY, 	iZ ),		// 000
				getCornerPSF(iX+1, 	iY, 	iZ ), 		// 001
				getCornerPSF(iX, 	iY+1, 	iZ ),  		// 010
				getCornerPSF(iX+1, 	iY+1, 	iZ ),  		// 011
				getCornerPSF(iX, 	iY, 	iZ+1 ),  	// 100
				getCornerPSF(iX+1, 	iY, 	iZ+1 ),  	// 101
				getCornerPSF(iX, 	iY+1, 	iZ+1 ),  	// 110
				getCornerPSF(iX+1, 	iY+1, 	iZ+1 )  	// 111
		};
		
		for(int i=0; i < cornerPSFs.length; i++){
			if(cornerPSFs[i] == null || cornerPSFs[i].isEmpty()){
				return null; //we can't deal with non-existing corners, because the others
							//will bias the image (x,y) and the edges will all fold to single points/lines
			}
		}

		double fx = ((x - x0) % dx) / dx;
		double fy = ((y - y0) % dy) / dy;
		double fz = ((z - z0) % dz) / dz;
		
		double coeffs[] = new double[]{
				(1-fz)	*(1-fy)	*(1-fx),			
				(1-fz)	*(1-fy)	*fx, 
				(1-fz)	*fy		*(1-fx), 
				(1-fz)	*fy		*fx, 
				fz		*(1-fy)	*(1-fx), 
				fz		*(1-fy)	*fx, 
				fz		*fy		*(1-fx), 
				fz		*fy		*fx,
		};
		
		PointSpreadFunction psf =  createLikePSF(cornerPSFs);
		
		psf.combine(cornerPSFs, coeffs);
		
		return psf;
		
	}

	/** Gets the PSF at a grid cell corner. It is doesn't exist or is invalid, it is created 
	 * from the average of the extrapolation of all neighbouring cell pairs. */
	private PointSpreadFunction getCornerPSF(int iX, int iY, int iZ) {
		
		PointSpreadFunction psf = psfGrid.get(iX, iY, iZ);
		
		//if it exists and is valid, just return it
		if(psf != null && !psf.isEmpty())
			return psf;
		
		if(!extrapolateMissingCorners)
			return null;
		
		//if(true)return null;
		//try to generate this point by extrapolating some others,
		//Look in each of the 6 directions away from this point for two neighbours.
		//If found, extrapolate from them.
		//Take average of all extrapolated values found;
		
		PointSpreadFunction extrapPSFs[] = new PointSpreadFunction[12];
		double coeffs[] = new double[12];
		
		if( !(
			fillPSFExtrapolationEntries(extrapPSFs, coeffs,  0,  iX, iY, iZ,   iX-1, iY, iZ,   iX-2, iY, iZ) |
			fillPSFExtrapolationEntries(extrapPSFs, coeffs,  2,  iX, iY, iZ,   iX+1, iY, iZ,   iX+2, iY, iZ) |
			fillPSFExtrapolationEntries(extrapPSFs, coeffs,  4,  iX, iY, iZ,   iX, iY-1, iZ,   iX, iY-2, iZ) |
			fillPSFExtrapolationEntries(extrapPSFs, coeffs,  6,  iX, iY, iZ,   iX, iY+1, iZ,   iX, iY+2, iZ) |
			fillPSFExtrapolationEntries(extrapPSFs, coeffs,  8,  iX, iY, iZ,   iX, iY, iZ-1,   iX, iY, iZ-2) |
			fillPSFExtrapolationEntries(extrapPSFs, coeffs, 10,  iX, iY, iZ,   iX, iY, iZ+1,   iX, iY, iZ+2) ) )
				return null; // no valid extrapolations found
		
		psf = createLikePSF(extrapPSFs);
		psf.combine(extrapPSFs, coeffs);
		
		//The distribution properties are extrapolated, but the intensity should be near-zero
		psf.multiplyIntensity(1e-20);
		
		return psf;
	}
	
	/** Fill in combination entries for extrapolation from the two given offsets of the given current node */
	private boolean fillPSFExtrapolationEntries(
			PointSpreadFunction psfs[],
			double coeffs[], 
			int idx, 
			int iX, int iY, int iZ, 
			int iXA, int iYA, int iZA,
			int iXB, int iYB, int iZB) {
		
		//check the request is in grid range
		if(iXA < 0 || iYA < 0 || iZA < 0 || 
				iXB < 0 || iYB < 0 || iZB < 0 ||
				iXA >= nx || iYA >= ny || iZA >= nz ||  
				iXB >= nx || iYB >= ny || iZB >= nz)
			return false;
		
		//get the grid PSFs
		PointSpreadFunction psfA =  psfGrid.get(iXA, iYA, iZA);
		PointSpreadFunction psfB =  psfGrid.get(iXB, iYB, iZB);
		
		//check they are valid
		if(psfA == null || psfB == null || psfA.isEmpty() || psfB.isEmpty())
			return false;
		
		double fx = iX - iXA;
		double fy = iY - iYA;
		double fz = iZ - iZA;

		double f = fx + fy + fz; //not mathematically correct, but only 1 will ever be non-zero here
		
		
		//linear extrapolation
		psfs[idx] = psfA;
		coeffs[idx] = (1.0 - f);

		psfs[idx+1] = psfB;
		coeffs[idx+1] = f;
		
		return true;
	}
	
	/** create a PSF of the same type as those in the given array */
	private PointSpreadFunction createLikePSF(PointSpreadFunction psfs[]) {
		PointSpreadFunction psf = null;
		try {
			for(int i=0; i < psfs.length; i++){
				if(psfs[i] != null){
					psf = psfs[i].getClass().newInstance();
					break;
				}
			}
			if(psf == null)
				return null;
		} catch (InstantiationException e) {
			throw new RuntimeException(e);
		} catch (IllegalAccessException e) {
			throw new RuntimeException(e);
		}
		
		return psf;
	}
}
