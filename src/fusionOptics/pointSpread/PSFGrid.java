package fusionOptics.pointSpread;

import java.util.Arrays;
import java.util.List;

import uk.co.oliford.cache.common.Cache;
import uk.co.oliford.cache.randomAccessCache.RACacheService;
import uk.co.oliford.jolu.BinaryMatrixWriter;

/** Handles saving and retrieval of PSFs on a regular 3D grid from the cache */
public class PSFGrid {
	private Cache psfCache = RACacheService.instance().getCache("optics.PSF");
	private String setName;
	
	private Class psfClass;
	
	private double x0, x1, dx;
	private double y0, y1, dy;
	private double z0, z1, dz;
	private int nx, ny, nz;
	
	public PSFGrid(String setName) {
		this.setName = setName;
		double gridDef[][] = (double[][]) psfCache.get(setName, "grid");
		if(gridDef == null)
			throw new RuntimeException("No PSF grid definition in cache set '"+setName+"'. Has the cache been generated?");
		setGridMembers(gridDef);
		psfClass = (Class) psfCache.get(setName, "class");
	}

	public PSFGrid(Cache psfCache, String setName) {
		this.psfCache = psfCache;
		this.setName = setName;
		setGridMembers((double[][]) psfCache.get(setName, "grid"));
		psfClass = (Class) psfCache.get(setName, "class");
	}
	
	public PSFGrid(String setName, double gridDef[][], Class psfClass) {
		this.setName = setName;
		setGridDefinition(gridDef);
		this.psfClass = psfClass;
		psfCache.put(setName, "class", psfClass);
	}
	
	public PSFGrid(Cache psfCache, String setName, double gridDef[][], Class psfClass) {
		this.psfCache = psfCache;
		this.setName = setName;
		setGridDefinition(gridDef);
		this.psfClass = psfClass;
		psfCache.put(setName, "class", psfClass);
	}

	private void setGridMembers(double gridDef[][]){
		x0 = gridDef[0][0]; x1 = gridDef[0][1]; nx = (int)gridDef[0][2]; dx = (x1 - x0) / (nx-1);
		y0 = gridDef[1][0]; y1 = gridDef[1][1]; ny = (int)gridDef[1][2]; dy = (y1 - y0) / (ny-1);
		z0 = gridDef[2][0]; z1 = gridDef[2][1]; nz = (int)gridDef[2][2]; dz = (z1 - z0) / (nz-1);
		if(nx <= 1){ dx=0; x0=(x0+x1)/2; x1=x0; }
		if(ny <= 1){ dy=0; y0=(y0+y1)/2; y1=y0; }
		if(nz <= 1){ dz=0; z0=(z0+z1)/2; z1=z0; }
	}
	
	/** Returns a list of all source positions in the PSF cache. */
	public double[][] getSourcePositions(){
		List<Object> keys = psfCache.getKeys(setName);
		
		double sourcePos[][] = new double[keys.size()][];
		int i=0;
		for(Object key : keys){
			if(key instanceof int[]){
				int idx[] = (int[])key;
				sourcePos[i] = gridPos(idx[0], idx[1], idx[2]);
				i++;
			}
		}
		return Arrays.copyOf(sourcePos, i);
	}

	public PointSpreadFunction get(int iX, int iY, int iZ) {
		return get(new int[]{ iX, iY, iZ });
	}
	
	public PointSpreadFunction get(int idx[]) {
		Object o = psfCache.get(setName, idx);
		
		if(o instanceof PointSpreadFunction) //old cache type
			return (PointSpreadFunction)o;
		
		//now we store the data
		double data[][] = (double[][]) o;
		if(data == null)
			return null;
		PointSpreadFunction psf;
		try {
			psf = (PointSpreadFunction) psfClass.newInstance();
		} catch (InstantiationException e) {
			throw new RuntimeException(e);
		} catch (IllegalAccessException e) {
			throw new RuntimeException(e);
		}
		psf.setMeanSourceDir(data[0]);
		psf.setMeanSourceUp(data[1]);
		psf.setCharacterisationData(data[2]);
		
		return psf; 
	}

	public PointSpreadFunction get(double[] sourcePos) {
		return get(new int[]{ (int)( (sourcePos[0]-x0)/dx + 0.5 ), 
					(int)( (sourcePos[1]-y0)/dy + 0.5 ), 
					(int)( (sourcePos[2]-z0)/dz + 0.5 ) });
	}
	
	public double[][] getGridDefinition() {
		return (double[][]) psfCache.get(setName, "grid");
	}
	
	public void setGridDefinition(double gridDef[][]) {
		psfCache.put(setName, "grid", gridDef);
		setGridMembers(gridDef);
	}

	/** Returns array of source positions and stats for all PSFs in the collection */
	public double[][] getAllData() {
		List<Object> keys = psfCache.getKeys(setName);
		
		double allData[][] = new double[keys.size()][];
		
		int i=0;
		for(Object key : keys){
			if(!(key instanceof int[]))
				continue;
			
			int idx[] = (int[])key;
			
			double sourcePos[] = gridPos(idx[0], idx[1], idx[2]);
			
			PointSpreadFunction psf = (PointSpreadFunction) get(idx);
			
			if(!psf.isEmpty()) {
				double data[] = psf.getCharacterisationData();
				double meanStartDir[] = psf.getMeanSourceDir();
				double meanStartUp[] = psf.getMeanSourceUp();
				
				allData[i] = new double[9 + data.length];
				System.arraycopy(sourcePos, 0, allData[i], 0, 3);
				System.arraycopy(meanStartDir, 0, allData[i], 3, 3);
				System.arraycopy(meanStartUp, 0, allData[i], 6, 3);
				System.arraycopy(data, 0, allData[i], 9, data.length);
			
				i++;
			}
		}
		
		return Arrays.copyOf(allData, i);
	}
	
	/** Returns array of source positions and stats for all PSFs in the collection */
	public void dumpAllData(String fileName) {
		BinaryMatrixWriter out = null; 
		List<Object> keys = psfCache.getKeys(setName);
		
		long t0 = System.currentTimeMillis();
		int i=0, j=0, n = keys.size();
		for(Object key : keys){
			if(!(key instanceof int[]))
				continue;
			
			int idx[] = (int[])key;
			
			double sourcePos[] = gridPos(idx[0], idx[1], idx[2]);
			
			PointSpreadFunction psf = (PointSpreadFunction) psfCache.get(setName, idx);
			
			if(!psf.isEmpty()) {
				double data[] = psf.getCharacterisationData();
				double meanStartDir[] = psf.getMeanSourceDir();
				double meanStartUp[] = psf.getMeanSourceUp();
				
				if(out == null)
					out = new BinaryMatrixWriter(fileName, 9 + data.length);
			
				out.writeRow(sourcePos, meanStartDir, meanStartUp, data);
				
				i++;
			}
			
			if((System.currentTimeMillis() - t0) > 1000){
				System.out.println(j + " / " + n + " " + (j*100.0 / n) + "%");
				t0 = System.currentTimeMillis();
			}
			j++;
		}
		
		if(out != null)
			out.close();
	}
	
	public void put(int iX, int iY, int iZ, PointSpreadFunction psf) {
		//if(!psf.isEmpty())
		//psfCache.put(setName, new int[]{ iX, iY, iZ }, psf); //direct in
		psfCache.put(setName, new int[]{ iX, iY, iZ }, new double[][]{ //store the data
				psf.getMeanSourceDir(),
				psf.getMeanSourceUp(),
				psf.getCharacterisationData() }); //*/
	}
	
	public void put(double[] sourcePos, PointSpreadFunction psf) {
		// TODO Auto-generated method stub
		put((int)( (sourcePos[0]-x0)/dx + 0.5 ), 
				(int)( (sourcePos[1]-y0)/dy + 0.5 ), 
				(int)( (sourcePos[2]-z0)/dz + 0.5 ),
				psf);
		
	}

	public final int getNX() { return nx; }
	public final int getNY() { return ny; }
	public final int getNZ() { return nz; }

	public final double[] gridPos(int iX, int iY, int iZ) {
		return new double[]{ x0+dx*iX, y0+dy*iY, z0+dz*iZ };
	}

	public String getSetName() { return setName; }
	
	public double getX0(){ return x0; }
	public double getX1(){ return x1; }
	public double getY0(){ return y0; }
	public double getY1(){ return y1; }
	public double getZ0(){ return z0; }
	public double getZ1(){ return z1; }
	
}
