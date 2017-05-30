package fusionOptics.drawing;

import fusionOptics.types.Element;
import fusionOptics.types.RaySegment;

/** Interface for anything which accepts optics drawings and ray paths, and outputs to somewhere */ 
public interface RayDrawer {
	
	public void setDrawPolarisationFrames(boolean enable);
	
	public void setDrawIntersectionNormals(boolean enable);
	
	public void drawRay(RaySegment ray);
	
	/** Draw elements, optics and surfaces */
	public void drawElement(Element elem);
		
	public void setSkipRays(int nSkipRays);
	
	public void startGroup(String name);
	
	public void endGroup();
	
	public void destroy();
	
}
