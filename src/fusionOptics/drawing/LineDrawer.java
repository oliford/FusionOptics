package fusionOptics.drawing;

import net.jafama.FastMath;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import fusionOptics.Util;
import fusionOptics.types.Element;
import fusionOptics.types.Optic;
import fusionOptics.types.RaySegment;
import fusionOptics.types.Surface;



/** Base class and common code for simple ray drawing 
 * to formats which just take a series of lines (like SVG).
 * 
 * The base provides some code for drawRay(), which splits things 
 *  
 */ 
public abstract class LineDrawer implements RayDrawer {
	private boolean drawWholePaths;
	private double smallLineLength;
	private boolean drawPolarisationFrames = false;
	private boolean drawIntersectionNormals = false;
	private int nSkipRays = 0;
	private int nSkipped = 0;
	
	public LineDrawer(boolean drawWholePaths, double smallLineLength) {
		this.drawWholePaths = drawWholePaths;
		this.smallLineLength = smallLineLength;
	}
	
	@Override
	public void setDrawPolarisationFrames(boolean enable) { this.drawPolarisationFrames = enable; }
	@Override
	public void setDrawIntersectionNormals(boolean enable) { this.drawIntersectionNormals = enable; }

	/** Length of 'a small line'. Used to draw the polarisation states, and to draw rays going to infinity. */
	public void setSmallLineLength(double smallLineLength) { this.smallLineLength = smallLineLength; }
	
	@Override
	public void setSkipRays(int nSkipRays){ this.nSkipRays = nSkipRays; }
	
	@Override
	public void drawRay(RaySegment ray) {
		if(nSkipped < nSkipRays){ nSkipped++; return; } else { nSkipped = 0; }
		drawRay(ray, new LinkedList<double[]>());
	}
		
	public void drawRay(RaySegment ray, LinkedList<double[]> lineData) {
		
		//add the start of this segement
		addPoint(lineData, ray.startPos, ray);
		
		while(ray.endHit != null && !ray.endHit.isEnd()) { 
			//while there are outgoing rays...

			//get them, with the strongest first
			List<RaySegment> nextSegs = ray.endHit.getSortedOutgoingRays();

			//for each new branch... (every outgoing ray excluding the strongest) 
			for(int i=1; i < nextSegs.size(); i++){
				LinkedList<double[]> newLineData;
				
				if(drawWholePaths){
					//if drawing whole paths, we need to start the new ray with the path already collected
					newLineData = new LinkedList<double[]>(lineData);
				}else{
					//otherwise the new ray starts from here
					newLineData = new LinkedList<double[]>();
					addPoint(lineData, ray.startPos, ray);
				}
					
				drawRay(nextSegs.get(i), newLineData);	
			}
			
			//continue with strongest
			ray = nextSegs.get(0);
			addPoint(lineData, ray.startPos, ray);
		}

		//if the segement doesn't end (i.e. isn't absorped but goes out to infinity)
		//then draw a small line
		if(ray.endHit != null){
			addPoint(lineData, ray.endHit.pos, ray);
		}else{
			lineData.add(new double[]{ 
					ray.startPos[0] + smallLineLength * ray.dir[0],
					ray.startPos[1] + smallLineLength * ray.dir[1],
					ray.startPos[2] + smallLineLength * ray.dir[2],
				});
		}
		
		drawSinglePath(lineData);
	}
	
	private final void addPoint(LinkedList<double[]> lineData, double pos[], RaySegment ray){
		lineData.add(new double[]{
				pos[0], pos[1], pos[2], ray.startIntensity(), (ray.endHit != null) ? ray.endIntensity() : ray.startIntensity()
		});
		if(drawIntersectionNormals && ray.startHit != null && ray.startHit.normal != null){
			drawSinglePath(new double[][]{
					{ pos[0] - 0.5*smallLineLength * ray.startHit.normal[0], pos[0] + 1.5*smallLineLength * ray.startHit.normal[0] },
					{ pos[1] - 0.5*smallLineLength * ray.startHit.normal[1], pos[1] + 1.5*smallLineLength * ray.startHit.normal[1] },
					{ pos[2] - 0.5*smallLineLength * ray.startHit.normal[2], pos[2] + 1.5*smallLineLength * ray.startHit.normal[2] },
			});
		}
		if(drawPolarisationFrames){
			throw new UnsupportedOperationException();
		
		/*	double right[] = Util.cross(ray.dir, ray.up);
			double linPol = FastMath.sqrt(ray.stokes[1]*ray.stokes[1] + ray.stokes[2]*ray.stokes[2]);
			
			double cos2q = ray.stokes[1] / linPol;
			double sin2q = ray.stokes[2] / linPol;
			
			//double q = Math.atan2(sin2q, cos2q) / 2;
			//double cosq = Math.cos(q);
			//double sinq = Math.sin(q);
			
			double cosq = FastMath.sqrt((cos2q + 1) / 2);
			double sinq = (cosq == 0) ? FastMath.sqrt((1 - cos2q) / 2) : sin2q  / (2 * cosq);
			
			
			lineData.add(new double[]{
					pos[0] + smallLineLength * (cosq * ray.up[0] + sinq * right[0]),
					pos[1] + smallLineLength * (cosq * ray.up[1] + sinq * right[1]),
					pos[2] + smallLineLength * (cosq * ray.up[2] + sinq * right[2]),
					ray.stokes[0]
			});
			lineData.add(pos);*/
						
			
		}
	}
	
	public abstract void drawSinglePath(List<double[]> path);
	
	public abstract void drawSinglePath(double path[][]);
	
	
	/** Base version of this just draws each line of an optic's drawing independently */
	public void drawElement(Element elem){
	
		if(elem instanceof Optic){
			Optic optic = (Optic)elem;
			startGroup("optic_"+optic.getName());
			
			for(Surface surface : optic.getSurfaces()) {
				drawElement(surface);
			}
			for(Optic subOptic : optic.getSubOptics()) {
				drawElement(subOptic);
			}
			endGroup();	
		}else{ //surface
			Surface surface = ((Surface)elem);
			startGroup("surface_"+surface.getName());
			for(double line[][] : surface.draw())				
				drawSinglePath(line);
			endGroup();
		}
	}
	
	public abstract void destroy();

	
}
