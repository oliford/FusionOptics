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
package fusionOptics.drawing;

import net.jafama.FastMath;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import fusionOptics.Util;
import fusionOptics.interfaces.Absorber;
import fusionOptics.interfaces.Reflector;
import fusionOptics.types.Element;
import fusionOptics.types.Optic;
import fusionOptics.types.Pol;
import fusionOptics.types.RaySegment;
import fusionOptics.types.Surface;
import uk.co.oliford.jolu.OneLiners;



/** Draws rays into a VRML file.
 * 
 *  Builds the polarisation and intensity info into the object names.
 *  Draws polarisation information if requested.
 * 
 *  This doesn't use the LineDrawer base, because we draw as we go
 *  rather than wanting the precompiled List<double[]> lineData.
 *  
 */ 
public class VRMLDrawer implements RayDrawer {
	private boolean drawPolarisationFrames = false;
	private boolean drawIntersectionNormals = false; 
	private double smallLineLength;
	private boolean drawOnlyStrongest = false;
	
	private FileOutputStream fos;
	private PrintWriter writer;
	private String fileName;
	private int nRays = 0;
	private double currentColour[] = new double[]{ 0, 0, 0 };
	private double rayStartIntensity = 0;
	private int nSkipRays = 0;
	private int nSkipped = 0;
	
	private double globalRotMat[][] = null;
	
	public DecimalFormat fmt = new DecimalFormat("00.000000");
	private static final DecimalFormat polFmt = new DecimalFormat("#.##");
	
	public VRMLDrawer(String fileName, double smallLineLength) {
		this(fileName);
		this.smallLineLength = smallLineLength;
	}
	
	/**
	 * @param fileName VRML files to write to.
	 * @param smallLineLength Length of 'a small line'. Used to draw the polarisation states, and to draw rays going to infinity.
	 */
	public VRMLDrawer(String fileName) {
		this.fileName = fileName;
		DecimalFormatSymbols dfs = polFmt.getDecimalFormatSymbols();
		dfs.setNaN("NaN");
		dfs.setInfinity("Inf");
		
		dfs.setDecimalSeparator('d');
		polFmt.setDecimalFormatSymbols(dfs);
		
		dfs.setDecimalSeparator('.');
		fmt.setDecimalFormatSymbols(dfs);
				
		OneLiners.makePath(fileName);
		try {
			fos = new FileOutputStream(fileName);
			writer = new PrintWriter(fos);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
		writer.println("#VRML V2.0 utf8\n");
		
		addMat("rayLine", "0 0 0" , 0);
		addMat("polEllipPolarised0", "1 0 0" , 0);
		addMat("polEllipPolarised1", "0 0 1" , 0);
		addMat("polEllipPolarised2", "0 1 0" , 0);
		addMat("polEllipPolarised3", "1 0 1" , 0);
		addMat("polEllipIntensity", "0.5 0.5 0.5" , 0);
		addMat("absorbingSurface", "1 1 0" , 0.7);
		addMat("reflectiveSurface", "0 1 1" , 0.7);
		addMat("generalSurface", "0 1 0" , 0.9);
		
		writer.println("DEF rayTracing Group { children [");
		
	}
	
	/** Add a material of the given name, to be used for drawing things later */
	public void addMat(String name, String diffuseColor, double transparency){
		writer.println(
				"DEF "+name+" Appearance { \n" + 
				"     material Material { \n" +
				"      diffuseColor "+diffuseColor+" \n" +
				"      ambientIntensity 1.0 \n" +
				"      specularColor 0 0 0 \n" +
				"      shininess 0.0 \n" +
				"      transparency "+fmt.format(transparency)+" \n" +
				"     } \n" +
				"   } ");
	}
	
	/** Flush and close the file */
	public void destroy() {
		writer.println("] }");
				
		try {
			writer.flush();
			fos.flush();
			fos.close();
		} catch (IOException e) {
			System.err.println("WARNING: Error closing file '"+fileName+"'.");
		}
	}	
	
	@Override
	/** Draw the polarisation information for each ray */
	public void setDrawPolarisationFrames(boolean enable) { this.drawPolarisationFrames = enable; }
	
	/** Draw normal vectors are each ray-surface intersection */
	public void setDrawIntersectionNormals(boolean enable) { this.drawIntersectionNormals = enable; }
	
	/** Draw only the single path of the strongest ray segments for each ray tree */
	public void setDrawOnlyStrongest(boolean enable) { this.drawOnlyStrongest = enable;	}
	
	@Override
	/** Process only every n'th ray drawRay() call */
	public void setSkipRays(int nSkipRays){ this.nSkipRays = nSkipRays; }
	
	/** Length of 'a small line'. Used to draw the polarisation states, and to draw rays going to infinity. */
	public void setSmallLineLength(double smallLineLength) { this.smallLineLength = smallLineLength; }
	
	/** Draw the given ray (tree) in the given RGB color
	 *  * @param colour[r/g/b]  0.0 <= x <=  1.0  */
	public void drawRay(RaySegment ray, double color[]) {
		this.currentColour = color;
		drawRay(ray);
	}
	
	// Ray VRML string builders
	private StringBuilder rayVerticesStr = new StringBuilder();
	private StringBuilder rayColorsStr = new StringBuilder();
	private StringBuilder rayCoordIndicesStr = new StringBuilder();
	private StringBuilder rayColorIndicesStr = new StringBuilder();
	private int nVertices, nColors, nPols;
	
	@Override
	/** Draw the given ray (tree) in the current color */
	public void drawRay(RaySegment ray) {
		if(nSkipped < nSkipRays){ nSkipped++; return; } else { nSkipped = 0; }
		
		// polarisation ellipses are written as we go, but the ray
		// data is collected
		rayVerticesStr.setLength(0);
		rayColorsStr.setLength(0);
		rayCoordIndicesStr.setLength(0);
		rayColorIndicesStr.setLength(0);
		
		writer.println("DEF ray_"+nRays+" Group { children [ \n"); 
				
		nVertices = 0;
		nColors = 0;
		nPols = 0;
		
		if(drawPolarisationFrames || drawIntersectionNormals)
			writer.println("DEF rayPolsAndNorms_"+nRays+" Group { children [ \n");
			
		if(drawPolarisationFrames)
			drawPolarisationFrame(ray, ray.startPos, ray.E0, ray.startIntensity());
		
		
		rayStartIntensity = ray.startIntensity();
		processRaySegment(ray);
		
		if(drawPolarisationFrames || drawIntersectionNormals)
			writer.println("] } \n");
		
		
		//and now the actual ray
		writer.println("DEF rayLine_"+nRays+" Shape { \n" + 
						" appearance USE rayLine \n" +
						" geometry IndexedLineSet { \n" +
						"  coord  Coordinate { point [ \n" +
						rayVerticesStr + " ] } \n" +
						" color Color { color [ \n" + 
						rayColorsStr + " ] } \n" +  
						" colorPerVertex TRUE\n" + 
						" coordIndex [ " + rayCoordIndicesStr + " ] \n" +
						" colorIndex [ " + rayColorIndicesStr + " ] \n" +
						"} } ] }");
		
		nRays++;
	}
	
	/** Semi-recursive ray path descending */
	private void processRaySegment(RaySegment ray) {
				
		//add the start of this segement
		addRaySegment(ray);
		
		while(ray.endHit != null && !ray.endHit.isEnd()) { 
			//while there are outgoing rays...

			//get them, with the strongest first
			List<RaySegment> nextSegs = ray.endHit.getSortedOutgoingRays();

			//for each new branch... (every outgoing ray excluding the strongest)
			if(!drawOnlyStrongest) {
				for(int i=1; i < nextSegs.size(); i++)
					processRaySegment(nextSegs.get(i));	
			}
			//continue with strongest
			ray = nextSegs.get(0);
			addRaySegment(ray);
			
		}
		

	}
	
	private final double[] rotVec(double v[]){
		if(globalRotMat == null)
			return v;
		
		return new double[]{
				globalRotMat[0][0] * v[0] + globalRotMat[0][1] * v[1] + globalRotMat[0][2] * v[2],
				globalRotMat[1][0] * v[0] + globalRotMat[1][1] * v[1] + globalRotMat[1][2] * v[2],
				globalRotMat[2][0] * v[0] + globalRotMat[2][1] * v[1] + globalRotMat[2][2] * v[2], 
		};
	}
	
	/** Adds the ray segment to the drawing, and draws the polarisation state */
	private final void addRaySegment(RaySegment ray){
		
		if(Double.isNaN(ray.startPos[0]) || Double.isNaN(ray.dir[0]) || 
				(ray.endHit != null && Double.isNaN(ray.endHit.pos[0])) )
			return;
		
		double s[] = rotVec(ray.startPos);
		
		rayVerticesStr.append(fmt.format(s[0]));
		rayVerticesStr.append(" ");
		rayVerticesStr.append(fmt.format(s[1]));
		rayVerticesStr.append(" "); 
		rayVerticesStr.append(fmt.format(s[2]));
		rayVerticesStr.append(",\n");
				
		double end[];
		double endIntensity; 
		if(ray.endHit != null){ //end of the line, was absorbed
			end = ray.endHit.pos;
			endIntensity = ray.endIntensity();
		}else{
			end = new double[]{ 
					ray.startPos[0] + smallLineLength * ray.dir[0],
					ray.startPos[1] + smallLineLength * ray.dir[1],
					ray.startPos[2] + smallLineLength * ray.dir[2],
				};
			endIntensity = ray.startIntensity();
		}

		end = rotVec(end);
		rayVerticesStr.append(fmt.format(end[0]));
		rayVerticesStr.append(" ");
		rayVerticesStr.append(fmt.format(end[1]));
		rayVerticesStr.append(" "); 
		rayVerticesStr.append(fmt.format(end[2]));
		rayVerticesStr.append(",\n");
		double cf = ray.startIntensity() / rayStartIntensity;
		if(cf > 1) cf = 1;
		if(cf < 0) cf = 0;
		if(Double.isNaN(cf) || Double.isNaN(currentColour[0])){
			rayColorsStr.append("0 0 0,\n");			
		}else{			
			rayColorsStr.append(fmt.format(currentColour[0] * cf));
			rayColorsStr.append(" ");
			rayColorsStr.append(fmt.format(currentColour[1] * cf));
			rayColorsStr.append(" ");
			rayColorsStr.append(fmt.format(currentColour[2] * cf));
			rayColorsStr.append(",\n");
		}
		
		cf = endIntensity / rayStartIntensity;
		if(cf > 1) cf = 1;
		if(cf < 0) cf = 0;
		if(Double.isNaN(cf) || Double.isNaN(currentColour[0])){
			rayColorsStr.append("0 0 0,\n");			
		}else{			
			rayColorsStr.append(fmt.format(currentColour[0] * cf));
			rayColorsStr.append(" ");
			rayColorsStr.append(fmt.format(currentColour[1] * cf));
			rayColorsStr.append(" ");
			rayColorsStr.append(fmt.format(currentColour[2] * cf));
			rayColorsStr.append(",\n");
		}
		
				
		rayCoordIndicesStr.append(nVertices+ ","+(nVertices+1)+",-1, ");
		nVertices+=2;
		rayColorIndicesStr.append(nColors+ ","+(nColors+1)+",-1, ");
		nColors+=2;
		

		if(drawIntersectionNormals && ray.endHit != null){
			drawPolygonEdge(new double[][]{
					{ end[0] - 0.5*smallLineLength * ray.endHit.normal[0], end[0] + 1.5*smallLineLength * ray.endHit.normal[0] },
					{ end[1] - 0.5*smallLineLength * ray.endHit.normal[1], end[1] + 1.5*smallLineLength * ray.endHit.normal[1] },
					{ end[2] - 0.5*smallLineLength * ray.endHit.normal[2], end[2] + 1.5*smallLineLength * ray.endHit.normal[2] },
			});
			
		}
		
		//draw starting pol at 1/4 distance, and ending pol at 3/4 distance so
		//that they arn't
		//also one at the end of segments that die on an absorper
		boolean absorberHit = ray.endHit != null && ray.endHit.isEnd() && 
				(ray.endHit.surface == null || ray.endHit.surface.getInterface() instanceof Absorber);
		
		if(drawPolarisationFrames) {
			
			
			drawPolarisationFrame(ray, new double[]{
					ray.startPos[0] + 1*(end[0] - ray.startPos[0]) / 6,
					ray.startPos[1] + 1*(end[1] - ray.startPos[1]) / 6,
					ray.startPos[2] + 1*(end[2] - ray.startPos[2]) / 6 }, 
						ray.E0, ray.startIntensity());
			
			//we don't need this if we're going to draw one at the ending surface
			//and also don't need it if the ray never ends
			if(ray.endHit != null && !absorberHit) { 
				drawPolarisationFrame(ray, new double[]{
						ray.startPos[0] + 5*(end[0] - ray.startPos[0]) / 6,
						ray.startPos[1] + 5*(end[1] - ray.startPos[1]) / 6,
						ray.startPos[2] + 5*(end[2] - ray.startPos[2]) / 6 }, 
							ray.E1, ray.endIntensity());
			}
		}
		
		if(absorberHit && drawPolarisationFrames)
			drawPolarisationFrame(ray, ray.endHit.pos, ray.E1, ray.endIntensity());
		
	}
	
	/** Draws the given polarisation state at specified position using the given ray's sense of 'up' */
	public void drawPolarisationFrame(RaySegment ray, double drawPos[], double E[][], double intensity){
		if(Double.isNaN(intensity))
			return;
		
		int nAngs = 16;
		double dAng = 2 * Math.PI / (nAngs-1);
		
		double right[] = Util.cross(ray.dir, ray.up);
		
		for(int j=0; j < E.length; j++){
			if(Double.isNaN(E[j][Pol.uIm]))
				continue;
			
			String name = "pol-n_"+nPols+"_"+j +"-psi_"+polFmt.format(Pol.psi(E[j])*180/Math.PI) +"-chi_"+polFmt.format(Pol.chi(E[j])*180/Math.PI);
			writer.println("DEF "+name+" Shape { \n" + 
					" appearance USE polEllipPolarised"+(j % 4)+" \n" +
					" geometry IndexedLineSet { \n" +
					"  coord  Coordinate { point [");
			
			for(int i=0; i < nAngs; i++){
				double phi = i * dAng;
				double EuOsc = E[j][Pol.uRe]*FastMath.cos(phi) - E[j][Pol.uIm]*FastMath.sin(phi);
				double ErOsc = E[j][Pol.rRe]*FastMath.cos(phi) - E[j][Pol.rIm]*FastMath.sin(phi);
	
				writer.println(
						fmt.format(drawPos[0] + smallLineLength * (EuOsc*ray.up[0] + ErOsc*right[0])) + " " +
						fmt.format(drawPos[1] + smallLineLength * (EuOsc*ray.up[1] + ErOsc*right[1])) + " " +
						fmt.format(drawPos[2] + smallLineLength * (EuOsc*ray.up[2] + ErOsc*right[2])) + ",");
			}
			
			writer.print("] }\n coordIndex [");
			for(int i=0; i < nAngs; i++)
				writer.print(i + ",");
			writer.println("-1 ] } }");
		}
		
		// Vertical
		if(true){
			if(Double.isNaN(drawPos[0]))
				System.out.println("Seriously, wtf?");
			
			writer.println("DEF vert Shape { \n" + 
					" appearance USE rayLine \n" +
					" geometry IndexedLineSet { \n" +
					"  coord  Coordinate { point [");
			
			writer.println(
					fmt.format(drawPos[0]) + " " +
					fmt.format(drawPos[1]) + " " +
					fmt.format(drawPos[2]) + ",");
			

			writer.println(
					fmt.format(drawPos[0] + smallLineLength * 0) + " " +
					fmt.format(drawPos[1] + smallLineLength * 0) + " " +
					fmt.format(drawPos[2] + smallLineLength * 1) + ",");
			
			writer.print("] }\n coordIndex [");
			for(int i=0; i < 2; i++)
				writer.print(i + ",");
			writer.println("-1 ] } }");
		}
		
		//The easiest way to do this, is to spin the phase and draw the real parts of Eu and Er
		
		//		rayCentre, ray.up, right, smallLineLength*ray.intensity(), smallLineLength*ray.stokes[0], "");
		
		String name = "unpol-n_" + nPols + 
					"-I_" + polFmt.format(intensity);

		writer.println("DEF "+name+" Shape { \n" + 
				" appearance USE polEllipIntensity \n" +
				" geometry IndexedLineSet { \n" +
				"  coord  Coordinate { point [");
		
		double magE = FastMath.sqrt(intensity);
		for(int i=0; i < nAngs; i++){
			double phi = i * dAng;
			double EuOsc = magE * FastMath.cos(phi);
			double ErOsc = magE * FastMath.sin(phi);
			
			writer.println(
					fmt.format(drawPos[0] + smallLineLength * (EuOsc*ray.up[0] + ErOsc*right[0])) + " " +
					fmt.format(drawPos[1] + smallLineLength * (EuOsc*ray.up[1] + ErOsc*right[1])) + " " +
					fmt.format(drawPos[2] + smallLineLength * (EuOsc*ray.up[2] + ErOsc*right[2])) + ",");
			
		}
		
		writer.print("] }\n coordIndex [");
		
		for(int i=0; i < nAngs; i++)
			writer.print(i + ",");
			
		writer.println("-1 ] } }");
		
	  
		nPols++;
	}
		
	public void drawOptic(Optic optic) {		
		drawOptic(optic, null, null);
	}
	
	@Override
	/** Adds the given optical element to the VRML */
	public void drawElement(Element elem) {
		if(elem instanceof Optic)
			drawOptic((Optic)elem);
		else
			drawSurface((Surface)elem);
	}
	
	/** Make the given string VRML compatible */
	public String cleanString(String str){
		return str.replaceAll("[#\\.\\t\\s,/\\-\\\\']+", "_");
	}
	
	/** Adds the given optic to the VRML */
	/** Adds the given optic to the VRML, drawn in the specified material */
	public void drawOptic(Optic optic, String forceSurfaceType, double[] forceLineColor) {
		
		writer.println("DEF optic_"+cleanString(optic.getName())+" Group { children [ ");
		
		try{
			for(Surface surface : optic.getSurfaces()) {
				drawSurface(surface, forceSurfaceType, forceLineColor);
			}
			
			//draw sub optics
			for(Optic subOptic : optic.getSubOptics()) {
				drawOptic(subOptic, forceSurfaceType, forceLineColor);
			}
			
			writer.println("]}");
		}catch(RuntimeException err){
			System.err.println("RuntimeException while drawing " + optic.getName());
			throw(err);
		}
	}
	
	/** Adds the given optical surface to the VRML in the given colour */
	public void drawSurface(Surface surface){
		drawSurface(surface, null, null);		
	}
	
	public void drawSurface(Surface surface, String forceSurfaceType){
		drawSurface(surface, forceSurfaceType, null);
	}
	
	/** Adds the given optical surface to the VRML */
	public void drawSurface(Surface surface, String forceSurfaceType, double[] forceLineColor){
		try{
			writer.println("DEF surface_"+cleanString(surface.getName())+" Group { children [ ");
			
			//draw surfaces of this optic
			for(double line[][] : surface.draw()){
				drawPolygonEdge(line, forceLineColor);
				
				int n = line[0].length;
				
				//if the drawing line comes back to (almost) the same point,
				//draw a surface
				double dist2 = FastMath.sqrt(
						FastMath.pow2(line[0][0] - line[0][n-1]) +
						FastMath.pow2(line[1][0] - line[1][n-1]) +
						FastMath.pow2(line[2][0] - line[2][n-1]));
				
				if(dist2 < (smallLineLength/1000)) {
					//try to figure out if its a total mirror/absorber 
					if(surface.getInterface() instanceof Absorber)
						drawClosedPolygon(line, (forceSurfaceType != null) ? forceSurfaceType : "absorbingSurface");
					else if(surface.getInterface() instanceof Reflector)
						drawClosedPolygon(line, (forceSurfaceType != null) ? forceSurfaceType : "reflectiveSurface");
					else
						drawClosedPolygon(line, (forceSurfaceType != null) ? forceSurfaceType : "generalSurface");
					
				}
			}	
			writer.println("]}");
			
		}catch(RuntimeException err){
			System.err.println("RuntimeException while drawing " + surface.getName());
			throw(err);
		}
	}

	/** Starts a new VRML group. endGroup() needs to be called at some point for each startGroup()*/
	public void startGroup(String groupName){
		writer.println("DEF "+cleanString(groupName)+" Group { children [ ");
	}
	
	/** Ends a VRML group */
	public void endGroup(){
		writer.println("]}");
	}
	
	public void drawPolygonEdge(double[][] line) {
		drawPolygonEdge(line, null);
	}
	
	/** Adds an arbitary polygon/line to the VRML 
	 * @param line[x/y/z][vertexIdx] */
	public void drawPolygonEdge(double[][] line, double[] forceLineColor) {
		writer.println(
				"Shape { appearance USE rayLine \n" + 
				"geometry IndexedLineSet { coord  Coordinate { point [");
		
		int nMoves = 0;
		for(int i=0; i < line[0].length; i++) {
			if(Double.isNaN(line[0][i]+line[1][i]+line[2][i]))
				throw new RuntimeException("NaN in VRML");
			
			if(i > 0 && (line[0][i] == line[0][i-1]) && 
					(line[1][i] == line[1][i-1]) &&
					(line[2][i] == line[2][i-1])){
				//point didn't move
				//System.err.println("foad");
			}else{
	
				double l[] = rotVec(new double[]{ line[0][i], line[1][i], line[2][i] });
				writer.println(
						fmt.format(l[0]) + " " +
						fmt.format(l[1]) + " " +
						fmt.format(l[2]) + ",");
				nMoves++;
			}
		}
		if(nMoves < 2)
			throw new RuntimeException("Attempting to draw a polygon with less than 2 sides");
				
		writer.print("] }\n" + 
				"color Color { color [ ");
		if(forceLineColor != null){
			writer.print(fmt.format(forceLineColor[0]) + " " + 
						fmt.format(forceLineColor[1]) + " " + 
						fmt.format(forceLineColor[2]));
		}else{
			writer.print("0 1 0");
		}
		writer.print(" ] } \n" + 
				"colorPerVertex FALSE\n" +
				"coordIndex [ ");
		
		for(int i=0; i < nMoves; i++){
			writer.print(i + ",");
		}		
		writer.println("-1 ] } }");
	}

	/** Adds a filled closed polygon to the VRML
	 * @param line[x/y/z][vertexIdx] */
	public void drawClosedPolygon(double[][] line) {
		drawClosedPolygon(line, "generalSurface");
	}
	
	/** Adds a filled closed polygon to the VRML
	 * @param line[x/y/z][vertexIdx] 
	 * @param surfaceType	Surface style/type to use */
	public void drawClosedPolygon(double[][] line, String surfaceType) {
		
		writer.println(
				"Shape { appearance USE "+surfaceType+" \n" + 
				"geometry IndexedFaceSet { coord  Coordinate { point [");
		
		int nMoves = 0;
		for(int i=0; i < line[0].length; i++) {
			if(Double.isNaN(line[0][i]+line[1][i]+line[2][i]))
				throw new RuntimeException("NaN in VRML");
			
			if(i > 0 && (line[0][i] == line[0][i-1]) && 
						(line[1][i] == line[1][i-1]) &&
						(line[2][i] == line[2][i-1])){
				//didn't move
				//System.err.println("foad");
			}else{
				double l[] = rotVec(new double[]{ line[0][i], line[1][i], line[2][i] });
				
				writer.println(
						fmt.format(l[0]) + " " +
						fmt.format(l[1]) + " " +
						fmt.format(l[2]) + ",");
				nMoves++;
			}
		}
		if(nMoves < 2)
			throw new RuntimeException("Attempting to draw a polygon with less than 2 sides");
		
				
		writer.print("] }\n" + 
					   "coordIndex [ ");
		
		for(int i=0; i < nMoves; i++)
			writer.print(i + ",");				
		writer.println("-1, ");

		for(int i=(nMoves-1); i >= 0; i--)
			writer.print(i + ",");		
		writer.println("-1, ");

		writer.println("] } }");
	}

	/** Add some arbitrary VRML code to the file */
	public void addVRML(String vrmlString) {
		writer.println(vrmlString);
	}

	/** Draw the specified ray and optic into a VRML file with the given name */
	public static void dumpRay(String fileName, Optic optic, RaySegment ray) {
		dumpRay(fileName, optic, ray, 0.001, null, null);
	}
	
	public static void dumpRay(String fileName, Optic optic, RaySegment ray, double smallLineLength) {
		dumpRay(fileName, optic, ray, smallLineLength, null, null);
	}
	
	/** Draw the specified ray and optic into a VRML file with the given name 
	 * @param smallLineLength Length of small line used to draw polarisation states */
	public static void dumpRay(String fileName, Optic optic, RaySegment ray, double smallLineLength, String pre, String post) {
		VRMLDrawer vrmlOut = new VRMLDrawer(fileName, smallLineLength);
		if(pre != null)
			vrmlOut.addVRML(pre);
		vrmlOut.setDrawPolarisationFrames(true);
		vrmlOut.setDrawIntersectionNormals(true);
		if(optic != null)
			vrmlOut.drawOptic(optic);
		if(ray != null)
			vrmlOut.drawRay(ray);
		if(post != null)
			vrmlOut.addVRML(post);
		vrmlOut.destroy();
		
	}


	public void setRayColour(double[] colour) {
		this.currentColour = colour;
	}

	/** Sets a global rotation matrix to be applied to all drawing */ 
	public void setTransformationMatrix(double[][] rotationMatrix) {
		this.globalRotMat = rotationMatrix;
	}

	public String getFileName() { return this.fileName; }	
}
