package fusionOptics.optics;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import fusionOptics.interfaces.IsoIsoStdFresnel;
import fusionOptics.surfaces.Disc;
import fusionOptics.surfaces.Dish;
import fusionOptics.types.Interface;
import fusionOptics.types.Material;
import fusionOptics.types.Medium;
import fusionOptics.types.Optic;
import fusionOptics.types.Surface;

public class Window extends Optic {

	private double centre[];
	private Disc frontSurface, backSurface;
	private Medium windowMedium;
	
	
	public Window(String name, double centre[], double normal[], double radius, double width, Medium windowMedium, Interface iface) {
		super(name);
		this.centre = centre; 
		
		init(centre,  normal, radius, width, windowMedium, iface);
	}
	
	
	private void init(double centre[], double normal[], double radius, double width, Medium windowMedium, Interface iface) {
		this.windowMedium = windowMedium;
		
		double frontCentre[] = new double[3];
		double backCentre[] = new double[3];
		double frontNormal[] = new double[3];
		
		for(int i=0;i<3;i++){
			frontCentre[i] = centre[i] + width * normal[i];
			backCentre[i] = centre[i] - width * normal[i];
			frontNormal[i] = -normal[i];
		}
		
		frontSurface = new Disc(getName() + "-front", frontCentre, frontNormal, radius, windowMedium, null, iface);
		backSurface = new Disc(getName() + "-back", backCentre, normal, radius, windowMedium, null, iface);
			
		addElement(frontSurface);
		addElement(backSurface);
	}
		
}
