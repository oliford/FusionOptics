package fusionOptics.optimisationMulti;

import fusionOptics.Util;
import fusionOptics.materials.IsotropicFixedIndexGlass;
import fusionOptics.materials.IsotropicLinearDispersiveGlass;
import fusionOptics.types.Element;
import fusionOptics.types.Material;

public class RefractableMaterial extends Parameter {
	
	/** The element to move */
	private Material material;
	
	public RefractableMaterial(Material material, double min, double max) {
		super(min, max); //parameter is relative to the initial bounding sphere centre, so starts at 0
		this.material = material;
		this.initVal = get();
	}
	
	@Override
	public String toString() { 
		return "RefractableMaterial["+material.getClass().getSimpleName()+"]: init = " + initVal() + ", current="+get()+", min = " + min() + ", max = " + max() + ", scale = " + scale;
	}
	
	@Override
	public double get() { //distance from init pos along the given vector
		if(material instanceof IsotropicFixedIndexGlass){
			return ((IsotropicFixedIndexGlass)material).getRefractiveIndex(0, Double.NaN, Double.NaN); //wavelength and temp shouldn't matter
			
		}else if(material instanceof IsotropicLinearDispersiveGlass){
			return ((IsotropicLinearDispersiveGlass)material).getRefractiveIndexD();
			
		}else
			throw new IllegalArgumentException("Unsupported material type " + material.getClass().getSimpleName());
	}
	
	@Override
	public void set(double val) {
		if(material instanceof IsotropicFixedIndexGlass){
			((IsotropicFixedIndexGlass)material).setRefractiveIndex(val);
			
		}else if(material instanceof IsotropicLinearDispersiveGlass){
			((IsotropicLinearDispersiveGlass)material).setRefractiveIndexD(val);
			
		}else
			throw new IllegalArgumentException("Unsupported material type " + material.getClass().getSimpleName());
	}
	
	
}
