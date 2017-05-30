package fusionOptics.lenses;

import fusionOptics.optics.SequentialLensSeries;
import fusionOptics.types.Material;

/** (very roughly) Schneider-Kreuznach 'TV' Xenon 25mm F/0.95
 * 
 * Curvatures and positions taken from a really low res graphic 
 * and materials optimised to match the focal length and plane position.
 * 
 * Well, the optimisation works, but doesn't give a nice representation of the lens.
 * I wouldn't trust this too much. 
 * 
 * [ http://www.taunusreiter.de/Cameras/Biotar.html, 
      last image: http://www.taunusreiter.de/Cameras/SchneiderTVXenon.jpg ]
 *
 * [ http://www.marcocavina.com/articoli_fotografici/Rodenstock_De_Oude_f_0,75/05.gif 
 * 	  Different lenses, but might be handy to get general idea of index ordering ]
 * 
 * [ ??? - Case Diagram. I've lost where this came from, but it's in this dir anyway ]
 * 
 * Also, this isn't quite the same one as I have in the lab.
 * 
 */
public class SchneiderXenon25mmF095 extends SequentialLensSeries {
	
	private static final double n[] = new double[]{
		//copy from result of OptimiseIndices
		
		//Hook & Jeeves starting from hand adjusted, ok but image size (focal length) slightly wrong
		//1.3950073242187613,		1.734796142578112,		1.399438476562524,		1.6814160156250004,		1.64,		1.1832031250000008,		1.6098632812499878,		1.7474853515624793
		
		//full GA from random init, 4.05e-4
		//1.344955246204713,	1.8223064223269834,	1.4196847739684784,	1.85,	1.790740073265983,	1.372088890323628,	1.8347289435761394,	1.5872939433040245
		
		//again, with slightly wider allowed ranges, 2.3e-4
//		1.4685906099747106,	1.211826305277503,	1.2226551493473474,	1.16,	1.2614335378335026,	1.1889654879239566,	1.7449175873520546,	1.3151933530746553
//		1.1853812250649647	,1.4767971184888269	,1.3227281651221179	,1.1048122145086994	,1.2734687443013117	,1.0260244968945027	,1.4993048942935974	,1.426486295976534
		1.6903881755946377	,1.2	,1.6734942653889575	,1.2580556227921182	,1.3860690318919402	,1.2459668639846053	,1.5809045822925143	,1.7
	};
	
	public SchneiderXenon25mmF095() {
		super("SchneiderXenon25mmF095", new double[][] { 			 
				   {	  0, 		-32.6, 		0, 		0,		0 		}, //move lens plane so that focal plane is at +25mm
				   {  41.915682,    2.7850102,    15.75,   n[0],    0. },
				   {  191.89743,    0.8431682,    15.75,   0.0,    0. },
				   {  32.180921,    2.9894147,    15.2,    n[1],    0. },
				   {  79.883803,    0.4088088,    15.2,    0.0,    0. },
				   {  18.345297,    6.0043799,    13.1,    n[2],    0. },
				   {  11.510524,    11.5744  ,    10.0,    0.0,    0. },
				   { -12.417569,    0.830393 ,    10.0,    n[3],    0. },
				   {  53.196251,    8.7510643,    12.4,    n[4],    0. },
				   { -18.409173,    0.740966 ,    12.4,    0.0,    0. },
				   { -191.9102 ,    0.8048424,    13.1,    n[5],    0. },
				   {  39.360626,    4.484122 ,    13.1,    n[6],    0. },
				   { -44.445186,    1.0475727,    13.1,    0.0,    0. },
				   {  34.557122,    2.7850102,    11.8,    n[7],    0. },
				   { -107.48223,    13.54,        11.8,    0.0,    0. }
				   // Focal plane here, at +25mm
				}, 
				0.001 , //mm --> m
				0.024);
	}
	
	public SchneiderXenon25mmF095(double[] centre) {
		this();
		shift(centre);
	}

	public double getCaseRadius() { return 0.024; }
}
