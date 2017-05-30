package fusionOptics.lenses;

import algorithmrepository.Algorithms;
import fusionOptics.Util;
import fusionOptics.optics.SequentialLensSeries;

/** Nikon 50mm f/1.1 from US patent 2828671, via http://www.pierretoscani.com/echo_shortpres.html */
public class Nikon50mmF11 extends SequentialLensSeries {
	
	public Nikon50mmF11() { this(1.0); }
	
	public Nikon50mmF11(double scale) {
		super("Nikon50mmF11", new double[][] { 			 
					{  0, 		-36.6, 		0, 		0,		0 		}, //move lens plane so that focal plane is at +50mm
					{  89.089, 	4.62, 	28, 	1.6073,	59.5 	},
					{ 227.667, 	0.32, 	28, 	0, 		0 		},
					{  49.435, 	4.52, 	25, 	1.6073,	59.5 	},
					{  87.070, 	0.74, 	25, 	0, 		0 		},
					{  28.438, 	10.31,	22,		1.7170,	47.9 	}, 
					{ 463.573, 	2.18, 	22,		1.5927, 35.4 	}, 
					{  17.701, 	13.4, 	14.6,		0, 		0		},
					{ -22.644, 	2.87, 	14.6,		1.6483, 33.8 	}, 
					{  72.133, 	10.9, 	17, 	1.7170, 47.9 	},
					{ -30.884, 	0.32, 	17,		0, 		0		},
					{  61.820, 	5.85, 	16,		1.7170, 47.9 	}, 
					{       0, 	0.32, 	16,		0, 		0		},
					{  75.641, 	1.54, 	17,		1.6259, 35.6 	}, 
					{  51.508, 	4.15, 	17,		1.6385, 55.5 	}, 
					{ 275.348, 	0, 		17,		0, 		0		}					
				}, 
				0.001 * scale, //mm --> m
				0.035);
	}

	public Nikon50mmF11(double[] centre) { this(centre, 1.0, null); }
	public Nikon50mmF11(double[] centre, double scale) { this(centre, 1.0, null); }
	public Nikon50mmF11(double[] centre, double scale, double axis[]) {
		this(scale);
		shift(centre);
		if(axis != null){
			double y[] = Util.createPerp(axis);
			double z[] = Util.reNorm(Util.cross(axis, y));
			rotate(centre, Algorithms.rotationMatrix(new double[][]{ axis, y, z }));
		}
	}

	public double getCaseRadius() { return 0.035; }
}
