package fusionOptics.lenses;

import fusionOptics.interfaces.Absorber;
import fusionOptics.optics.SequentialLensSeries;
import fusionOptics.surfaces.Iris;

/** Wide angle lens (120°) from Zemax data base (F_004.ZMX) f/2, f=30.98562 mm) */
public class WideAngle120deg31mmF2 extends SequentialLensSeries {
	
	public WideAngle120deg31mmF2() { this(1.0); }
	
	public WideAngle120deg31mmF2(double scale) {
		super("WideAngle120deg31mmF2", new double[][] { 			 
		
				/** each entry of data[][] is:
				 * radius.curv,	distance to next, radius, index_d, abbe_d */ 
					{	0			,	-304		,	0		,	0			,	0	},  //move lens plane so that focal plane is at appr. +31mm
					{	0			,	17.5514	,	88.265	,	1.56885954	,	60	},
					{	61.83122	,	22.352	,	54.61	,	0			,	0	},
					{	68.6943		,	32.2072	,	51.816	,	1.71742875	,	60	},
					{	-219.8395	,	18.3134	,	51.816	,	1.61275053	,	60	},
					{	31.95574	,	17.907	,	25.654	,	0			,	0	},
					{	-70.29704	,	17.4498	,	25.654	,	1.71742875	,	60	},
					{	-55.27294	,	17.5514	,	26.924	,	1.61275053	,	60	},
					{	-111.4019	,	5.588	,	26.924	,	0			,	0	},
					{	368.7978	,	14.4272	,	22.86	,	1.61275053	,	60	},
					{	-39.50208	,	6.3246	,	22.86	,	1.71742875	,	60	},
					{	-132.4915	,1.5748+17.5514,22.86	,	0			,	0	}, // Additional thickness of the STOP in Zemax file
					{	-263.7917	,	11.8364	,	28.829	,	1.61275053	,	60	},
					{	-65.7225	,	16.8656	,	28.829	,	0			,	0	},
					{	2414.638	,	8.3058	,	34.29	,	1.61275053	,	60	},
					{	-129.0218	,	7.2136	,	34.29	,	0			,	0	},
					{	-55.9308	,	6.1468	,	33.655	,	1.71742875	,	60	},
					{	-214.6478	,	5.8166	,	38.481	,	0			,	0	},
					{	290.6497	,	8.001	,	47.625	,	1.71742875	,	60	},
					{	117.762		,	1.651	,	45.212	,	0			,	0	},
					{	128.9939	,	25.3746	,	47.625	,	1.61275053	,	60	},
					{	-88.86952	,	1.5494	,	47.625	,	0			,	0	},
					{	81.93024	,	20.5994	,	47.625	,	1.61275053	,	60	},
					{	603.2144	,	0		,	47.625	,	0			,	0	}
		
		}, 
				0.001 * scale, //mm --> m
				0.2 * scale); // Take a bit larger case 
		
		/** Add the iris provided by the Zemax file */
		this.addElement(new Iris("IrisZemax", 
				new double[]{0.001*scale*(-304 + 17.5514 + 22.352 + 32.2072 + 0.0 + 18.3134 + 17.907 + 17.4498 + 0.0 + 17.5514 + 5.588 + 14.4272 + 0.0 + 6.3246 + 1.5748), 0,0},
				new double[]{1,0,0},
				0.001*scale * 200.0,
				0.001*scale * 0.5*36.06057, 
				Absorber.ideal()));
	}

	public WideAngle120deg31mmF2(double[] centre) { this(centre, 1.0); }
	public WideAngle120deg31mmF2(double[] centre, double scale) {
		this(scale);
		shift(centre);
	}

	/***
	 * Returns lens case radius (without scaling factor, i.e. along to original data sheet
	 * @return
	 */
	public double getCaseRadius() { return 0.2; }
}
