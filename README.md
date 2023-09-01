# Short-Term-SAUPCA-Toolkit
Code to analyze short-term conjunctions using semi-analytical uncertainty propagation.

To initiate the code, the only scripts that need to be edited are: 'Main\_Code.m' and 'constantsAndInitialState.m'. The states and uncertainties of the two objects in conjunction and some constant parameters need to be defined. The input parameters that are required to run the code and that can be changed without needing code updates include:


            Main\_Code.m  & -----------------------------------------------\\
            runMonteCarlo: & When toggled to 1, this flag allows the Monte Carlo probability of collision calculation process to run, using a chosen number of points.\\
            runGMMSTTMethod: & When toggled to 1, this flag allows the semi-analytical GMM-STT probability of collision calculation process to run, using a chosen number of points.\\
            plotThings: & Flag to plot probability of collision results. \\
            saveResultsToText:  & Save results to an output file. \\
            constantsAndInitialState.m  & ----------------------------------------------- \\
            
            [a e i O w M]:   &  [km and rad] Classical orbital elements (COEs) defining the epoch states of the two objects in conjunction.  \\
			P and P2:        & [km and km/s] Cartesian covariance of both objects in conjunction.   \\
			constants.case\_flag: & Pre-coded sample test cases can be accessed using this variables.                  \\
			constants.points:	 & Number of Monte Carlo points currently being evaluated for each object. \\
			constants.JMAX:	&	List of number of GMM components to split each object into. The only limitation is that this script cannot split the distribution into less than 15 components. \\
            constants.testFig2:	&  When toggled to 1, this test flag plots the GMM mean spread for object 1 at the nominal TCA using 15 components. `runMonteCarlo' and `runGMMSTTMethod' must be toggled to 1 to run this. \\
            constants.plot\_GMM\_ell:&	When toggled to 1, this test flag plots the GMM covariances, in addition to their means. `runMonteCarlo' and `runGMMSTTMethod' must be toggled to 1 to run this. \\
            constants.testFig4:	&	When toggled to 1, this test flag plots the weights against the means of the GMM distribution. `runGMMSTTMethod' must be toggled to 1 to run this. \\
		  constants.par: & Allow parallel runs in MATLAB.\\
            constants.nodes: & Number of nodes allowed for the parallel MATLAB runs.\\
            usePredefinedCases: & When toggled to 1, this flag allows the predefined test cases to be accessed, using constants.case\_flag.\\
            constants.P\_prop: & Define Nominal TCA. \\
            ICState.obj2.COE: & Define the COE for object 2.\\
            r1: & Hard body radius of object 1.\\
            r2: & Hard body radius of object 2.\\
            constants.rho: & Reflectivity of the objects. \\
            constants.Aom:  & Area over mass ratio of objects. 
