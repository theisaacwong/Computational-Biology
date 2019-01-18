



public class LauncherWS {


	public static void main(String[] args) {

		//args parser
		int omega = 10000; //original = 10000
		int tmax = 100;
		int popsize = 5000; //original = 1000
		int cnumber = 100;
		double gamma = 0.001;
		int pmax = 100000000;
		int minGam = -8;
		int maxGam = -2;
		long seed = (long)(Integer.MAX_VALUE * Math.random());
		long gam_test_max = Long.MAX_VALUE;
		long maxIters = Long.MAX_VALUE;
		double dumMean = -1;
		double dumVar = -1;
		double dumSize = -1;
		
		System.out.println("seed: " + seed);
		
		if(args.length == 0) {
			//do nothing
		}else if(args.length == 1 && args[0].equals("--help")) {
			System.out.println("--omega\t\tmax CN");
			System.out.println("--tmax\t\tmax generations");
			System.out.println("--popsize\tpopulation size");
			System.out.println("--cnumber\tinitial copy number");
			System.out.println("--gamma\t\trecombination rate");
			System.out.println("--pmax\t\tnumber of gammas in posterior");
			System.out.println("--minGam\t\tminimum gamma");
			System.out.println("--maxGam\t\tmaximum gamma");
			System.out.println("--seed\t\tseed");
			System.out.println("--maxRuns\t\tmaximum gammas to test");
			System.out.println("--maxIters\t\tmax resampling for unfit individuals");
			System.out.println("--dMean\t\tdummy mean");
			System.out.println("--dVar\t\tdummy variance");
			System.out.println("--dSize\t\tdummy size");
			return;
		} else if(args.length %2 != 0) {System.out.println("illegal arg format, type '-help' for options");return;}
		else {
			for(int i = 0; i < args.length; i+=2) {
				if(args[i].equalsIgnoreCase("--gamma")) {
					gamma = Double.parseDouble(args[i+1]);
				} else if(args[i].equalsIgnoreCase("--omega")) {
					omega = Integer.parseInt(args[i+1]);
				} else if(args[i].equalsIgnoreCase("--tmax")) {
					tmax = Integer.parseInt(args[i+1]);
				} else if(args[i].equalsIgnoreCase("--popsize")) {
					popsize = Integer.parseInt(args[i+1]);
				} else if(args[i].equalsIgnoreCase("--pmax")) {
					pmax = Integer.parseInt(args[i+1]);
				} else if(args[i].equalsIgnoreCase("--minGam")) {
					minGam = Integer.parseInt(args[i+1]);
				} else if(args[i].equalsIgnoreCase("--maxGam")) {
					maxGam = Integer.parseInt(args[i+1]);
				} else if(args[i].equalsIgnoreCase("--seed")) {
					seed = Long.parseLong(args[i+1]);
				} else if(args[i].equalsIgnoreCase("--maxRuns")) {
					gam_test_max = Long.parseLong(args[i+1]);
				} else if(args[i].equalsIgnoreCase("--maxIters")) {
					maxIters = Long.parseLong(args[i+1]);
				} else if(args[i].equalsIgnoreCase("--dMean")) {
					dumMean = Double.parseDouble(args[i+1]);
				} else if(args[i].equalsIgnoreCase("--dVar")) {
					dumVar = Double.parseDouble(args[i+1]);
				} else if(args[i].equalsIgnoreCase("--dSize")) {
					dumSize = Double.parseDouble(args[i+1]);
				}
				else {
					System.out.println("couldnt understand: " + args[i]);
				}
			}
		}

		//need to make gamarray static so all threads access it and remove the right gamma value after sampling

		ForwardPopGenSimulationWS intialFPGS = new ForwardPopGenSimulationWS(omega, tmax, popsize, cnumber, gamma, pmax, minGam, maxGam, seed);
		intialFPGS.makeDummyPop();
		double[] initialZTrans = intialFPGS.zTransDummy.clone();
		int[] initialDummyDist = intialFPGS.dummy_dist.clone();

		dumMean = Math.min(dumMean, intialFPGS.dummyMean) == -1 ? intialFPGS.dummyMean : dumMean;
		dumVar = Math.min(dumVar, intialFPGS.dummyVariance) == -1 ? intialFPGS.dummyVariance : dumVar;
		dumSize = Math.min(dumSize, intialFPGS.dummySize) == -1 ? intialFPGS.dummySize : dumSize;
		
		System.out.println("omega: " + omega);
		System.out.println("tmax: " + tmax);
		System.out.println("popsize: " + popsize);
		System.out.println("cnumber: " + cnumber);
		System.out.println("gamma: " + gamma);
		System.out.println("pmax: " + pmax);
		System.out.println("minGam: " + minGam);
		System.out.println("maxGam: " + maxGam);
		System.out.println("seed: " + seed);
		System.out.println("gam_test_max: " + gam_test_max);
		System.out.println("maxIters: " + maxIters);
		
		
		
		int processors = Runtime.getRuntime().availableProcessors();
		for (int i = 0; i < processors; i++) 
		{ 
			Thread object = new Thread(new ForwardPopGenSimulationWS(
					dumMean, 
					dumVar,
					dumSize,
					omega,
					tmax, 
					popsize, 
					cnumber, 
					gamma, 
					pmax, 
					minGam,
					maxGam, 
					gam_test_max,
					maxIters)); 
			object.start(); 

		}



	}


}
