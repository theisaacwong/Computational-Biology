

/**
 * Uses t test, weighted sampling, non-uniform prior sampling
 */

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Random;

public class ForwardPopGenSimulationWS implements Runnable{

	//parameters to generate initial dummy population
	public int omega; 				// maximum
	public int tmax; 				//generation
	public int popsize;				//population size
	public int cnumber; 				// initial copy number
	public double gam;					// recombination rate
	public int[] parents; 	// initial array
	public int[] daughters;	// also an initial array

	//Parameters for Bayesian sampling
	public int p;					//corresponds to position in posterior vector; keeps track of how many gammas we have in posterior
	public int pmax;					//how many gammas we want in the posterior
	//public double[] gamarray;	//generate uniform prior distrib. of gamma //deprecated
	public long gam_test;				//this is a counter for how many values of gamma that we have tested

	public int minGam;
	public int maxGam;
	public long maxIters = Long.MAX_VALUE;

	public double[] zTransDummy;
	public int[] dummy_dist;

	public long randSeed = 1;
	public long seed; 
	public Random r = new Random(randSeed);
	public long gam_test_max = Long.MAX_VALUE;
	
	public double originalGam = 0;
	
	public static String tPath = "t_values.csv";
	public ArrayList<Double> tValues = new ArrayList<>();
	public ArrayList<Double> gValues = new ArrayList<>();
	
	public double dummyMean;
	public double dummyVariance;
	public double dummySize;

	/**
	 * @deprecated
	 * use the constructors with arguments instead
	 */
	@Deprecated
	public ForwardPopGenSimulationWS() {}

	/**
	 * @deprecated
	 * use the constructors with all arguments instead
	 */
	@Deprecated
	public ForwardPopGenSimulationWS(double[] ztrans, int[] dummy) {}

	//for Bayesian sampling, with use for an already made dummy population
	public ForwardPopGenSimulationWS(double[] ztrans, int[] dummy, int omega, int tmax, int popsize, int cnumber, double gamma, int pmax, int minGam, int maxGam, long gam_test_max, long maxiters) {
		this.omega = omega;
		this.tmax = tmax;
		this.popsize = popsize;
		this.cnumber = cnumber;
		this.gam = gamma;
		this.originalGam = gamma;
		this.parents = array(popsize, cnumber);
		this.daughters = array(popsize, cnumber);
		this.minGam = Math.min(minGam, maxGam);
		this.maxGam = Math.max(minGam, maxGam);

		this.p = 1;
		this.pmax = pmax;

		this.zTransDummy = ztrans.clone();
		this.dummy_dist = dummy.clone();
		this.r = new Random(seed);
		this.seed = (long)((double)System.currentTimeMillis()*Math.random() + (long)(Math.random() * Integer.MAX_VALUE));
		this.gam_test_max = gam_test_max;
		this.maxIters = maxiters;
		
		this.minGam = Math.abs(this.minGam);
		this.maxGam = Math.abs(this.maxGam);
		
	}
	
	//for Bayesian sampling, with use for an already made dummy population, for t-testing - don't need original dummy pop, just mean and variance
	public ForwardPopGenSimulationWS(double dumMean, double dumVar, double dumSize, int omega, int tmax, int popsize, int cnumber, double gamma, int pmax, int minGam, int maxGam, long gam_test_max, long maxiters) {
		this.omega = omega;
		this.tmax = tmax;
		this.popsize = popsize;
		this.cnumber = cnumber;
		this.gam = gamma;
		this.originalGam = gamma;
		this.parents = array(popsize, cnumber);
		this.daughters = array(popsize, cnumber);
		this.minGam = Math.min(minGam, maxGam);
		this.maxGam = Math.max(minGam, maxGam);

		this.p = 1;
		this.pmax = pmax;

		this.dummyMean = dumMean;
		this.dummySize = dumSize;
		this.dummyVariance = dumVar;
		
		this.r = new Random(seed);
		this.seed = (long)((double)System.currentTimeMillis()*Math.random() + (long)(Math.random() * Integer.MAX_VALUE));
		this.gam_test_max = gam_test_max;
		this.maxIters = maxiters;
		
		this.minGam = Math.abs(this.minGam);
		this.maxGam = Math.abs(this.maxGam);
	}

	//for generating an initial dummy population
	public ForwardPopGenSimulationWS(int omega, int tmax, int popsize, int cnumber, double gamma, int pmax, int minGam, int maxGam, long seed) {
		this.omega = omega;
		this.tmax = tmax;
		this.popsize = popsize;
		this.cnumber = cnumber;
		this.gam = gamma;
		this.originalGam = gamma;
		this.parents = array(popsize, cnumber);
		this.daughters = array(popsize, cnumber);
		this.minGam = Math.min(minGam, maxGam);
		this.maxGam = Math.max(minGam, maxGam);

		this.seed = seed;
		this.r = new Random(this.seed);

		this.p = 1;
		this.pmax = pmax;

		this.gam_test = 0;
		this.zTransDummy = new double[0];
		
		this.minGam = Math.abs(this.minGam);
		this.maxGam = Math.abs(this.maxGam);
	}


	public void run() {
		try
		{ 
			// Displaying the thread that is running 
			System.out.println ("Thread " + Thread.currentThread().getId() + " is running. " + "Begin Bayesian sampling. @time: " + System.currentTimeMillis() + "\tseed: " + this.seed); 
			this.BayesianSampling();
		} 
		catch (Exception e) 
		{ 
			// Throwing an exception
			System.out.println ("Exception is caught");
			e.printStackTrace();
		} 
	}

	public String toString() {
		
		String s = "Thread: " + Thread.currentThread().getId() + 
		
		"\nomega: " + omega + 
		"\ntmax: " + tmax + 
		"\npopsize: " + popsize + 
		"\ncnumber: " + cnumber + 
		"\ngamma: " + originalGam + 
		"\npmax: " + pmax + 
		"\nminGam: " + minGam + 
		"\nmaxGam: " + maxGam + 
		"\nseed: " + seed + 
		"\ngam_test_max: " + gam_test_max + 
		"\nmaxIters: " + maxIters + 
		"\np: " + p;
		
		return s;
	}
	
	public void nextPopulation(boolean isBayesian) {
		parents = array(popsize, cnumber);
		daughters = array(popsize, cnumber);
		if(isBayesian) {
			double minimum = Math.pow(10, -1*minGam);
			double value = 1 + (10 - 1) * r.nextDouble();
			value = Math.random(); // either seems to be fine, just off by one order of magnitude
			double power = -1*( maxGam + r.nextInt((minGam - maxGam)));
			gam = minimum + value * Math.pow(10, power);
			//gam = (minGam + (maxGam - minGam) * Math.random());	
		}
		for(int t = 0; t < tmax; t++) { //iterate over tmax generations
			for(int ind = 0; ind < popsize; ind++) {  //iterate over all daughters
				double cond = r.nextDouble(); //generate a random value for recombination
				if(cond <= gam) {
					int parent1 = sample(parents);
					int parent2 = sample(parents);
					double normcon = constant(parent1, parent2);  
					double[] prob = qijk(normcon, parent1, parent2);
					int[] recomb = seq(1, parent1 + parent2 - 1); 

					//avoid indexoutofbounds exceptions
					if(recomb.length == 0 || prob.length == 0) {
						recomb = new int[] {1};
						prob = new double[] {1};
					}

					daughters[ind] = sample(recomb, prob);

					//NO FITNESS!!!!
//					//create a  bounded random decimal value based on min and max possible fitness values
//					int minCN = 1; 
//					int maxCN = parent1 + parent2 - 1;
//					double maxWI = Math.max((wi((1/omega), minCN)), (wi((1/omega), maxCN)));
//					double minWI = Math.min((wi((1/omega), minCN)), (wi((1/omega), maxCN)));
//					double r = minWI + (maxWI - minWI) * Math.random() - .001;
//
//					long iterCount = 0;
//					//doesnt work for some reason?
//					//ArrayList<Integer> tempRecomb = new ArrayList<Integer>(Arrays.asList(recomb));
//					ArrayList<Integer> tempRecomb = new ArrayList<Integer>();
//					for(int i : recomb) {tempRecomb.add(i);}
//					ArrayList<Double> tempProb = new ArrayList<Double>();
//					for(double d : recomb) {tempProb.add(d);}
//
//					//sample without replacement, as no need to retry the same copy number if we already know it works
//					while(wi((1/omega), daughters[ind]) < r && tempRecomb.size() > 2 && recomb.length != 1 && iterCount < maxIters) {
//						daughters[ind] = sample(tempRecomb, tempProb, false); 
//						iterCount++;
//					}

				} else {
					int parent1 = sample(parents);
					daughters[ind] = parent1;
					if(daughters[ind] > omega) {
						daughters[ind] = omega;
						System.err.println("daughter had more than omega CN!!!!!");
					}
				}
			}
			parents = daughters.clone();
		}
	}

	/**
	 * makes initial dummy population
	 */
	public void makeDummyPop() {

		nextPopulation(false);	

		dummy_dist = parents.clone(); //possible pointer error; later I also overwrite daughters; need to come back and re-write for sure
		zTransDummy = new double[0];
		if(variance(dummy_dist) == 0) {
			System.out.println("The variance of the dummy pop is zero");
		} else {
			zTransDummy = zScale(dummy_dist);
		}

		this.dummyMean = mean(dummy_dist);
		this.dummyVariance = var(dummy_dist);
		this.dummySize = this.popsize;
		
		System.out.println("dummy_dist mean: " + mean(dummy_dist));
		System.out.println("dummy_dist standard deviation: " + sd(dummy_dist));

		System.out.println("zTransDummy mean: " + mean(zTransDummy));
		System.out.println("zTransDummy standard deviation: " + sd(zTransDummy));


	}

	public void BayesianSampling() {
		
		ArrayList<Integer> percents = new ArrayList<>();
		for(int i = 1; i < 100; i++) {
			percents.add(i);
		}
		
		while( p <= pmax) {
			p++;
			
			nextPopulation(true);

			/**
			 * weighted sampling
			 */
			var test_distrib = parents.clone();
			
			double sampleMean = mean(test_distrib);
			double sampleVar = var(test_distrib);
			
			double t = (dummyMean - sampleMean) / Math.sqrt(dummyVariance/dummySize + sampleVar/popsize);
			
			tValues.add(t);
			gValues.add(gam);
			
			if(p%1000 == 0) {System.out.println(p + " / " + pmax);}
			
			System.out.println(t + "," + gam);
			
		}
		
		//write t values
		try{
			BufferedWriter output = null;
			File file = new File(tPath + seed + ".txt");
			output = new BufferedWriter(new FileWriter(file));
			
			output.write("t,gamma\n");
			for(int i = 0; i < tValues.size(); i++) {
				output.write(tValues.get(i) + "," + gValues.get(i) + "\n");
			}

			output.close();

		}catch(Exception e){
			System.out.println("could not create file");
		}

	}




	/**
	 * #EQUATION 1#
	 * #BELOW IS A FITNESS FCN THAT IS NOT INVOLVED IN MODEL#
	 * #"DECREASING FCN"# 
	 */
	public double wi(double s, double i) {
		return 1 - s * (i - 1);
	}

	/**
	 * #EQUATION 2B#
	 * #THE CONSTANT c IS INVOLVED IN CALCULATION OF THE PROBABILITY OF EXCHANGE#
	 * #OUTPUT DEPENDS ON PARENT SUM = EVEN OR ODD#
	 */
	public double constant(double j, double k) {
		if((j+k)%2==0){
			double even = 2/(j+k);
			return even;
		}else{
			double odd = (2*(j+k))/(Math.pow((j+k), 2)-1);
			return odd;
		}
	}

	/**
	 * #EQUATION 2A#
	 * #CALCULATES PROB OF EXCHANGE#
	 * #CONFUSED BECAUSE THE PAPER SAYS GAMMA*Qijk IS THE PROB OF EXCHANGE, WHERE GAMMA IS THE RECOMB RATE#
	 * #CREATES VECTOR OF DIFF COPY NUMBER VALUES FOR DAUGHTER = THE NUMBER I#
	 *
	 * j + k must be greater than 1
	 * returns an array of size (j+k-1)
	 * 
	 * @param normcon
	 * @param j
	 * @param k
	 * @return
	 */
	public double[] qijk(double normcon, double j, double k) {
		double[] honk = arrayd(Math.max((int)(j+k-1), 0), 0);

		if(honk.length == 0) {
			return new double[] {0.0};
		}

		for(int i = 0; i < honk.length; i++) {
			honk[i] = normcon * (1-Math.abs(((2*(i + 1))/(j+k))-1));
		}

		return honk.clone();
	}
	

	/**
	 * @param len - length of array
	 * @param val - value of each element
	 * @return an int[] array
	 */
	public int[] array(int len, int val) {
		int[] rval = new int[len];
		for(int i = 0 ; i < len; i++) {
			rval[i] = val;
		}
		return rval;
	}

	/**
	 * (inclusive, exclusive)
	 */
	public double[] arrayd(int len, int val) {
		double[] rval = new double[len];
		for(int i = 0 ; i < len; i++) {
			rval[i] = val;
		}
		return rval;
	}

	/**
	 * (inclusive, inclusive), returns of length (end-start + 1)
	 */
	public int[] seq(int start, int end) {
		int[] rval = new int[Math.max(end, start) - Math.min(end, start) + 1];
		for(int i = 0; i < rval.length; i++) {
			rval[i] = start + i;
		}
		return rval;
	}

	public double variance(int[] d) {
		double mean = mean(d);
		double sum = 0;
		for(int i : d) {
			sum += Math.pow(mean - (double)i, 2);
		}
		sum /= d.length;
		return sum;
	}

	public double variance(double[] d) {
		double mean = mean(d);
		double sum = 0;
		for(double i : d) {
			sum += Math.pow(mean - (double)i, 2);
		}
		sum /= d.length;
		return sum;
	}

	public double var(int[] d) {
		double mean = mean(d);
		double sum = 0;
		for(int i : d) {
			sum += Math.pow(mean - (double)i, 2.0);
		}
		sum /= d.length;
		return sum;
	}

	public double mean(int[] d) {
		double sum = 0;
		for(int i : d) {
			sum += (double)i;
		}
		sum /= d.length;
		return sum;
	}

	public double mean(double[] d) {
		double sum = 0;
		for(double i : d) {
			sum += (double)i;
		}
		sum /= d.length;
		return sum;
	}

	public double sd(int[] d) {
		return Math.sqrt(variance(d));
	}

	public double sd(double[] d) {
		return Math.sqrt(variance(d));
	}

	public int sample(int[] d) {
		return d[r.nextInt(d.length)];
	}

	public double sample(double[] d) {
		return d[r.nextInt(d.length)];
	}

	public double[] zScale(int[] d) {
		double[] dd = new double[d.length];

		double mean = mean(d);
		double sd = sd(d);

		for(int i = 0; i < d.length; i++) {
			dd[i] = (d[i] - mean)/sd;
		}

		return dd;
	}

	public int min(int[] a) {
		int min = Integer.MAX_VALUE;
		for(int i : a) {
			if(i < min) {
				min = i;
			}
		}
		return min;
	}

	public int max(int[] a) {
		int min = Integer.MIN_VALUE;
		for(int i : a) {
			if(i > min) {
				min = i;
			}
		}
		return min;
	}

	public int sample(ArrayList<Integer> pop, ArrayList<Double> weights, boolean replace) {
		if(replace == true) {
			double totalWeight = 0.0;
			for(Double d : weights) {
				totalWeight += d;
			}
			double weightedRandom = r.nextDouble() * totalWeight;
			double weightedCount = 0.0;

			for(int i = 0; i < pop.size(); i++) {
				weightedCount += weights.get(i);
				if(weightedCount >= weightedRandom) {
					return pop.get(i);
				}
			}

			System.out.println("Error code 729472");
			return pop.get(pop.size()-1);

		} else {
			double totalWeight = 0.0;
			for(Double d : weights) {
				totalWeight += d;
			}
			double weightedRandom = r.nextDouble() * totalWeight;
			double weightedCount = 0.0;

			for(int i = 0; i < pop.size(); i++) {
				weightedCount += weights.get(i);
				if(weightedCount >= weightedRandom) {

					int tempInd = pop.get(i);

					weights.remove(i);
					pop.remove(i);

					return tempInd;
				}
			}

			System.out.println("error code 38264922");
			return pop.get(pop.size()-1);
		}
	}

	public int sample(int[] population, double[] weights) {
		double totalWeight = 0.0;
		for(double d : weights) {
			totalWeight += d;
		}

		double weightedRandom = r.nextDouble() * totalWeight;
		double weightedCount = 0.0;

		for(int i = 0; i < population.length; i++) {
			weightedCount += weights[i];
			if(weightedCount >= weightedRandom) {
				return population[i];
			}
		}
		System.err.println("Error code: 2397");
		return population[population.length];

	}


}
