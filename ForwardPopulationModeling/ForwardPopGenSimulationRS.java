

/**
 * Uses KS test, rejection sampling, 
 */

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Random;
import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest;

public class ForwardPopGenSimulationRS implements Runnable{

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
	public int pp;					//
	public int pmax;					//how many gammas we want in the posterior
	//public double[] gamarray;	//generate uniform prior distrib. of gamma //deprecated
	public double[] posterior; //array containing gammas in the posterior
	public ArrayList<Double> prior;		//array containing gammas in the prior, might get too big
	public long gam_test;				//this is a counter for how many values of gamma that we have tested

	public double minGam;
	public double maxGam;
	public long maxIters = Long.MAX_VALUE;

	public double[] zTransDummy;
	public int[] dummy_dist;

	public long randSeed = 1;
	public long seed; 
	public Random r = new Random(randSeed);
	public long gam_test_max = Long.MAX_VALUE;
	
	public double originalGam = 0;
	
	public ArrayList<int[]> goodDaughters = new ArrayList<>();
	public ArrayList<int[]> badDaughters = new ArrayList<>();

	/**
	 * @deprecated
	 * use the constructors with arguments instead
	 */
	@Deprecated
	public ForwardPopGenSimulationRS() {
		this.omega = 10000;
		this.tmax = 100;
		this.popsize = 1000;
		this.cnumber = 100;
		this.gam = 0.001;
		this.parents = array(popsize, cnumber);
		this.daughters = array(popsize, cnumber);
		this.minGam = 10e-7;
		this.maxGam = 10e-2;

		this.p = 1;
		this.pp = 1;
		this.pmax = 4;

		this.posterior = new double[pmax];


		this.prior = new ArrayList<>();

		this.gam_test = 0;

		this.zTransDummy = new double[0];
	}

	/**
	 * @deprecated
	 * use the constructors with all arguments instead
	 */
	@Deprecated
	public ForwardPopGenSimulationRS(double[] ztrans, int[] dummy) {
		this.omega = 10000;
		this.tmax = 100;
		this.popsize = 1000;
		this.cnumber = 100;
		this.gam = 0.001;
		this.parents = array(popsize, cnumber);
		this.daughters = array(popsize, cnumber);
		this.p = 1;
		this.pp = 1;
		this.pmax = 20;
		this.minGam = 10e-7;
		this.maxGam = 10e-2;

		this.posterior = new double[pmax];
		this.prior = new ArrayList<>();

		this.gam_test = 0;

		this.zTransDummy = ztrans.clone();
		this.dummy_dist = dummy.clone();
		this.r = new Random(System.currentTimeMillis());
	}

	//for Bayesian sampling, with use for an already made dummy population
	public ForwardPopGenSimulationRS(double[] ztrans, int[] dummy, int omega, int tmax, int popsize, int cnumber, double gamma, int pmax, double minGam, double maxGam, long gam_test_max, long maxiters) {
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
		this.pp = 1;
		this.pmax = pmax;

		this.posterior = new double[pmax];
		this.prior = new ArrayList<>();

		this.zTransDummy = ztrans.clone();
		this.dummy_dist = dummy.clone();
		this.r = new Random(seed);
		this.seed = (long)((double)System.currentTimeMillis()*Math.random() + (long)(Math.random() * Integer.MAX_VALUE));
		this.gam_test_max = gam_test_max;
		this.maxIters = maxiters;
	}

	//for generating an initial dummy population
	public ForwardPopGenSimulationRS(int omega, int tmax, int popsize, int cnumber, double gamma, int pmax, double minGam, double maxGam, long seed) {
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
		this.pp = 1;
		this.pmax = pmax;

		this.posterior = new double[pmax];
		this.prior = new ArrayList<>();

		this.gam_test = 0;
		this.zTransDummy = new double[0];
	}


	//	public static void main(String[] args) {
	//
	//		System.out.println("Isaac's ABC");
	//		System.out.println("hello world");
	//
	//		ForwardPopGenSimulation f = new ForwardPopGenSimulation();
	//		f.makeDummyPop();
	//		f.BayesianSampling();
	//	}

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
			gam = (minGam + (maxGam - minGam) * Math.random());	
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

		System.out.println("dummy_dist mean: " + mean(dummy_dist));
		System.out.println("dummy_dist standard deviation: " + sd(dummy_dist));

		System.out.println("zTransDummy mean: " + mean(zTransDummy));
		System.out.println("zTransDummy standard deviation: " + sd(zTransDummy));


		goodDaughters.add(parents.clone());
		badDaughters.add(parents.clone());

	}

	public void BayesianSampling() {
		while( p + 1 <= pmax) {
			
			nextPopulation(false);
			//nextPopulation(true);

			/**
			 * perform ks test
			 */

			var test_distrib = parents.clone();
			if(variance(this.dummy_dist) == 0 && variance(test_distrib) == 0) {
				posterior[p] = gam;
				p++;
				gam_test++;
				double acceptanceRate = (double)((double)p)/((double)gam_test);
				System.out.println("gam: " + gam + "\tp: " + p + "\tgam_test: " + gam_test + "\tgam_accept: " + acceptanceRate + "\tacceptancec_rate: " + acceptanceRate + "\ttime: " + System.currentTimeMillis() + "\tthread: " + Thread.currentThread().getId());
				break;
			} else {
				var zTransTest = zScale(test_distrib);
				double delta = ksTest(this.zTransDummy, zTransTest);

				//test whether to save value in posterior
				if(delta >= 0.05) {
					goodDaughters.add(parents.clone());
					posterior[p] = gam;
					p++;
					gam_test++;
					double acceptanceRate = (double)((double)p)/((double)gam_test);
					System.out.println("gam: " + gam + "\tp: " + p + "\tgam_test: " + gam_test + "\tgam_accept: " + acceptanceRate + "\tacceptancec_rate: " + acceptanceRate + "\ttime: " + System.currentTimeMillis() + "\tthread: " + Thread.currentThread().getId());
				} else {
					badDaughters.add(parents.clone());
					gam_test++;
					if(gam_test >= this.gam_test_max) {
						System.out.println("max gammas tested: " + this.gam_test_max);
						System.out.println("posterior: " + java.util.Arrays.toString(posterior));
						System.out.println("prior: " + prior.toString());
						return;
					}
					if(gam_test % 1000000 == 0) {
						saveDaughters();
						//should it print out which gammas it has tried? 
						System.out.println("gamma's tested: " + gam_test + " @thread: " + Thread.currentThread().getId() + " @time: " + System.currentTimeMillis());
						System.out.println(this.toString());
						//System.out.println("prior: " + prior.toString());
					}
					prior.add(gam);
					if(prior.size() > 10000000) {
						//for excessive memory usage
						//an object has a min memory of 16bytes
						//a double takes 8 bytes
						//an arraylist of size Integer.MAX_VALUE would take 52GB

						//now, with max size 10,000,000
						//max memory usage is about 240MB
						prior.remove(0);
					}
					pp++;
				}
			}
		}

		saveDaughters();
		
		System.out.println("posterior: " + java.util.Arrays.toString(posterior));
		System.out.println("prior: " + prior.toString());

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
	
	public void saveDaughters() {
		
		String goodDaughtersPath = "daugthers_good_gam_" + this.originalGam + "_popsize_" + this.popsize + "_" + Thread.currentThread().getId() + ".csv";
		String badDaughtersPath = "daughters_badd_gam_" + this.originalGam + "_popsize_" + this.popsize + "_" + Thread.currentThread().getId() + ".csv";
		
		//write bad daughters
		try{
			BufferedWriter output = null;
			File file = new File(badDaughtersPath);
			output = new BufferedWriter(new FileWriter(file));

			for(int i = 0; i < badDaughters.get(0).length; i++) {
				for(int k = 0; k < badDaughters.size(); k++) {
					output.write(Integer.toString(badDaughters.get(k)[i]));
					if(k != badDaughters.size() -1) {
						output.write(",");
					}
				}
				if(i != badDaughters.get(0).length - 1) {
					output.write("\n");
				}
			}

			output.close();

		}catch(Exception e){
			System.out.println("could not create file");
		}

		//write good daughters
		try{
			BufferedWriter output = null;
			File file = new File(goodDaughtersPath);
			output = new BufferedWriter(new FileWriter(file));

			for(int i = 0; i < goodDaughters.get(0).length; i++) {
				for(int k = 0; k < goodDaughters.size(); k++) {
					output.write(goodDaughters.get(k)[i]);
					if(k != goodDaughters.size() -1) {
						output.write(",");
					}
				}
				if(i != goodDaughters.get(0).length - 1) {
					output.write("\n");
				}
			}

			output.close();

		}catch(Exception e){
			System.out.println("could not create file");
		}
		
		

		
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

	public double ksTest(double[] n, double[] h) {
		return new KolmogorovSmirnovTest().kolmogorovSmirnovTest(n, h);
	}


	@Deprecated
	public void BayesianSampling_OLD() {
		while( p + 1 <= pmax) {
			parents = array(popsize, cnumber);
			daughters = array(popsize, cnumber);
			gam = (minGam + (maxGam - minGam) * Math.random());
			for(int t = 0; t < tmax; t++) {
				for(int ind = 0; ind < popsize; ind++) {
					double cond = r.nextDouble();
					if(cond <= gam) {
						int parent1 = sample(parents);
						int parent2 = sample(parents);
						double normcon = constant(parent1, parent2);  
						double[] prob = qijk(normcon, parent1, parent2);
						int[] recomb = seq(1, parent1 + parent2 - 1); 

						if(recomb.length == 0) {
							recomb = new int[] {1};
							prob = new double[] {1};
						}

						daughters[ind] = sample(recomb, prob);

						int minCN = 1; 
						int maxCN = parent1 + parent2 - 2;
						double maxWI = Math.max((wi((1/omega), minCN)), (wi((1/omega), maxCN)));
						double minWI = Math.min((wi((1/omega), minCN)), (wi((1/omega), maxCN)));
						double r = minWI + (maxWI - minWI) * Math.random() - .001;

						long iterCount = 0;
						//doesnt work for some reason?
						//ArrayList<Integer> tempRecomb = new ArrayList<Integer>(Arrays.asList(recomb));
						ArrayList<Integer> tempRecomb = new ArrayList<Integer>();
						for(int i : recomb) {tempRecomb.add(i);}
						ArrayList<Double> tempProb = new ArrayList<Double>();
						for(double d : recomb) {tempProb.add(d);}

						//sample without replacement for faster generation
						while(wi((1/omega), daughters[ind]) < r && tempRecomb.size() > 2 && recomb.length != 1 && iterCount < maxIters) {
							daughters[ind] = sample(tempRecomb, tempProb, false); 
							iterCount++;
						}

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

			/**
			 * perform ks test
			 */

			var test_distrib = parents.clone();
			if(variance(this.dummy_dist) == 0 && variance(test_distrib) == 0) {
				goodDaughters.add(parents.clone());
				posterior[p] = gam;
				p++;
				gam_test++;
				System.out.println("gam: " + gam + "\tp: " + p + "\tgam_test: " + gam_test + "\ttime: " + System.currentTimeMillis() + "\tthread: " + Thread.currentThread().getId());
				break;
			} else {
				var zTransTest = zScale(test_distrib);
				double delta = ksTest(this.zTransDummy, zTransTest);

				//test whether to save value in posterior
				if(delta >= 0.05) {
					posterior[p] = gam;
					p++;
					gam_test++;
					System.out.println("gam: " + gam + "\tp: " + p + "\tgam_test: " + gam_test + "\ttime: " + System.currentTimeMillis() + "\tthread: " + Thread.currentThread().getId());
				} else {
					gam_test++;
					if(gam_test >= this.gam_test_max) {
						System.out.println("max gammas tested: " + this.gam_test_max);
						System.out.println("posterior: " + java.util.Arrays.toString(posterior));
						System.out.println("prior: " + prior.toString());
						return;
					}
					if(gam_test % 1000000 == 0) {
						//should it print out which gammas it has tried? 
						System.out.println("gamma's tested: " + gam_test + " @thread: " + Thread.currentThread().getId() + " @time: " + System.currentTimeMillis());
						System.out.println("prior: " + prior.toString());
					}
					prior.add(gam);
					if(prior.size() > 10000000) {
						//for excessive memory usage
						//an object has a min memory of 16bytes
						//a double takes 8 bytes
						//an arraylist of size Integer.MAX_VALUE would take 52GB

						//now, with max size 10,000,000
						//max memory usage is about 240MB
						prior.remove(0);
					}
					pp++;
				}
			}
		}

		System.out.println("posterior: " + java.util.Arrays.toString(posterior));
		System.out.println("prior: " + prior.toString());

	}

	
	@Deprecated
	public void makeDummyPop_OLD() {

		parents = array(popsize, cnumber);
		daughters = array(popsize, cnumber);
		gam = (minGam + (maxGam - minGam) * Math.random());
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
					if(recomb.length == 0) {
						recomb = new int[] {1};
						prob = new double[] {1};
					}

					daughters[ind] = sample(recomb, prob);

					//create a  bounded random decimal value based on min and max possible fitness values
					int minCN = 1; 
					int maxCN = parent1 + parent2 - 2;
					double maxWI = Math.max((wi((1/omega), minCN)), (wi((1/omega), maxCN)));
					double minWI = Math.min((wi((1/omega), minCN)), (wi((1/omega), maxCN)));
					double r = minWI + (maxWI - minWI) * Math.random() - .001;

					long iterCount = 0;
					//doesnt work for some reason?
					//ArrayList<Integer> tempRecomb = new ArrayList<Integer>(Arrays.asList(recomb));
					ArrayList<Integer> tempRecomb = new ArrayList<Integer>();
					for(int i : recomb) {tempRecomb.add(i);}
					ArrayList<Double> tempProb = new ArrayList<Double>();
					for(double d : recomb) {tempProb.add(d);}

					//sample without replacement, as no need to retry the same copy number if we already know it works
					while(wi((1/omega), daughters[ind]) < r && tempRecomb.size() > 2 && recomb.length != 1 && iterCount < maxIters) {
						daughters[ind] = sample(tempRecomb, tempProb, false); 
						iterCount++;
					}

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

		dummy_dist = parents.clone(); //possible pointer error; later I also overwrite daughters; need to come back and re-write for sure
		zTransDummy = new double[0];
		if(variance(dummy_dist) == 0) {
			System.out.println("The variance of the dummy pop is zero");
		} else {
			zTransDummy = zScale(dummy_dist);
		}

		System.out.println("dummy_dist mean: " + mean(dummy_dist));
		System.out.println("dummy_dist standard deviation: " + sd(dummy_dist));

		System.out.println("zTransDummy mean: " + mean(zTransDummy));
		System.out.println("zTransDummy standard deviation: " + sd(zTransDummy));



	}
}
