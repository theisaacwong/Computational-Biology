package dspr_359_1pt688;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.awt.*;
/**
 * 
 * @author Isaac
 * 
 * For the input files, run Condense.foo5pt2, then run Align.foo6pt10, need to make sure headers line up.have correct labels. 
 * 
 * takes the sam output from bowtie2, though it might also work for bwa sam output
 * for filtering and sorting a sam file based on relation to 1.688 loci (or any given loci)
 * and for recalculating thr MAPQ score for the adjusted reads
 * 
 * for SAM file specifications, see : https://en.wikipedia.org/wiki/SAM_(file_format)
 * for bowtie2  specifications, see : http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
 * 
 * output: named as appropriatley as possible
 * AS_LOCI_X_vs_Loci_NOT_X looks at teh second best loci of everything
 *
 */
public class FilteredSam {

	public static String samFilePath = "";
	public static String lociStartPath = "";
	public static String lociEndPath = "";

	public static String AS_Loci_X_vs_AS_Loci_NOT_X = "";
	
	public static String lociAS_Dist_Path = "";

	public static ArrayList<Integer> lociStart = new ArrayList<Integer>();
	public static ArrayList<Integer> lociEnd = new ArrayList<Integer>();

	public static HashMap<String, ArrayList<String[]>> map = new HashMap<>();

	public static void main(String[] args) throws IOException {
		//serverMainMethod(args);
		localMainMethod();
	}

	public static void localMainMethod() throws IOException {

		long startTime = System.currentTimeMillis();

		int threshold = 0;

		System.out.println("threshold = " + threshold + "\n");
		System.out.println("step\t\ttotal time (ms)\tunpaired reads\tpaired reads");

		samFilePath = "D:\\College\\lab\\DGRP\\RUN_4\\DGRP_to_dmel_3.5.18_2.sam";
		lociStartPath = "D:\\College\\lab\\run_1\\lociStartAlign_OUT.csv";
		lociEndPath = "D:\\College\\lab\\run_1\\lociEndAlign_OUT.csv";

		AS_Loci_X_vs_AS_Loci_NOT_X = "D:\\College\\lab\\DGRP\\RUN_6\\AS_Loci_X_vs_AS_Loci_NOT_X.csv";
		
		lociAS_Dist_Path = "D:\\College\\lab\\DGRP\\RUN_6\\LOCI_AS_DIST.csv";

		String outputPath = "D:\\College\\lab\\DGRP\\RUN_6\\base";

		parseIndFile(lociStartPath, lociStart);
		parseIndFile(lociEndPath, lociEnd);

		String movingFileBase = ".poster";
		movingFileBase = filterOutXReads(movingFileBase);
		movingFileBase = secondarySort(movingFileBase);
		movingFileBase = addLociLabels(movingFileBase);
		movingFileBase = metric(movingFileBase, threshold);

		long endTime   = System.currentTimeMillis();
		long totalTime = endTime - startTime;
		System.out.println("\ntotal time: " + totalTime + "\n");

		java.awt.Toolkit.getDefaultToolkit().beep();  
		
	}

	public static void serverMainMethod(String[] args) throws IOException {
		long startTime = System.currentTimeMillis();

		


		samFilePath = args[0];
		lociStartPath = args[1];
		lociEndPath = args[2];

		int threshold = Integer.parseInt(args[3]);
		
		lociAS_Dist_Path = args[4];
		AS_Loci_X_vs_AS_Loci_NOT_X = args[5];
		
//		String outputPath = args[6];


		System.out.println("threshold = " + threshold + "\n");
		System.out.println("step\t\ttotal time (ms)\tunpaired reads\tpaired reads");
		
		parseIndFile(lociStartPath, lociStart);
		parseIndFile(lociEndPath, lociEnd);

		String movingFileBase = ".inter";
		movingFileBase = filterOutXReads(movingFileBase);
		movingFileBase = secondarySort(movingFileBase);
		movingFileBase = addLociLabels(movingFileBase);
		movingFileBase = metric(movingFileBase, threshold);

		long endTime   = System.currentTimeMillis();
		long totalTime = endTime - startTime;
		System.out.println("\ntotal time: " + totalTime + "\n");
		
	}



	/**
	 * step one, filter out the reads which don't align to the X chromosome
	 */
	public static String filterOutXReads(String s) {
		System.out.print("filtering \t");
		long startTime = System.currentTimeMillis();

		String rval = samFilePath + s + "_1";

		try{
			BufferedWriter output = null;
			File file = new File(rval);
			output = new BufferedWriter(new FileWriter(file));

			FileInputStream inputStream = null;
			Scanner sc = null;
			try {inputStream = new FileInputStream(samFilePath);	sc = new Scanner(inputStream, "UTF-8");
			String line = "";	line = sc.nextLine();
			while (sc.hasNextLine()) {
				try{line = sc.nextLine();}	catch(Exception e){	/* lol */}

				String[] linee = line.split("\t");
				if(linee[2].equalsIgnoreCase("X") && linee[0].contains("@") == false) {
					output.write(line);
					output.write("\n");
				}

			}
			// note that Scanner suppresses exceptions
			if (sc.ioException() != null) {throw sc.ioException();}
			} finally {
				if (inputStream != null) {inputStream.close();}
				if (sc != null) {sc.close();}
			}			

			output.close();
		}catch(Exception e){System.out.println("could not create file");}

		long endTime   = System.currentTimeMillis();
		long totalTime = endTime - startTime;
		System.out.print("total time: " + totalTime + "\n");

		return rval;

	}

	/**
	 * step 1.5, sorting secondary on position, primary on name
	 */
	public static String secondarySort(String ss) throws IOException{
		System.out.print("2nd sorting \t");
		long startTime = System.currentTimeMillis();

		String rval = ss + "_2";

		int errorCounter = 0;
		int correCounter = 0;

		ArrayList<String[]> temp = new ArrayList<>();
		temp.add(new String[]{"@null", "0", "0", "0", "0"});

		try{
			BufferedWriter output = null;
			File file = new File(rval);
			output = new BufferedWriter(new FileWriter(file));

			FileInputStream inputStream = null;
			Scanner sc = null;
			try {inputStream = new FileInputStream(ss);	sc = new Scanner(inputStream, "UTF-8");
			String line = "";	

			String[] linee;
			try{line = sc.nextLine();}	catch(Exception e){	/* lol */}
			linee = line.split("\t");
			temp.add(Arrays.copyOf(linee, linee.length));

			while (sc.hasNextLine()) {

				if(temp.size() == 1) {
					if(sc.hasNextLine()) {
						try{line = sc.nextLine();}	catch(Exception e){	/* lol */}
						linee = line.split("\t");		
					}
				} else {
					temp.add(Arrays.copyOf(linee, linee.length));
					if(sc.hasNextLine()) {
						try{line = sc.nextLine();}	catch(Exception e){	/* lol */}
						linee = line.split("\t");		
					}
				}

				while(temp.get(temp.size()-1)[0].equalsIgnoreCase(linee[0]) && sc.hasNextLine()) {
					temp.add(Arrays.copyOf(linee, linee.length));
					if(sc.hasNextLine()) {
						try{line = sc.nextLine();}	catch(Exception e){	/* lol */}
						linee = line.split("\t");						
					} 
				}
				if(sc.hasNextLine() == false) {
					temp.add(Arrays.copyOf(linee, linee.length));
				}

				if(temp.size()%2 == 1) {
					errorCounter++;
					temp = new ArrayList<String[]>();
				} else {
					correCounter++;
				}

				Collections.sort(temp ,new Comparator<String[]>() {
					public int compare(String[] strings, String[] otherStrings) {
						return Integer.parseInt(strings[3]) - Integer.parseInt(otherStrings[3]);
					}
				});

				//				for(String[] s : temp) {
				//					for(String t : s) {
				//						output.write(t);
				//						if(! s[s.length-1].equals(t)) {
				//							output.write("\t");
				//						}
				//					}
				//					output.write("\n");
				//				}

				for(int i = 0; i < temp.size(); i++) {
					for(int k = 0; k < temp.get(i).length; k++) {
						output.write(temp.get(i)[k]);
						if(k != temp.get(i).length - 1) {
							output.write("\t");
						}
					}
					output.write("\n");
				}

				temp = new ArrayList<String[]>();



			}
			// note that Scanner suppresses exceptions
			if (sc.ioException() != null) {throw sc.ioException();}
			} finally {
				if (inputStream != null) {inputStream.close();}
				if (sc != null) {sc.close();}
			}			

			output.close();
		}catch(Exception e){System.out.println("could not create file");}

		long endTime   = System.currentTimeMillis();
		long totalTime = endTime - startTime;
		System.out.print("total time: " + totalTime +  "\t" + errorCounter + "\t" + correCounter + "\n");

		return rval;
	}

	/**
	 * step two, labels which loci each read aligns to
	 */
	public static String addLociLabels(String ss) throws IOException{
		System.out.print("labeling loci \t");
		long startTime = System.currentTimeMillis();

		String rval = ss + "_3";

		boolean control = true;

		try{
			BufferedWriter output = null;
			File file = new File(rval);
			output = new BufferedWriter(new FileWriter(file));

			FileInputStream inputStream = null;
			Scanner sc = null;
			try {inputStream = new FileInputStream(ss);	sc = new Scanner(inputStream, "UTF-8");
			String line = "";	line = sc.nextLine();
			while (sc.hasNextLine()) {
				try{line = sc.nextLine();}	catch(Exception e){	/* lol */}

				String[] linee = line.split("\t");

				for(int i = 0; i < lociStart.size(); i++) {

					//if(Math.abs(Integer.parseInt(linee[3]) + linee[9].length() - lociStart.get(i)) <= 100 && Math.abs(Integer.parseInt(linee[3]) - lociEnd.get(i)) <= 100 ){ // linee[3]) + linee[9].length() >= lociStart.get(i) && Integer.parseInt(linee[3]) <= lociEnd.get(i)) {
					if(Integer.parseInt(linee[3]) + linee[9].length() + 400 >= lociStart.get(i) && Integer.parseInt(linee[3]) <= lociEnd.get(i) + 400
							|| Integer.parseInt(linee[3]) + linee[9].length() >= lociStart.get(i) && Integer.parseInt(linee[3]) <= lociEnd.get(i)
							|| Integer.parseInt(linee[3]) >= lociStart.get(i) && Integer.parseInt(linee[3]) <= lociEnd.get(i)) {
						output.write(line + "\t" + (i+1));
						//output.write(linee[0] + "\t" + linee[3] + "\t" + (i+1));
						output.newLine();
						control = false;
						break;
					}
				}
				if(control) {
					//output.write(linee[0] + "\t" + linee[3] + "\t" + (0));
					output.write(line + "\t" + (0));
					output.newLine();
				}

				control = true;


			}
			// note that Scanner suppresses exceptions
			if (sc.ioException() != null) {throw sc.ioException();}
			} finally {
				if (inputStream != null) {inputStream.close();}
				if (sc != null) {sc.close();}
			}			

			output.close();
		}catch(Exception e){System.out.println("could not create file");}

		long endTime   = System.currentTimeMillis();
		long totalTime = endTime - startTime;
		System.out.print("total time: " + totalTime + "\n");

		return rval;
	}


	/**
	 * step three
	 * 
	 * the meat of the algorithm,
	 * takes into account mate pair reads, in which the forward and reverse strands will have different reverse and forward XS values 
	 * 
	 * if Q == 1, looks at second highest AS value, 
	 * 		if all top AS values are for the same loci
	 * 			if highest AS not in locus is lower than XS by --threshhold should threshold be a number or a percent?, 
	 * 				write one of the reads of AS == XS, chosen randomly where multiple exist from top AS reads
	 * 					only need one of the F/R to meet the criteria for both to be written
	 * 					only writes one pair of F/R mate pair reads, so one F and one R, and F is the mate of R
	 * if Q == 255
	 * 		write the read?, doesn't seem to matter?
	 * 
	 * if 1 < Q < 255
	 * 		if highest AS value map to the same locus
	 * 			if the highest AS value not mapping to the above locus is lower than XS by a --threshold
	 * 				write the best read, chosen randomly between ties
	 * 				keep Q value
	 * 
	 * will not rewrite Q values or recalculate yet. 
	 */
	public static String metric(String ss, int threshold) {
		System.out.print("metric-ing \t\n");
		long startTime = System.currentTimeMillis();

		HashMap<Integer, ArrayList<Integer>> loci_AS_Dist = new HashMap<Integer, ArrayList<Integer>>();
		TreeSet<String> readRecord = new TreeSet<String>();
		TreeSet<String> keptReads = new TreeSet<String>();


		//this hashmap is for creating the graph described i lab meeting 3/22/18
		//keys are the loci number (1 through 250)
		//values are the coordinates, stored as tuples in an arrayList
		HashMap<Integer, ArrayList<Tuple>> lociAS_vs_notLociAS_graph = new HashMap<Integer, ArrayList<Tuple>>();

		
		String rval  = ss.concat("_4");

		int errorCounter = 0;
		int correCounter = 0;

		ArrayList<String[]> temp = new ArrayList<>();
		temp.add(new String[]{"@null", "-1", "null", "-1", "-1", "null", "null", "-1", "-1", "null", "null", "null","null","null",});

		try{
			BufferedWriter output = null;
			File file = new File(rval);
			output = new BufferedWriter(new FileWriter(file));

			FileInputStream inputStream = null;
			Scanner sc = null;
			try {inputStream = new FileInputStream(ss);	sc = new Scanner(inputStream, "UTF-8");
			String line = "";	

			String[] linee;
			try{line = sc.nextLine();}	catch(Exception e){	/* lol */}
			linee = line.split("\t");
			temp.add(Arrays.copyOf(linee, linee.length));


			while (sc.hasNextLine()) {
				if(temp.size() == 1) {
					if(sc.hasNextLine()) {
						try{line = sc.nextLine();}	catch(Exception e){	/* lol */}
						linee = line.split("\t");		
					}
				} else {
					temp.add(Arrays.copyOf(linee, linee.length));
					if(sc.hasNextLine()) {
						try{line = sc.nextLine();}	catch(Exception e){	/* lol */}
						linee = line.split("\t");		
					}
				}

				while(temp.get(temp.size()-1)[0].equalsIgnoreCase(linee[0]) && sc.hasNextLine()) {
					temp.add(Arrays.copyOf(linee, linee.length));
					if(sc.hasNextLine()) {
						try{line = sc.nextLine();}	catch(Exception e){	/* lol */}
						linee = line.split("\t");						
					}
				}
				if(sc.hasNextLine() == false) {
					temp.add(Arrays.copyOf(linee, linee.length));
				}

				if(temp.size()%2 == 1) {
					errorCounter++;
					temp = new ArrayList<String[]>();
					temp.add(new String[]{"@null", "-1", "null", "-1", "-1", "null", "null", "-1", "-1", "null", "null", "null","null","null",});
				} else {
					correCounter++;
				}

				/*
				 * Meat goes here, AL<String[]> temp stores all the lines for one loci as a String[]
				 */

				/**
				 * if(Q == ?){
				 * 	sort sub array based on AS value?
				 * 	put those AS values in order to an array?
				 * }
				 */


				/**
				 *   if Q == 1, looks at second highest AS value, 
				 * 		if all top AS values are for the same loci
				 * 			--need to check that it does have an XS score, by looking directly at index which should have it, not by looking at just the size of the line/ String[] array
				 * 			-get XS score
				 * 			--compare XS score to all AS scores
				 * 			--if XS score only matches to the same loci
				 * 			if highest AS not in locus is lower than XS by --threshold should threshold be a number or a percent?, 
				 * 				-- find highest AS not in XS loci
				 * 				write one of the reads of AS == XS, chosen randomly where multiple exist from top AS reads
				 * 					only need one of the F/R to meet the criteria for both to be written
				 * 					only writes one pair of F/R mate pair reads, so one F and one R, and F is the mate of R

				 */
				//Q == 1
				//*****************has not yet been protected against edge cases**********
				if(temp.get(0)[4].equals("1")) {


					boolean control1 = true;

					for(int i = 0; i < temp.size(); i++){
						if(!temp.get(i)[12].contains("XS:i:")) {
							control1 = false;
						}
					} 

					if(control1) {

						//F is even , R is odd
						int XS_F = Integer.parseInt(temp.get(0)[12].substring(5));
						int XS_R = Integer.parseInt(temp.get(1)[12].substring(5));

						int temp_AS = 0;
						int second_highest_AS_F = 0;
						int second_highest_AS_R = 0;

						int bestLoci_F = 0;
						int bestLoci_R = 0;

						int bestLociIndex_F = 0;
						int secondLociIndex_F = 0;

						int bestLociIndex_R = 1;
						int secondLociIndex_R = 1;

						//when choosing the second highest loci, will chose F over R if F and R are from different mate pairs; need to compare F to F and R to R
						for(int i = 0; i < temp.size(); i++) {
							temp_AS = Integer.parseInt(temp.get(i)[11].substring(5));

							//even, F
							if((i%2)==0) {
								//found the best mapping, to a non zero locus
								if(XS_F == temp_AS && Integer.parseInt(temp.get(i)[temp.get(i).length - 1]) != 0) {
									//store best loci, and index of best Loci
									bestLoci_F = Integer.parseInt(temp.get(i)[temp.get(i).length - 1]);
									bestLociIndex_F = i;
									//deal with second best mapping
								} else {
									//see if current locus has a higher AS than previously found
									if(Integer.parseInt(temp.get(i)[11].substring(5)) > second_highest_AS_F && bestLoci_F != Integer.parseInt(temp.get(i)[temp.get(i).length - 1])) {
										secondLociIndex_F = i;
										second_highest_AS_F = Integer.parseInt(temp.get(i)[11].substring(5));
									} else {
										//do nothing
									}
								}

								//odd, R
							} else {
								//found the best mapping, to a nonzero locus
								if(XS_R == temp_AS  && Integer.parseInt(temp.get(i)[temp.get(i).length - 1]) != 0) {
									//store best loci, and index of best Loci
									bestLoci_R = Integer.parseInt(temp.get(i)[temp.get(i).length - 1]);
									bestLociIndex_R = i;
									//deal with second best mapping
								} else {
									//see if current locus has a higher AS than previously found
									if(Integer.parseInt(temp.get(i)[11].substring(5)) > second_highest_AS_R && bestLoci_R != Integer.parseInt(temp.get(i)[temp.get(i).length - 1])) {
										secondLociIndex_R = i;
										second_highest_AS_R = Integer.parseInt(temp.get(i)[11].substring(5));
									} else {
										//do nothing
									}
								}
							}
						}

						//System.out.println(temp.size() + "\t" + bestLociIndex_F + "\t" + bestLociIndex_R);

						if(XS_F - second_highest_AS_F > threshold && Integer.parseInt(temp.get(bestLociIndex_F)[temp.get(bestLociIndex_F).length - 1]) != 0) {

							//write using best F
							for(String toWriteF : temp.get(bestLociIndex_F)) {
								output.write(toWriteF + "\t");
							} output.write("\n");
							for(String toWriteR : temp.get(bestLociIndex_F + 1)) {
								output.write(toWriteR + "\t");
							} output.write("\n");

							keptReads.add(temp.get(bestLociIndex_F)[0]);
							

							//lociAS_vs_notLociAS_graph
							//get current loci - bedtLociIndexF
							//see if an AL (value) exists for this locus (key(
							//add the tuple of <XS_F, second_highest_AS_F> to the arraylist for this locus
							if(lociAS_vs_notLociAS_graph.get(Integer.parseInt(temp.get(bestLociIndex_F)[temp.get(bestLociIndex_F).length - 1])) != null){
								lociAS_vs_notLociAS_graph.get(Integer.parseInt(temp.get(bestLociIndex_F)[temp.get(bestLociIndex_F).length - 1])).add(new Tuple(XS_F, second_highest_AS_F));
							} else {
								lociAS_vs_notLociAS_graph.put(Integer.parseInt(temp.get(bestLociIndex_F)[temp.get(bestLociIndex_F).length - 1]), new ArrayList<Tuple>());
								lociAS_vs_notLociAS_graph.get(Integer.parseInt(temp.get(bestLociIndex_F)[temp.get(bestLociIndex_F).length - 1])).add(new Tuple(XS_F, second_highest_AS_F));
							}

						} else if(XS_R - second_highest_AS_R > threshold && Integer.parseInt(temp.get(bestLociIndex_R)[temp.get(bestLociIndex_R).length - 1]) != 0) {

							//write using R
							for(String toWriteR : temp.get(bestLociIndex_R)) {
								output.write(toWriteR + "\t");
							} output.write("\n");
							for(String toWriteF : temp.get(bestLociIndex_R - 1)) {
								output.write(toWriteF + "\t");
							} output.write("\n");

							keptReads.add(temp.get(bestLociIndex_R)[0]);


							//lociAS_vs_notLociAS_graph
							//get current loci - bedtLociIndexF
							//see if an AL (value) exists for this locus (key(
							//add the tuple of <XS_R, second_highest_AS_R> to the arraylist for this locus
							if(lociAS_vs_notLociAS_graph.get(Integer.parseInt(temp.get(bestLociIndex_R)[temp.get(bestLociIndex_R).length - 1])) != null){
								lociAS_vs_notLociAS_graph.get(Integer.parseInt(temp.get(bestLociIndex_R)[temp.get(bestLociIndex_R).length - 1])).add(new Tuple(XS_R, second_highest_AS_R));
							} else {
								lociAS_vs_notLociAS_graph.put(Integer.parseInt(temp.get(bestLociIndex_R)[temp.get(bestLociIndex_R).length - 1]), new ArrayList<Tuple>());
								lociAS_vs_notLociAS_graph.get(Integer.parseInt(temp.get(bestLociIndex_R)[temp.get(bestLociIndex_R).length - 1])).add(new Tuple(XS_R, second_highest_AS_R));
							}
							
							
						} else {
							//mapping is not unique, do not write anything
						}
					}
					//Q == 255
				} else if(temp.get(0)[4].equals("255") && Integer.parseInt(temp.get(0)[temp.get(0).length - 1]) != 0) {
					//					for(String[] s : temp) {
					//						for(String st : s) {
					//							output.write(st + "\t");
					//						} output.write("\n");
					//					}

					// 1 < Q < 255
				} else if(!temp.get(0)[4].equals("null") && !temp.get(0)[4].equals("-1")) {
					boolean control1 = true;

					for(int i = 0; i < temp.size(); i++){
						if(!temp.get(i)[12].contains("XS:i:")) {
							control1 = false;
						}
					} 

					if(control1) {

						//F is even , R is odd
						int XS_F = Integer.parseInt(temp.get(0)[12].substring(5));
						int XS_R = Integer.parseInt(temp.get(1)[12].substring(5));

						//for recalcutaing the copied code
						int XS_F_REAL = XS_F;
						int XS_R_REAL = XS_R;

						XS_F = Integer.parseInt(temp.get(0)[11].substring(5));
						XS_R = Integer.parseInt(temp.get(1)[11].substring(5));


						for(int i = 0; i < temp.size(); i++) {
							if(i%2 == 0) {
								if(XS_F < Integer.parseInt(temp.get(0)[11].substring(5))) {
									XS_F = Integer.parseInt(temp.get(0)[11].substring(5));
								}
							} else {
								if(XS_R < Integer.parseInt(temp.get(0)[11].substring(5))) {
									XS_R = Integer.parseInt(temp.get(0)[11].substring(5));
								}
							}
						}
						//END re-calculations

						int temp_AS = 0;
						int second_highest_AS_F = 0;
						int second_highest_AS_R = 0;

						int bestLoci_F = 0;
						int bestLoci_R = 0;

						int bestLociIndex_F = 0;
						int secondLociIndex_F = 0;

						int bestLociIndex_R = 1;
						int secondLociIndex_R = 1;

						//when choosing the second highest loci, will chose F over R if F and R are from different mate pairs; need to compare F to F and R to R
						for(int i = 0; i < temp.size(); i++) {
							temp_AS = Integer.parseInt(temp.get(i)[11].substring(5));

							//even, F
							if((i%2)==0) {
								//found the best mapping, to a non zero locus
								if(XS_F == temp_AS && Integer.parseInt(temp.get(i)[temp.get(i).length - 1]) != 0) {
									//store best loci, and index of best Loci
									bestLoci_F = Integer.parseInt(temp.get(i)[temp.get(i).length - 1]);
									bestLociIndex_F = i;
									//deal with second best mapping
								} else {
									//see if current locus has a higher AS than previously found
									if(Integer.parseInt(temp.get(i)[11].substring(5)) > second_highest_AS_F && bestLoci_F != Integer.parseInt(temp.get(i)[temp.get(i).length - 1])) {
										secondLociIndex_F = i;
										second_highest_AS_F = Integer.parseInt(temp.get(i)[11].substring(5));
									} else {
										//do nothing
									}
								}

								//odd, R
							} else {
								//found the best mapping, to a nonzero locus
								if(XS_R == temp_AS  && Integer.parseInt(temp.get(i)[temp.get(i).length - 1]) != 0) {
									//store best loci, and index of best Loci
									bestLoci_R = Integer.parseInt(temp.get(i)[temp.get(i).length - 1]);
									bestLociIndex_R = i;
									//deal with second best mapping
								} else {
									//see if current locus has a higher AS than previously found
									if(Integer.parseInt(temp.get(i)[11].substring(5)) > second_highest_AS_R && bestLoci_R != Integer.parseInt(temp.get(i)[temp.get(i).length - 1])) {
										secondLociIndex_R = i;
										second_highest_AS_R = Integer.parseInt(temp.get(i)[11].substring(5));
									} else {
										//do nothing
									}
								}
							}
						}

						//System.out.println(temp.size() + "\t" + bestLociIndex_F + "\t" + bestLociIndex_R);

						if(XS_F - second_highest_AS_F > threshold && Integer.parseInt(temp.get(bestLociIndex_F)[temp.get(bestLociIndex_F).length - 1]) != 0) {

							//write using best F
							for(String toWriteF : temp.get(bestLociIndex_F)) {
								output.write(toWriteF + "\t");
							} output.write("\n");
							for(String toWriteR : temp.get(bestLociIndex_F + 1)) {
								output.write(toWriteR + "\t");
							} output.write("\n");

							keptReads.add(temp.get(bestLociIndex_F)[0]);
							

							//lociAS_vs_notLociAS_graph
							//get current loci - bedtLociIndexF
							//see if an AL (value) exists for this locus (key(
							//add the tuple of <XS_F, second_highest_AS_F> to the arraylist for this locus
							if(lociAS_vs_notLociAS_graph.get(Integer.parseInt(temp.get(bestLociIndex_F)[temp.get(bestLociIndex_F).length - 1])) != null){
								lociAS_vs_notLociAS_graph.get(Integer.parseInt(temp.get(bestLociIndex_F)[temp.get(bestLociIndex_F).length - 1])).add(new Tuple(XS_F, second_highest_AS_F));
							} else {
								lociAS_vs_notLociAS_graph.put(Integer.parseInt(temp.get(bestLociIndex_F)[temp.get(bestLociIndex_F).length - 1]), new ArrayList<Tuple>());
								lociAS_vs_notLociAS_graph.get(Integer.parseInt(temp.get(bestLociIndex_F)[temp.get(bestLociIndex_F).length - 1])).add(new Tuple(XS_F, second_highest_AS_F));
							}

						} else if(XS_R - second_highest_AS_R > threshold && Integer.parseInt(temp.get(bestLociIndex_R)[temp.get(bestLociIndex_R).length - 1]) != 0) {

							//write using R
							for(String toWriteR : temp.get(bestLociIndex_R)) {
								output.write(toWriteR + "\t");
							} output.write("\n");
							for(String toWriteF : temp.get(bestLociIndex_R - 1)) {
								output.write(toWriteF + "\t");
							} output.write("\n");

							keptReads.add(temp.get(bestLociIndex_R)[0]);

							//lociAS_vs_notLociAS_graph
							//get current loci - bedtLociIndexF
							//see if an AL (value) exists for this locus (key(
							//add the tuple of <XS_R, second_highest_AS_R> to the arraylist for this locus
							if(lociAS_vs_notLociAS_graph.get(Integer.parseInt(temp.get(bestLociIndex_R)[temp.get(bestLociIndex_R).length - 1])) != null){
								lociAS_vs_notLociAS_graph.get(Integer.parseInt(temp.get(bestLociIndex_R)[temp.get(bestLociIndex_R).length - 1])).add(new Tuple(XS_R, second_highest_AS_R));
							} else {
								lociAS_vs_notLociAS_graph.put(Integer.parseInt(temp.get(bestLociIndex_R)[temp.get(bestLociIndex_R).length - 1]), new ArrayList<Tuple>());
								lociAS_vs_notLociAS_graph.get(Integer.parseInt(temp.get(bestLociIndex_R)[temp.get(bestLociIndex_R).length - 1])).add(new Tuple(XS_R, second_highest_AS_R));
							}


						} else {
							//mapping is not unique, do not write anything
						}

						//Q == 255
					}
				}





				//Q == 255
				if(temp.get(0)[4].equals("255") && Integer.parseInt(temp.get(0)[temp.get(0).length - 1]) != 0 && temp.size() == 2) {

					int control2 = 0;

					for(String[] s : temp) {

						if(control2 >= 2) {
							break;
						}

						for(String st : s) {
							output.write(st + "\t");
						} output.write("\n");

						control2++;
					}

					keptReads.add(temp.get(0)[0]);

				}


				/**
				 * code for counting reads kept/not kept and for recording loci AS distribution
				 */
				int temp_AS_m;
				int temp_loci;
				for(String[] s : temp) {
					if(s[11].contains("AS:i:")) {
						temp_AS_m = Integer.parseInt(s[11].substring(5));
						temp_loci = Integer.parseInt(s[s.length-1]);
						if(loci_AS_Dist.get(temp_loci) != null) {
							loci_AS_Dist.get(temp_loci).add(temp_AS_m);
						} else {
							loci_AS_Dist.put(temp_loci, new ArrayList<Integer>());
							loci_AS_Dist.get(temp_loci).add(temp_AS_m);
						}
						readRecord.add(s[0]);
					}
				}


				temp = new ArrayList<String[]>();
			}
			// note that Scanner suppresses exceptions
			if (sc.ioException() != null) {throw sc.ioException();}
			} finally {
				if (inputStream != null) {inputStream.close();}
				if (sc != null) {sc.close();}
			}			

			output.close();
		}catch(Exception e){System.out.println("could not create file"); e.printStackTrace();}

		/**
		 * write the locus AS distribution
		 */
		System.out.print("writing Loci AS dist file\n");
		try{
			BufferedWriter output = null;
			File file = new File(lociAS_Dist_Path);
			output = new BufferedWriter(new FileWriter(file));

			int biggestDistribution = 0;
			ArrayList<Integer> im_bad_at_java = new ArrayList<Integer>();

			/**
			 * I'm bad at java, this arraylist is a sorted list of the keys in the hashmap, which are integers,
			 * so really I'm just sorting a list of integers the most complicated way possible, there's definitely better ways to do this
			 * then I'm iterating through this arraylist and using the index as the key to the map
			 * 
			 * I was originally worried that a loci might be missing from the map so iterating through every number might not work since some might be skipped
			 * 
			 */
			
			//errors for loci = 0 or -1
			for(Integer iter : loci_AS_Dist.keySet()) {
				if(iter != 0 && iter != -1) {
					biggestDistribution = Math.max(loci_AS_Dist.get(iter).size(), biggestDistribution);
					Collections.sort(loci_AS_Dist.get(iter));
					im_bad_at_java.add(iter);
				}
			}

			//sorted keys / loci
			//write loci # , which is the 'name' of the loci
			Collections.sort(im_bad_at_java);
			for(int i = 0; i < im_bad_at_java.size(); i++) {
				output.write("" + im_bad_at_java.get(i));
				if(i != im_bad_at_java.size() -1) {
					output.write(",");
				}
			}
			output.write("\n");
			for(int k = 0; k < biggestDistribution; k++) {
				for(int i = 0; i < im_bad_at_java.size(); i++) {

					if(loci_AS_Dist.get(im_bad_at_java.get(i)).size() > k ) {
						output.write("" + loci_AS_Dist.get(im_bad_at_java.get(i)).get(k));
					}
					if(i != im_bad_at_java.size() - 1) {
						output.write(",");
					}
				}
				output.write("\n");
			}


			output.close();

		}catch(Exception e){
			System.out.println("could not create file");
			e.printStackTrace();
		}
		

		/**
		 * write the locus Loci AS vs !Loci AS distribution distribution
		 */
		System.out.print("writing Loci AS cs !lociAS file file\n");
		try{
			BufferedWriter output = null;
			File file = new File(AS_Loci_X_vs_AS_Loci_NOT_X);
			output = new BufferedWriter(new FileWriter(file));

			output.write("LociNumber,Best_AS,NextBestASDifLoci\n");

			

//			int biggestDistribution = 0;
			ArrayList<Integer> im_bad_at_java = new ArrayList<Integer>();
//			lociAS_vs_notLociAS_graph
			
			//errors for loci = 0 or -1
			for(Integer iter : lociAS_vs_notLociAS_graph.keySet()) {
				if(iter != 0 && iter != -1) {
					im_bad_at_java.add(iter);
				}
			}

			//sorted keys / loci
			Collections.sort(im_bad_at_java);
			
			for(int i = 0; i < im_bad_at_java.size(); i++){
//				output.write(im_bad_at_java.get(i) + ",");
				for(int k = 0; k < lociAS_vs_notLociAS_graph.get(im_bad_at_java.get(i)).size(); k++){
					output.write(im_bad_at_java.get(i) + "," + lociAS_vs_notLociAS_graph.get(im_bad_at_java.get(i)).get(k).toString() + "\n");
				}
//				output.write("\n");
			}
			
			output.write("\n");
			
			output.close();

		}catch(Exception e){
			System.out.println("could not create file");
			e.printStackTrace();
		}
		
		

		long endTime   = System.currentTimeMillis();
		long totalTime = endTime - startTime;
		System.out.print("total time: " + totalTime + "\t" + errorCounter + "\t" + correCounter + "\t\t");

		System.out.print("kept / total reads : " + keptReads.size() + " / " + readRecord.size() + "\n");

		return rval;
	}

	public static void parseSamFile(String s) throws IOException {
		FileInputStream inputStream = null;
		Scanner sc = null;
		try {inputStream = new FileInputStream(s);	sc = new Scanner(inputStream, "UTF-8");
		String line = "";	line = sc.nextLine();
		while (sc.hasNextLine()) {
			try{line = sc.nextLine();}	catch(Exception e){	/* lol */}

			String[] linee = line.split("\t");





		}
		// note that Scanner suppresses exceptions
		if (sc.ioException() != null) {
			throw sc.ioException();
		}
		} finally {
			if (inputStream != null) {
				inputStream.close();
			}
			if (sc != null) {
				sc.close();
			}
		}



	}

	public static void parseIndFile(String s, ArrayList<Integer> l) throws IOException {
		FileInputStream inputStream = null;
		Scanner sc = null;
		try {

			inputStream = new FileInputStream(s);
			sc = new Scanner(inputStream, "UTF-8");

			String line = "";
			line = sc.nextLine();

			while (sc.hasNextLine()) {

				try{line = sc.nextLine();}
				catch(Exception e){	/* lol */}

				l.add(Integer.parseInt(line.split(",")[14]));

			}
			// note that Scanner suppresses exceptions
			if (sc.ioException() != null) {
				throw sc.ioException();
			}
		} finally {
			if (inputStream != null) {
				inputStream.close();
			}
			if (sc != null) {
				sc.close();
			}
		}
	}

}
