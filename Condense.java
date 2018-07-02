package dspr_359_1pt688;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;

/**
 * 
 * @author Isaac
 *
 * foo5pt2 takes as input the gff file from repeat masker
 * it then creates a csv with the start and end of each loci, and count/copy number
 * the rows are meaningless across columns as they are UNALIGNED
 * the coordinates of one column are not a one-to-one mapping to coordinates of another column, each column is an individual, different individuals have different coordinates for the "same" loci 
 * the integer argument is how far away two repeats may be and still be considered the 'same' "Array"
 *
 */

public class Condense {

	public static ArrayList<ArrayList<Integer>> fociCountAL; 
	public static ArrayList<ArrayList<Integer>> fociStartAL; 
	public static ArrayList<ArrayList<Integer>> fociEndAL; 
	public static ArrayList<ArrayList<Integer>> alignStartAL; 
	public static ArrayList<ArrayList<Integer>> alignEndAL;
	public static ArrayList<ArrayList<Integer>> alignCountAL;

	public static int control = 0;


	public static void main(String[] args) throws IOException{
		args = new String[]{"C:\\Users\\Isaac\\Documents\\College\\Junior\\Lab_Junior\\bluehive\\rmgff\\all_merged.gff", 
				"C:\\Users\\Isaac\\Documents\\College\\Junior\\Lab_Junior\\bluehive\\rmgff\\OUTPUT_V2.2.csv",
		"100"};

		args = new String[]{"C:\\Users\\Isaac\\Documents\\College\\Junior\\Lab_Junior\\bluehive\\rmgff\\all_merged.gff","400"};

		args = new String[]{"C:\\Users\\Isaac\\Documents\\College\\Junior\\Lab_Junior\\prop\\R_data\\fociCountA1R.csv", 
				"C:\\Users\\Isaac\\Documents\\College\\Junior\\Lab_Junior\\prop\\R_data\\fociStartA1R.csv", 
				"C:\\Users\\Isaac\\Documents\\College\\Junior\\Lab_Junior\\prop\\R_data\\fociEndA1R.csv", 
				"C:\\Users\\Isaac\\Documents\\College\\Junior\\Lab_Junior\\prop\\R_data\\alignStartA1R.csv",
				"C:\\Users\\Isaac\\Documents\\College\\Junior\\Lab_Junior\\prop\\R_data\\alignEndA1R.csv",
		"C:\\Users\\Isaac\\Documents\\College\\Junior\\Lab_Junior\\prop\\R_data\\alignCountCSVR.csv"};

		
		args = new String[] {"D:\\College\\lab\\DSPR_ALL_MERGED_WITH_REF.gff"};
		
		if(args.length == 5){
			foo5pt2(args[0], args[1], args[2], args[3], args[4]);
		} else if( args.length == 2) {
			foo5pt2(args[0], "FOCI_COUNT_" + args[1] + ".csv", "FOCI_START_" + args[1] + ".csv", "FOCI_END_" + args[1] + ".csv", args[1]);
		} else if( args.length == 1){
			foo5pt2(args[0], "D:\\College\\lab\\run_1\\all_loci_count.csv", "D:\\College\\lab\\run_1\\all_loci_start.csv", "D:\\College\\lab\\run_1\\all_loci_end.csv", "800");
		} else if (args.length == 6){
			foo6pt10(args[0], args[1], args[2], args[3], args[4], args[5]);
		}

	}

	public static void foo6pt10(String fociCountCSV, String fociStartCSV, String fociEndCSV, String alignStartCSV, String alignEndCSV, String alignCountCSV) throws IOException {
		int[][] fociCount = parseFile6pt10(244, fociCountCSV);
		int[][] fociStart = parseFile6pt10(244, fociStartCSV);
		int[][] fociEnd = parseFile6pt10(244, fociEndCSV);
		int[][] alignStart = parseFile6pt10(234, alignStartCSV);
		int[][] alignEnd = parseFile6pt10(234, alignEndCSV);
		
		int[][] alignCount = new int[alignStart.length][alignStart[0].length];
		
		for(int i = 0; i < alignCount.length; i++){
			for(int k = 0; k < alignCount[0].length; k++){
				alignCount[i][k] = -1;
			}
		}
		
		for(int i = 0; i < alignCount.length; i++){
			for(int k = 0; k < alignCount[0].length; k++){
				int a = fociStart[i][k];
				int b = fociEnd[i][k];
				
				for(int r = 0; r < alignCount.length; r++){
					int x = alignStart[r][k];
					int y = alignEnd[r][k];
					
					if((a != 0 || b != 0 || x != 0 || y != 0) && (a<=x && b>=y || a>=x && b<=y || a<=x && b<=y && b>=x || a>=x && b>=y && a<=y)){
						alignCount[i][k] = fociCount[r][k];
					}
				}
				
				if(alignCount[i][k] == -1){
					while(alignCount[i][k] == -1){
						int permissiveScore = 1;
						for(int r = 0; r < alignCount.length; r++){
							int x = alignStart[r][k];
							int y = alignEnd[r][k];
							
							if((b-x + permissiveScore) > 0 || (a-y - permissiveScore) < 0 ){
					              alignCount[i][k] = fociCount[r][k];
					           }
						}
						permissiveScore += 100;
					}
				}
			}
			System.out.println(i + "/" + alignCount.length);
		}
		
		System.out.println("writing to file");
		
		try{
			BufferedWriter output = null;
			File file = new File(alignCountCSV);
			output = new BufferedWriter(new FileWriter(file));

			
			for(int i = 0; i < alignCount.length; i++){
				for(int k = 0; k < alignCount[0].length; k++){
					output.write(alignCount[i][k]);
					if(k < alignCount[0].length - 1){
						output.write(",");	
					}
				}
				output.newLine();
			}

			output.close();

		}catch(Exception e){
			System.out.println("could not create file");
		}

	
	}

	public static int[][] parseFile6pt10(int r, String PATH) throws IOException{
		
		//hard coding will give problems later
		int[][] rval = new int[r][14];

		FileInputStream inputStream = null;
		Scanner sc = null;
		try {

			inputStream = new FileInputStream(PATH);
			sc = new Scanner(inputStream, "UTF-8");

			String line = "";


			int row = 0;
			while (sc.hasNextLine()) {

				try{line = sc.nextLine();}
				catch(Exception e){	/* lol */}

				String[] lineAR = line.split(",");

				boolean hasZero = false;
				
				for(int i = 0; i < lineAR.length; i++){
					if(lineAR[i].equals("0")){
						hasZero = true;
					}
				}
				
				if(lineAR[0].contains("A") || lineAR[0].contains("a") || hasZero){
					//do nothing
				} else {

					for(int i = 0; i < lineAR.length; i++){
						if(lineAR[i].equals("")){
							lineAR[i] = "-1";
						}
						rval[row][i] = Integer.parseInt(lineAR[i]);
					}
					row++;	
				}
				
				hasZero = false;
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

		return rval;
	}

	public static void foo6pt9(String fociCountCSV, String fociStartCSV, String fociEndCSV, String alignStartCSV, String alignEndCSV, String alignCountCSV) throws IOException {
		fociCountAL = parseFileFoo6pt9(fociCountCSV);
		fociStartAL = parseFileFoo6pt9(fociStartCSV); 
		fociEndAL = parseFileFoo6pt9(fociEndCSV);
		alignStartAL = parseFileFoo6pt9(alignStartCSV);
		alignEndAL = parseFileFoo6pt9(alignEndCSV);
		alignCountAL = new ArrayList<ArrayList<Integer>>();


		for(int i = 0; i < fociCountAL.size(); i++){
			alignCountAL.add(new ArrayList<Integer>());
		}


		System.out.println(fociCountAL.size());
		System.out.println(fociStartAL.size());
		System.out.println(fociEndAL.size());
		System.out.println(alignStartAL.size());
		System.out.println(alignEndAL.size());
		System.out.println(alignCountAL.size());

		for(int i = 0; i < alignCountAL.size(); i++){
			System.out.print(fociCountAL.get(i).size() + " ");
		} System.out.println("\n");

		for(int i = 0; i < fociStartAL.get(0).size(); i++){
			for(int k = 0; k < fociStartAL.size(); k++){
				setLociFoo6pt9(k, i);
			}
		}

		System.out.println("nonoverlaps: " + control);

		try{
			BufferedWriter output = null;
			File file = new File(alignCountCSV);
			output = new BufferedWriter(new FileWriter(file));

			for(int i = 0; i < alignCountAL.get(0).size(); i++){
				for(int k = 0; k < alignCountAL.size(); k++){
					output.write(alignCountAL.get(k).get(i));
					if(k < alignCountAL.size() - 1){
						output.write(",");
					}
				}
				output.newLine();
			}

			output.close();

		}catch(Exception e){
			System.out.println("could not create file");
		}



	}

	public static void setLociFoo6pt9(int DSPRindividual, int b4LociIndex){

		boolean CONTROL = false;

		for(int i = 0; i < fociCountAL.get(DSPRindividual).size(); i++){

			for(int k = 0; k < fociCountAL.get(b4LociIndex).size(); k++){

				int a = fociStartAL.get(DSPRindividual).get(i);
				int b = fociEndAL.get(DSPRindividual).get(i);
				int x = alignStartAL.get(DSPRindividual).get(k);
				int y = alignEndAL.get(DSPRindividual).get(k);

				if( a <= x && b >= y
						|| a >= x && b <= y
						|| a <= x && b <= y && b >= x
						|| a >= x && b >= y && a <= y){

					CONTROL = true;

					alignCountAL.get(DSPRindividual).add(fociCountAL.get(DSPRindividual).get(i));

					return;
				}
			}
		}

		//TODO write functionality to grab nearest loci since one wasn't found
		//CONTROL is probably false by here, don't think it had much function
		//grab nearby Loci
		//current implementation is a janky patch until then

		System.out.println("grabbing nearby loci, no overlap found");
		control++;

		alignCountAL.get(DSPRindividual).add(fociCountAL.get(DSPRindividual).get(DSPRindividual) % 260);

	}

	public static ArrayList<ArrayList<Integer>> parseFileFoo6pt9(String PATH) throws IOException{
		ArrayList<ArrayList<Integer>> rAL = new ArrayList<ArrayList<Integer>>();
		for(int i = 0; i < 14; i++){

			ArrayList<Integer> tempALINT = new ArrayList<Integer>();

			FileInputStream inputStream = null;
			Scanner sc = null;
			try {

				inputStream = new FileInputStream(PATH);
				sc = new Scanner(inputStream, "UTF-8");

				String line = "";


				while (sc.hasNextLine()) {

					try{line = sc.nextLine();}
					catch(Exception e){

					}


					if(control == 0){
						line = line.substring(1);
						control++;
					}

					line = line.replaceAll(",,", ",-1,");
					if(line.endsWith(",")){
						line = line.concat(",-1");
					}

					ArrayList<String> tempALSTR = new ArrayList<String>(Arrays.asList(line.split(",")));
					//System.out.println(PATH);
					//System.out.println(tempALSTR.size());
					try{
						tempALINT.add(Integer.parseInt(tempALSTR.get(i)));
					} catch(NumberFormatException e){
						tempALINT.add(new Integer(0));
					}



					//					if(control == 0){
					//						line = line.substring(1);
					//						control++;
					//					}
					//System.out.println(line);

					//									ArrayList<String> tempALSTR = new ArrayList<String>(Arrays.asList(line.split(",")));
					//				//System.out.println(tempALSTR);
					//				ArrayList<Integer> tempALINT = new ArrayList<Integer>();
					//				for(String s : tempALSTR){
					//					try{
					//					tempALINT.add(new Integer(Integer.parseInt(s)));
					//					} catch(NumberFormatException e){
					//						tempALINT.add(new Integer(0));
					//						//System.out.println("encountered \" \"");
					//					}
					//				}
					//				rAL.add(tempALINT);

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

			rAL.add(tempALINT);

		}
		return rAL;
	}

	//takes as input a gff file and outputs a file with the start/end and count of the loci
	public static void foo5pt2( String inputPath, String outputPath, String secondOutputPath, String thirdOutputPath, String maxLen) throws IOException {
		String INPUT_PATH = inputPath;
		String OUTPUT_FILE = outputPath;
		int maxLength = Integer.parseInt(maxLen);

		FileInputStream inputStream = null;
		Scanner sc = null;

		//fill individual list
		ArrayList<String> individualls = new ArrayList<String>();
		try {
			inputStream = new FileInputStream(INPUT_PATH);
			sc = new Scanner(inputStream, "UTF-8");

			String line = "";

			while (sc.hasNextLine()) {

				try{line = sc.nextLine();}
				catch(Exception e){

				}

				String[] ind = line.split("\t");
				if(!individualls.contains(ind[0])){
					individualls.add(ind[0]);
				}

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


		for(String s : individualls){
			System.out.print(s  + "\t");
		}System.out.println("");



		ArrayList<ArrayList<Foci>> allResults = new ArrayList<>();

		int currentCount = 1;
		int previousLocation = 0;
		int startingLocation = 0;

		int currentIndividual = 0;

		ArrayList<Foci> tempFoc = new ArrayList<>();

		try {

			inputStream = new FileInputStream(INPUT_PATH);
			sc = new Scanner(inputStream, "UTF-8");

			String line = "";

			while (sc.hasNextLine()) {
				try{line = sc.nextLine();}
				catch(Exception e){
				}

				//code goes here
				String[] lineAR = line.split("\t");

				if(lineAR[0].equals(individualls.get(currentIndividual))){
					if(Math.abs(Integer.parseInt(lineAR[3]) - previousLocation) < maxLength){
						currentCount++;
						previousLocation = Integer.parseInt(lineAR[3]);
					} else {
						Foci tempFoci = new Foci();
						tempFoci.count = currentCount;
						tempFoci.start = startingLocation;
						tempFoci.end = previousLocation;
						tempFoci.name = lineAR[0];

						tempFoc.add(tempFoci);

						startingLocation = Integer.parseInt(lineAR[3]);
						previousLocation = Integer.parseInt(lineAR[4]);
						currentCount = 1;
					}
				} else {
					allResults.add(tempFoc);
					tempFoc = new ArrayList<>();
					currentIndividual++;
				}


				// note that Scanner suppresses exceptions
				if (sc.ioException() != null) {
					throw sc.ioException();
				}
			}
		}
		finally {
			if (inputStream != null) {
				inputStream.close();
			}
			if (sc != null) {
				sc.close();
			}
		}

		allResults.add(tempFoc);

		int largestSize = 0;
		
		//System.out.println(allResults.size());
		for(ArrayList<Foci> t : allResults){
			System.out.print(t.size() + "\t");
			if(t.size() > largestSize) {
				largestSize = t.size();
			}
		}
		
		System.out.println(allResults.size() + "\n");
		
//		for(int k = 0; k < allResults.size(); k++) {
//			System.out.print(allResults.get(k).get(1).name);
//			if(k < allResults.size() - 1) {
//				System.out.print(",");
//			}
//		}
//		System.out.println();
//		for(int i = 1; i < largestSize; i++) {
//			for(int k = 0; k < allResults.size(); k++) {
//				if(i >= allResults.get(k).size() ) {
//					System.out.print("");
//				}else {
//					System.out.print(allResults.get(k).get(i).count);
//				}
//				if(k < allResults.size() - 1) {
//					System.out.print(",");
//				}
//			}
//			System.out.println();
//		}

		//primary output file
		//writes count
		try{
			BufferedWriter output = null;
			File file = new File(OUTPUT_FILE);
			output = new BufferedWriter(new FileWriter(file));

			for(int k = 0; k < allResults.size(); k++) {
				output.write("" + allResults.get(k).get(1).name);
				if(k < allResults.size() - 1) {
					output.write(",");
				}
			}
			output.write("\n");
			for(int i = 1; i < largestSize; i++) {
				for(int k = 0; k < allResults.size(); k++) {
					if(i >= allResults.get(k).size() ) {
						output.write("");
					}else {
						output.write("" +allResults.get(k).get(i).count);
					}
					if(k < allResults.size() - 1) {
						output.write(",");
					}
				}
				output.write("\n");
			}
			
			
//			output.write("individual,");
//			for(int i = 1; i < allResults.get(0).size(); i++){
//				output.write(allResults.get(0).get(i).start + ",");
//			}
//
//			output.write("\n");
//
//			int ii = 0;
//
//			for(ArrayList<Foci> t : allResults){
//				output.write(t.get(ii).name + ",");
//				boolean ctrl = false;
//				for(Foci f : t){
//					if(ctrl){
//						output.write(f.count + ",");
//					} else {
//						ctrl = true;
//					}
//
//				}
//				ii++;
//				output.write("\n");
//			}


			output.close();

		}catch(Exception e){
			System.out.println("could not create file for");
		}

		//secondary output file
		//takes gff, outputs csv with locations instead of foci count
		try{
			BufferedWriter output = null;
			File file = new File(secondOutputPath);
			output = new BufferedWriter(new FileWriter(file));

			for(int k = 0; k < allResults.size(); k++) {
				output.write("" + allResults.get(k).get(1).name);
				if(k < allResults.size() - 1) {
					output.write(",");
				}
			}
			output.write("\n");
			for(int i = 1; i < largestSize; i++) {
				for(int k = 0; k < allResults.size(); k++) {
					if(i >= allResults.get(k).size() ) {
						output.write("");
					}else {
						output.write("" +allResults.get(k).get(i).start);
					}
					if(k < allResults.size() - 1) {
						output.write(",");
					}
				}
				output.write("\n");
			}


			output.close();

		}catch(Exception e){
			System.out.println("could not create file for");
		}


		//tertiary output file
		//takes gff, outputs csv with locations end instead of foci count
		try{
			BufferedWriter output = null;
			File file = new File(thirdOutputPath);
			output = new BufferedWriter(new FileWriter(file));

			for(int k = 0; k < allResults.size(); k++) {
				output.write("" + allResults.get(k).get(1).name);
				if(k < allResults.size() - 1) {
					output.write(",");
				}
			}
			output.write("\n");
			for(int i = 1; i < largestSize; i++) {
				for(int k = 0; k < allResults.size(); k++) {
					if(i >= allResults.get(k).size() ) {
						output.write("");
					}else {
						output.write("" +allResults.get(k).get(i).end);
					}
					if(k < allResults.size() - 1) {
						output.write(",");
					}
				}
				output.write("\n");
			}


			output.close();

		}catch(Exception e){
			System.out.println("could not create file for");
		}
	}


}
