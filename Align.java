package dspr_359_1pt688;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Scanner;

/**
 * 
 * @author Isaac
 *
 * after running Condense.foo5pt2, run Align.foo.6pt10
 * this method takes as input the following
 * 	1. the unaligned counts file from Condense.foo5pt2
 * 	2. the unaligned start file from Condense.foo5pt2
 *	3. the unaligned end file from Condense.foo5pt2
 *	4. the "aligned start coordinates" from mauve aligner, the file is from my altered version of mauve aligner - coordinates were aligned to the first column, assumed
 *	5. the "aligned end coordinates" from mauve aligner, the file is from my altered version of mauve aligner
 *
 *	outputs: the aligned counts file which has continuity across rows. one row is the same loci across all columns/indiivudals. 
 *
 *	the files will need to be altered a bit to be read by R's read.csv, including adding/removing the header first line
 *
 *	the algorithm for aligning works as follows
 *		based on the first column(should be the 'reference individual'), get the start and end
 *		for each other column, find a row/index/loci which has overlapping coordinates according to the aligned coordinates
 *			for example if my reference loci has a loci with coordinates [100,200] 
 *			and a different individual has another loci with coordinates [300,400] 
 *			and mauve says that reference [100,200] corresponds to [350,450] of this 
 *			particular individual, this is considered the same loci since we find a 
 *			loci with coordinates overlapping with where we were expecting
 *
 * 
 * 		once an orthologous loci has been found, we take track of all the statistics for this loci(start, end, count)
 * 		and save it to a new file, named accordingly
 */

public class Align {

	public static String LOCI_START_ALIGN = "D:\\College\\lab\\run_1\\lociStartAlign_OUT.csv";
	public static String LOCI_END_ALIGN = "D:\\College\\lab\\run_1\\lociEndAlign_OUT.csv";
	
	public static void main(String[] args) throws IOException{
		
		args = new String[]{"D:\\College\\lab\\run_1\\lociCount.csv", 
				"D:\\College\\lab\\run_1\\lociStart.csv", 
				"D:\\College\\lab\\run_1\\lociEnd.csv", 
				"D:\\College\\lab\\run_1\\lociStartAlignCoords.csv",
				"D:\\College\\lab\\run_1\\lociEndAlignCoords.csv",
		"D:\\College\\lab\\run_1\\lociCountAlign_OUT.csv"};
		foo6pt10(args[0], args[1], args[2], args[3], args[4], args[5]);
	}
	
	public static void foo6pt10(String fociCountCSV, String fociStartCSV, String fociEndCSV, String alignStartCSV, String alignEndCSV, String alignCountCSV) throws IOException {
		int[][] fociCount = parseFile6pt10(251, fociCountCSV);
		int[][] fociStart = parseFile6pt10(251, fociStartCSV);
		int[][] fociEnd = parseFile6pt10(251, fociEndCSV);
		int[][] alignStart = parseFile6pt10(251, alignStartCSV);
		int[][] alignEnd = parseFile6pt10(251, alignEndCSV);
		
		int[][] alignCount = new int[alignStart.length][alignStart[0].length];
		int[][] alignStartNew = new int[alignStart.length][alignStart[0].length];
		int[][] alignEndNew = new int[alignStart.length][alignStart[0].length];
		
		int numPermissed = 0;
		
		for(int i = 0; i < alignCount.length; i++){
			for(int k = 0; k < alignCount[0].length; k++){
				alignCount[i][k] = -1;
				alignStartNew[i][k] = -1;
				alignEndNew[i][k] = -1;
			}
		}
		
//		System.out.println(alignCount.length);
//		System.out.println(alignCount[0].length);
		
		int numWarns = 0;
		
		for(int i = 0; i < alignCount.length; i++){
			
			//breaker catches when everything is zero, meaning end of the array
			//also need something for when there is no alignment??
			int breakerCounter = 0;
			for(int breakr : alignStart[i]){
				if(breakr == 0){
					breakerCounter++;
				}
			}
			//if(breakerCounter >= alignStart[i].length
			if(breakerCounter > 3){
				//end of data, don't want to accidentally analyze empty data set...
				System.out.println("break! at "  + i);
				break;
			}
			
			for(int k = 0; k < alignCount[0].length; k++){
				int a = fociStart[i][k];
				int b = fociEnd[i][k];
				
				for(int r = 0; r < alignCount.length; r++){
					int x = alignStart[r][k];
					int y = alignEnd[r][k];
					
					if((a != 0 || b != 0 || x != 0 || y != 0) 
							&& (a<=x && b>=y && a<=y && b>=x 
							|| a>=x && b<=y && a<=y && b>=x
							|| a<=x && b<=y && b>=x && a<=y
							|| a>=x && b>=y && a<=y && b>=x)){
						alignCount[i][k] = fociCount[r][k];
						alignStartNew[i][k] = fociStart[r][k];
						alignEndNew[i][k] = fociEnd[r][k];
					}
				}
				
				
				//grab closest loci
				if(alignCount[i][k] == -1){
					//System.out.println("using permissive method");
					int permissiveScore = 1;
					while(alignCount[i][k] == -1){
						for(int r = 0; r < alignCount.length; r++){
							int x = alignStart[r][k];
							int y = alignEnd[r][k];
							
							if((b-x + permissiveScore) > 0 && a<=x && a<=y && b<=x && b<=y
									|| (a-y - permissiveScore) < 0 && a>=x && a>=y && b>=x && b>=y){
					              alignCount[i][k] = fociCount[r][k];
					              alignStartNew[i][k] = fociStart[r][k];
					              alignEndNew[i][k] = fociEnd[r][k];
					              numPermissed++;
					             // System.out.println("found " + permissiveScore + " away from position " + alignStart[i][k]);
					              double perc = ((double)permissiveScore) /((double)alignStart[i][k]);
					              //System.out.println(perc);
					              if(perc > 0.01) {
					            	  System.out.println("******** WARNING ********" + perc);
					            	  numWarns++;
					              }
					             // System.out.println("**********" +   alignCount[i][k]);
					              break;
					           }
						}
						//depending on how much you increase permissiveScore by each iteration, determines runtime
						//low increase increases accuracy but increases runtime by quite a lot
						permissiveScore +=1;
					}
				}
			}
			//double d = i / (alignCount.length*14) * 100;
			//System.out.println(d + "%");
		}
		double dudu = numPermissed/ (alignCount.length * 14) *100;
		System.out.println("num permissed: " + numPermissed + "/" + (alignCount.length * 14) + " = " + dudu + "%");
		System.out.println(numWarns);
		System.out.println("writing to file");
		
		try{
			BufferedWriter output = null;
			File file = new File(alignCountCSV);
			output = new BufferedWriter(new FileWriter(file));

			output.write("a1.X,a2.X,a3.X,a4.X,a5.X,a7.X,ab8.X,ab.X,b1.X,b2.X,b3.X,b4.X,b6.X,ore.X,dmel");
			output.newLine();
			
			for(int i = 0; i < alignCount.length; i++){
				for(int k = 0; k < alignCount[0].length; k++){
					output.write(alignCount[i][k] + "");
					//System.out.print(alignCount[i][k] + ",");
					if(k < alignCount[0].length - 1){
						output.write(",");	
					}
				}
				//System.out.println();
				output.newLine();
			}

			output.close();

		}catch(Exception e){
			System.out.println("could not create file");
		}
		
		
		/**
		 * not workgin yet, need to use alignStartNew[][] and alignEndNew[][]
		 * code fininished, neeed testing
		 */
		//new locations, start
		try{
			BufferedWriter output = null;
			File file = new File(LOCI_START_ALIGN);
			output = new BufferedWriter(new FileWriter(file));

			output.write("a1.X,a2.X,a3.X,a4.X,a5.X,a7.X,ab8.X,ab.X,b1.X,b2.X,b3.X,b4.X,b6.X,ore.X,dmel");
			output.newLine();
			
			for(int i = 0; i < alignStart.length; i++){
				for(int k = 0; k < alignStart[0].length; k++){
					output.write(alignStartNew[i][k] + "");
					//System.out.print(alignStart[i][k] + ",");
					if(k < alignStart[0].length - 1){
						output.write(",");	
					}
				}
				//System.out.println();
				output.newLine();
			}

			output.close();

		}catch(Exception e){
			System.out.println("could not create file");
		}
		
		//new locations, end
				try{
					BufferedWriter output = null;
					File file = new File(LOCI_END_ALIGN);
					output = new BufferedWriter(new FileWriter(file));

					output.write("a1.X,a2.X,a3.X,a4.X,a5.X,a7.X,ab8.X,ab.X,b1.X,b2.X,b3.X,b4.X,b6.X,ore.X,dmel");
					output.newLine();
					
					for(int i = 0; i < alignEnd.length; i++){
						for(int k = 0; k < alignEnd[0].length; k++){
							output.write(alignEndNew[i][k] + "");
							//System.out.print(alignStart[i][k] + ",");
							if(k < alignEnd[0].length - 1){
								output.write(",");	
							}
						}
						//System.out.println();
						output.newLine();
					}

					output.close();

				}catch(Exception e){
					System.out.println("could not create file");
				}

	
	}
	
	public static int[][] parseFile6pt10(int r, String PATH) throws IOException{
		
		//hard coding will give problems later
		int[][] rval = new int[r][15];

		FileInputStream inputStream = null;
		Scanner sc = null;
		try {

			inputStream = new FileInputStream(PATH);
			sc = new Scanner(inputStream, "UTF-8");

			String line = "";


			int row = 0;
			while (sc.hasNextLine() && row < 251) {

				try{line = sc.nextLine();}
				catch(Exception e){	/* lol */}

				String[] lineAR = line.split(",");

				
				//has zero means that one or more loci did not align to the thing, sp disrgard the whole line
				boolean hasZero = false;
				for(int i = 0; i < lineAR.length; i++){
					if(lineAR[i].equals("0") && PATH.contains("Align")){
						hasZero = true;
					}
				}
				
				//lineAR[0] is the first individual, used to make sure column labels aren't parsed
				if(lineAR[0].contains(".") || lineAR[0].contains(".") || hasZero){
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
	
	
	
}
