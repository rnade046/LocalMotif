/**
 * 
 */
package opt;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

public class CheckNetworkTopology {

	public static void main(String[] args) {
		
		/* Check distribution of shortest paths 
		 * 	Load dm; get top triangle; print values to file; histogram in R */
		String distanceMatrixFile = "C:\\Users\\Rachel\\Documents\\LESMoNlocal\\analysis\\corrNetTop2_distanceMatrix2.txt";
		String pathsFile = "C:\\Users\\Rachel\\Documents\\LESMoNlocal\\analysis\\corrNetTop2_listOfShortestPathsBinning.tsv";
		
		ArrayList<Double> allPaths = loadPaths(distanceMatrixFile);
		int[] pathsBinning = computeBinning(allPaths);
		printBinning(pathsBinning, pathsFile);
		
	}
	
	private static ArrayList<Double> loadPaths(String distanceMatrixFile) {
		
		/* Load top triangle of distance matrix */
		ArrayList<Double> allPaths = new ArrayList<>();
				
		try {
			InputStream in = new FileInputStream(new File(distanceMatrixFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));
			
			String line = input.readLine();
			int index = 1;
			while(line != null) {
				String[] col = line.split("\t");
				
				/* Add elements to list */
				
				for(int i=index; i<col.length; i++) {
					
					if(Double.parseDouble(col[i]) != 0) {
					
						allPaths.add(Double.parseDouble(col[i]));
					}
				}
				
				line = input.readLine();
				index++;
			}
	
			input.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		return allPaths;
 	}
	
	private static int[] computeBinning(ArrayList<Double> allPaths){
		
		int maxBin = (int) Math.ceil(Collections.max(allPaths)); // find max bin 
		int[] pathsBinning = new int[maxBin];

		for(Double path : allPaths) {
			int idx = (int) Math.ceil(path); // find 
			pathsBinning[idx-1] += 1; // increase number of occurrence by 1
		}
		
		return pathsBinning;
	}

	private static void printBinning(int[] pathsBinning, String outputFile) {
		
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			
			out.write("Bin\tOccurrence\n");
			
			for(int i=0; i<pathsBinning.length; i++) {
				out.write(i + "-" + (i+1) + "\t" + pathsBinning[i] + "\n");
			}
			
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
}
