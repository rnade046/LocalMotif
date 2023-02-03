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

public class CheckNetworkTopology {

	public static void main(String[] args) {
		
		/* Check distribution of shortest paths 
		 * 	Load dm; get top triangle; print values to file; histogram in R */
		String distanceMatrixFile = "C:\\Users\\Rachel\\Documents\\LESMoNlocal\\analysis\\corrNetTop2_removedOverConnectedProteins_300_distanceMatrix2.txt";
		String pathsFile = "C:\\Users\\Rachel\\Documents\\LESMoNlocal\\analysis\\corrNetTop2_removedOverConnectedProteins_300_paths.txt";
		
		System.out.println("Load distance matrix");
		ArrayList<Double> allPaths = loadPaths(distanceMatrixFile);
		
		System.out.println("Compute Binning");
		int[] pathsBinning = computeBinning(allPaths);
		
		System.out.println("Print binning");
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
				
				if(index%100 == 0) {
					System.out.print(index + ".");
				}
				if(index%1000 == 0) {
					System.out.println();
				}
				
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
			System.out.println("Done\n");
		} catch (IOException e) {
			e.printStackTrace();
		}
		return allPaths;
 	}
	
	private static int[] computeBinning(ArrayList<Double> allPaths){
		
		int maxBin = (int) Math.ceil(Collections.max(allPaths)); // find max bin 
		int[] pathsBinning = new int[maxBin];
		System.out.println("Max bin: " + maxBin);
		for(int i=0; i< allPaths.size(); i++) {
			
			if(i%100 == 0) {
				System.out.print(i + ".");
			} 
			if(i%1000 == 0) {
				System.out.println();
			}
			
			int idx = (int) Math.ceil(allPaths.get(i)); // find 
			pathsBinning[idx-1] += 1; // increase number of occurrence by 1
		}
		System.out.println("Done");
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
