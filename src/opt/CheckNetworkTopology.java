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

public class CheckNetworkTopology {

	public static void main(String[] args) {
		
		/* Check distribution of shortest paths 
		 * 	Load dm; get top triangle; print values to file; histogram in R */
		String distanceMatrixFile = "C:\\Users\\Rachel\\Documents\\LESMoNlocal\\analysis\\corrNetTop2_distanceMatrix2.txt";
		String pathsFile = "C:\\Users\\Rachel\\Documents\\LESMoNlocal\\analysis\\corrNetTop2_listOfShortestPaths.tsv";
		
		checkShortestPathDistribution(distanceMatrixFile, pathsFile);
		
	}
	
	private static void checkShortestPathDistribution(String distanceMatrixFile, String outputFile) {
		
		/* Load top triangle of distance matrix */
		
		try {
			InputStream in = new FileInputStream(new File(distanceMatrixFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));
			
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			
			String line = input.readLine();
			
			while(line != null) {
				String[] col = line.split("\t");
				
				/* Add elements to list */
				int index = 1;
				for(int i=index; i<col.length; i++) {
					out.write(col[i] + "\n");
				}
				
				out.flush();
				line = input.readLine();
				index++;
			}
	
			input.close();
			out.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		}

 	}
}
