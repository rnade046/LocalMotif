package utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;

public class AssessEnrichment {

	public static void assessSignificanceScores(String degenMotifsClusteringFilesPrefix, int numOfDegenMotifFiles, String outputFile) {
		
		/* initialize list */
		ArrayList<Double> significantScoresList = new ArrayList<>();
		
		for(int i=0; i<numOfDegenMotifFiles; i++) {
			String degenMotifClusteringFile = degenMotifsClusteringFilesPrefix + i;
			
			File f = new File(degenMotifClusteringFile);
			if(f.exists()) {
				getSignificantScores(degenMotifClusteringFile, significantScoresList);
			} else {
				System.out.println("No file: " + i + "\n");
			}
			
		} 
		
		printSignifcanceScores(significantScoresList, outputFile);
		System.out.println("Done");
	}
	
	private static void getSignificantScores(String inputFile, ArrayList<Double> significanceScoresList){
		
		try {
			InputStream in = new FileInputStream(new File(inputFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // no header
			
			while(line!=null) {
				double pval = Double.parseDouble(line.split("\t")[3]);
				significanceScoresList.add(pval);
				
				line = input.readLine();
			}
			
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private static void printSignifcanceScores(ArrayList<Double> significantScoresList, String outputFile) {
		
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			
			for(double pval: significantScoresList) {
				out.write(pval + "\n");
				out.flush();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
