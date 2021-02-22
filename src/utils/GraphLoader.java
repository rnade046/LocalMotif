package utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;

import graph.Interaction;

public class GraphLoader {

	public static void loadInteractionRepository(String inputFile, double fdr){

		
		HashMap<String, Double> interactionToSpectralCountMap = loadSaintReport(inputFile, fdr);
		
//		ArrayList<Interaction> interactionList = new ArrayList<>();
//
//		return interactionList;
	}

	public static HashMap<String,Double> loadSaintReport(String inputRepositoryFile, double fdr) {

		HashMap<String, Double> interactionToSpectralCountMap = new HashMap<String, Double>();// To contain interaction name and ID, and number of occurence
		try {
			/* Read BioGrid ; get all possible human protein-protein interactions */

			InputStream in = new FileInputStream(new File(inputRepositoryFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // header

			while (line != null) { // stops when there are no more lines in file (in)

				String[] col = line.split("\t"); // split line by tabs; obtain individual columns

				if(Double.parseDouble(col[15]) <= fdr) { // col[15] = BFDR
					/* Keep interactions that pass FDR thresholds */ 
					String interactor1 = col[0]; // get name of bait 
					String interactor2 = col[2]; // get name of prey
					
					interactionToSpectralCountMap.put(interactor1+"\t"+interactor2, Double.parseDouble(col[5])); //col[5] = Average Spectral count
				} 
			}
			line = input.readLine(); // read next line

			input.close(); // close BufferedReader
		} catch (IOException e) {
			e.printStackTrace();
		}
		return interactionToSpectralCountMap;
	}

}
