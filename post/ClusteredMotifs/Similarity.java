package ClusteredMotifs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

public class Similarity {

	/**
	 * Compute the similarity between all significant motifs
	 * 
	 * @param motifMapOfAnnotatedProteins
	 * @param outputFile
	 * @param distanceOutputFile
	 */
	public static void computeMotifSimilary(HashMap<String, String[]> motifMapOfAnnotatedProteins, String outputFile, String distanceOutputFile) {

		double maxDistance = 0;
		boolean containsMaxValue = false;
		
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			List<String> motifs = new ArrayList<>(motifMapOfAnnotatedProteins.keySet());

			for(int i=0; i<motifs.size(); i++) {
				for(int j=i+1; j<motifs.size(); j++) {

					double distance = computeDistanceSimilarity(motifMapOfAnnotatedProteins.get(motifs.get(i)), motifMapOfAnnotatedProteins.get(motifs.get(j)));
					
					if(distance == Double.MAX_VALUE) {
						containsMaxValue = true;
					}
					
					if(distance != Double.MAX_VALUE && distance > maxDistance) {
						maxDistance = distance;
					}

					out.write(motifs.get(i) + "\t" + motifs.get(j) + "\t" + distance + "\n");
					out.flush();
				}
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		if(containsMaxValue) {
			updateDistanceSimilarity(outputFile, maxDistance, distanceOutputFile);
		}
	}

	public static void updateDistanceSimilarity(String distanceMatrixFile, double maxDistance, String distanceOutputFile) {

		double newDistance = maxDistance + 1;
		
		try {
			InputStream in = new FileInputStream(new File(distanceMatrixFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(distanceOutputFile)));

			String line = input.readLine();
			while(line!=null) {
				
				String[] col = line.split("\t");
				if(Double.parseDouble(col[2]) == Double.MAX_VALUE){
					out.write(col[0] + "\t" +col[1] + "\t" + newDistance + "\n");
				}
				
				line = input.readLine();
			}
			
			out.close();
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	private static double computeDistanceSimilarity(String[] annotatedProteinList1, String[] annotatedProteinList2) {

		HashSet<String> proteinList2 = new HashSet<>(Arrays.asList(annotatedProteinList2)); 

		int countSharedProteins = 0;
		for(String protein: annotatedProteinList1) {
			if(proteinList2.contains(protein)) {
				countSharedProteins++;
			}
		}

		double similarity = countSharedProteins / (double) Math.min(annotatedProteinList1.length, annotatedProteinList2.length);
		double distance = Double.MAX_VALUE;

		if(similarity != 0) {
			distance = 1 / similarity - 1;
		}

		return distance;
	}

}
