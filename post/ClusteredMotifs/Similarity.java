package clusteredMotifs;

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
	public static void computeMotifSimilary(HashMap<String, String[]> motifMapOfAnnotatedProteins, String motifsOutputFile, String distanceOutputFile) {

		double[][] distanceMatrix = new double[motifMapOfAnnotatedProteins.size()][motifMapOfAnnotatedProteins.size()];

		double maxDistance = 0;
		boolean containsMaxValue = false;

		/* initialize distance matrix */
		List<String> motifs = new ArrayList<>(motifMapOfAnnotatedProteins.keySet());
		printMotifOrder(motifs, motifsOutputFile);
		for(int i=0; i<motifs.size(); i++) {
			
			if(i%10 == 0) {
				System.out.print(i + ".");
			}

			if(i%100 == 0) {
				System.out.println();
			}

			for(int j=i+1; j<motifs.size(); j++) {

				double distance = computeDistanceSimilarity(motifMapOfAnnotatedProteins.get(motifs.get(i)), motifMapOfAnnotatedProteins.get(motifs.get(j)));

				distanceMatrix[i][j] = distance;
				distanceMatrix[j][i] = distance;
				
				if(distance == Double.MAX_VALUE) {
					containsMaxValue = true;
				}

				if(distance != Double.MAX_VALUE && distance > maxDistance) {
					maxDistance = distance;
				}
				
			}	

		}
		
		if(containsMaxValue) {
			System.out.println("**Updating max distance**");
			
			double newDistance = maxDistance + 1;
			
			for(int i=0; i<motifs.size(); i++) {
				
				if(i%10 == 0) {
					System.out.print(i + ".");
				}

				if(i%100 == 0) {
					System.out.println();
				}
				
				for(int j=i+1; j<motifs.size(); j++) {
					if(distanceMatrix[i][j] == Double.MAX_VALUE) {
						distanceMatrix[i][j] = newDistance;
						distanceMatrix[j][i] = newDistance;
					}
				}
			}
		}

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(distanceOutputFile)));
			
			for(int i=0; i<motifs.size(); i++) {
				for(int j=0; j<motifs.size(); j++) {
					out.write(distanceMatrix[i][j] + "\t");
				}
			out.write("\n");
			out.flush();
			}

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
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

	private static void printMotifOrder(List<String> motifs, String outputFile) {
		
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			
			for(int i=0; i<motifs.size(); i++) {
				out.write(motifs.get(i) + "\n");
				out.flush();
			}

			out.close();
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
