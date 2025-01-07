

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

public class GetMotifFamiliesSequenceForBenchmarkComparison {

	public static void main(String[] args) {

		String wd = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/";
		String similarityMatrix = wd + "benchmark/corrNetTop2-400_coreTPD_p0.4_MCL_i2_similarity_coreProts_dec2024.tsv";
		double threshold = 0.5;
		
		String motifFamiliesPPMprefix = wd + "motifFamilies/corrNetTop2-400_coreTPD_p0.4_p6.94181641095822E-7/Groups_CoreProteins_h0.7/corrNetTop2-400_coreTPD_p0.4_coreProteins_h0.7_ppm_motifFamilyGroup";
		String outputPrefix = wd + "benchmark/fasta-mcli2-dec2024/corrNetTop2-400_coreTPD_p0.4_coreProts_TomTomFormatted_mclCluster";
	
		
		/* Iterate through similarity matrix - assess mcl clusters individually */
		InputStream input;
		try {
			input = new FileInputStream(new File(similarityMatrix));
			BufferedReader in = new BufferedReader(new InputStreamReader(input));
			
			String line = in.readLine(); // header
			line = in.readLine();
			int clusterCount = 1;
			while(line!=null) {
				
				System.out.println("Cluster: " + clusterCount);
				
				/* Identify motif families with similarity >= threshold */
				String[] col = line.split("\t"); // lines to ignore [0] = mcl-cluster, [n-1] = max-motifFamily, [n] = max-value
				int cluster = Integer.parseInt(col[0]);
				
				List<Integer> motifsOfInterest = new ArrayList<>();
				
				for(int i=1; i<col.length-2; i++) {
					if(Double.parseDouble(col[i]) >= threshold) {
						motifsOfInterest.add(i);
					}
				}
				
				/* If at least 1 family identified - print it's corresponding sequences to file */
				if(motifsOfInterest.size() >= 1) {
					printFormattedQueryMotifsForTomTom(motifFamiliesPPMprefix, motifsOfInterest, outputPrefix + cluster);
				}
				line = in.readLine();
				clusterCount++;
			}
			
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}


	private static void printFormattedQueryMotifsForTomTom(String familyPrefix, List<Integer> motifsOfInterest, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			out.write("MEME version 4\n\n" + "ALPHABET= ACGT\n\n" + 
					"Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n"); // header

			for(Integer i : motifsOfInterest) {
				System.out.print(i + ".");

				String familyFile = familyPrefix + i + ".tsv";
				/* Load the ppm (position probability matrix) for every family*/
				double[][] ppm = loadPPM(familyFile);
				out.write("MOTIF Family_" + i + "\n");
				out.write("letter-probability matrix: alength= 4 w= 8\n");

				for(int k=0; k<ppm[0].length; k++) {
					for(int j=0; j<ppm.length; j++) {
						out.write(ppm[j][k] + "\t");
					}
					out.write("\n");
					out.flush();
				}
				out.write("\n");
			}
			out.write("\n"); // end signal

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		System.out.println("Done");
	}

	private static double[][] loadPPM(String inputFile) {

		double[][] ppm = new double[4][8];

		try {
			InputStream in = new FileInputStream(new File(inputFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();
			int lineCount = 0;

			while(line != null) {
				String[] col = line.split("\t");

				for(int i=0; i<col.length; i++) {
					ppm[lineCount][i] = Double.parseDouble(col[i]);
				}

				line = input.readLine();
				lineCount++;
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return ppm;
	}
}
