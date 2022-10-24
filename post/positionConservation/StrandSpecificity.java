package positionConservation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;

import org.apache.commons.math3.distribution.NormalDistribution;

public class StrandSpecificity {

	public static void main(String[] args) {

		String wd = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/";

		String fwdFasta = wd + "/MotifPosition/corrNetTop2_3UTRlongestSequences.txt";
		String revCFasta = wd + "/MotifPosition/corrNetTop2_reverse-complement-sequences-3UTR.txt";

		String motifsFile = wd + "corrNetTop2-400_coreTPD_p0.4_coreProteins_h0.7_motifFamiliesInfo.tsv";

		String output = wd + "/MotifPosition/coreTPD0.4/StrandSpecificity_corrNet2-400_coreTPD_p0.4_3UTR_MotifFamilies_h0.7_oct.tsv";

		
		determineStrandSpecificity(fwdFasta, revCFasta, motifsFile, output);
		
	}

	private static int countMotifOccurrence(String fastaFile, HashSet<String> motifInstances) {

		int motifCount=0;
		int testLength=2000;
		int motifLength=8;
		/* Check for motif in all fasta of the filtered fasta file */
		InputStream in;
		try {
			in = new FileInputStream(new File(fastaFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();

			while(line != null) {

				/* load sequence; ignore lines that are identifiers*/
				if(line.startsWith(">")) {

					line = input.readLine();
					
					/* format sequence */
					String sequence = SequenceUtils.formatSequence(line, testLength);

					int half = (int) Math.floor(testLength / 2);
					/* search for motif instances in first 500 nucleotides */
					for(int i=0; i<half-motifLength; i++) {

						String substring = sequence.substring(i, i+motifLength);

						/* if current substring is a motif instance, increase motif position count */
						if(motifInstances.contains(substring)) {
							motifCount++;
						}
					}

					/* search for motif instances in last 500 nucleotides */
					for(int i=half; i<testLength-motifLength; i++) {

						String substring = sequence.substring(i, i+motifLength);

						/* if current substring is a motif instance, increase motif position count */
						/* if current substring is a motif instance, increase motif position count */
						if(motifInstances.contains(substring)) {
							motifCount++;
						}
					}
				}

				line = input.readLine();
			}
			input.close();

		} catch (IOException e) {
			e.printStackTrace();
		}


		return motifCount;
	}

	private static void determineStrandSpecificity(String fwdFasta, String revCFasta, String motifsFile, String outputFile) {

		//HashMap<String, Double> specificityMap = new HashMap<>();
		
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			out.write("motif\tmotifFamily\tFWD-occurrence\tRevC-Occurence\tStrand-Specificity\tProbability\n");
			
			/* load motifs to test */
			HashMap<String, Integer> motifsToTest = getMotifs(motifsFile);
			
			System.out.println("Motif to test: " + motifsToTest.size());
			int motifCount = 0;
			for(Entry<String, Integer> e: motifsToTest.entrySet()) {
				
				String motif = e.getKey();
				
				motifCount++;
				System.out.println(motifCount);
				
				/* determine possible motifs */
				HashSet<String> motifInstances = SequenceUtils.getPossibleMotifInstances(motif);

				/* Motif occurrence in FWD */
				int fwdMotifCount = countMotifOccurrence(fwdFasta, motifInstances);
				
				/* Motif occurrence in revC*/
				int revCMotifCount = countMotifOccurrence(revCFasta, motifInstances);

				/* Strand Specificity */
				double ss = fwdMotifCount/(double) revCMotifCount;
				
				/* PARAMS for binomial distribution */
				int trials = fwdMotifCount + revCMotifCount;
				double p = 0.5;
				
				/* normal distribution approximation */
				double mean = trials * p;
				double sdev = Math.sqrt(trials * p * (1-p));
				
				NormalDistribution nd = new NormalDistribution(mean, sdev);
				double probability = nd.probability(fwdMotifCount-1, trials);
				
				/* output info */
				out.write(motif + "\t" + e.getValue() + "\t" + fwdMotifCount + "\t" + revCMotifCount + "\t" + ss + "\t" + probability + "\n");
				out.flush();
				
				//specificityMap.put(motif, probability);
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		//return specificityMap;
	}
	
	

	@SuppressWarnings("unused")
	private static double sumOfMotifOccurences(String inputFile) {

		double sum = 0;
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inputFile))));

			String line = in.readLine();
			while(line!=null) {
				sum += Double.parseDouble(line.split("\t")[1]);
				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return sum;
	}

	@SuppressWarnings("unused")
	private static void printSrandSpecificities(String outputFile, double[] specificityAll, double[] specificityCore) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			out.write("motif#\tSS-all\tSS-core\n");
			for (int i=0;i < specificityAll.length; i++) {
				out.write((i+1) + "\t" + specificityAll[i] + "\t" + specificityCore[i] +"\n");
				out.flush();
			}

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
	
	private static HashMap<String, Integer> getMotifs(String motifFamilyFile){
		
		HashMap<String, Integer> motifsToTestMap = new HashMap<>();
		
		InputStream in;
		try {
			in = new FileInputStream(new File(motifFamilyFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // header
			line = input.readLine();

			while(line != null) {
				motifsToTestMap.put(line.split("\t")[0], Integer.parseInt(line.split("\t")[1]));

				line = input.readLine();
			}

			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return motifsToTestMap;
	}
	
	@SuppressWarnings("unused")
	private static void testND() {
		
		int trials = 1000;
		double p = 0.5;
		
		/* normal distribution approximation */
		double mean = trials * p;
		double sdev = Math.sqrt(trials * p * (1-p));
		
		int fwdCount = 700;
		
		NormalDistribution nd = new NormalDistribution(mean, sdev);
		double probability = 1- nd.probability(0, fwdCount);
	
		
		System.out.println("test : " + probability);
	}
	

}
