package positionConservation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

import org.apache.commons.math3.distribution.BinomialDistribution;

public class StrandSpecificity {
	
	public static void main(String[] args) {
	
		String wd = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/MotifPosition/coreTPD0.4/";
		
		String fwdFile = wd + "corrNet2-400_coreTPD_p0.4_3UTR_FWD_motifPositions_allSeq_Normalized_motif";
		String fwdFile2 = wd + "corrNet2-400_coreTPD_p0.4_3UTR_FWD_motifPositions_coreProts_Normalized_motif";
		String rvcFile = wd + "corrNet2-400_coreTPD_p0.4_3UTR_RevC_motifPositions_allSeq_Normalized_motif";
		int numMotifs = 31;
		
		String output = wd + "StrandSpecificity_corrNet2-400_coreTPD_p0.4_3UTR.tsv";
		
		double[] ss = determineStrandSpecificity(fwdFile, rvcFile, numMotifs);
		double[] ss2 = determineStrandSpecificity(fwdFile2, rvcFile, numMotifs);
		
		
		/* print results */
		printSrandSpecificities(output, ss, ss2);
		
	}

	private static double[] determineStrandSpecificity(String fwdFile, String revComplFile, int numMotifs) {
		
		/* for each motif; determine strand specificity */
		
		double[] strandSpecificty = new double[numMotifs];
		double [] significance = new double[numMotifs];
		
		for(int i=1; i<= numMotifs; i++) {
			System.out.println("motif : " + i);

			/* compute sum of occurrence in FWD strand */
			double fwdSum = sumOfMotifOccurences(fwdFile + i);
		
			/* compute sum of occurrence in RevCompl Strand */
			double rvcSum = sumOfMotifOccurences(revComplFile + i);
			
			/* determine specificity */
			double specificity = fwdSum/rvcSum;
			strandSpecificty[i-1] = specificity;
			
			/* determine significance */
			double trials = fwdSum + rvcSum;
			BinomialDistribution bd = new BinomialDistribution((int)trials, 0.5);
			
		}
		
		return strandSpecificty;
		
	}
	
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
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
}
