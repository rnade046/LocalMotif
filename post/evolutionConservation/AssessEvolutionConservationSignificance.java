package evolutionConservation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

//import org.apache.commons.math3.distribution.NormalDistribution;

public class AssessEvolutionConservationSignificance {

	public static void main(String[] args) {

		String motifFreqFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/evolutionConservation/coreTPD0.4_FeatureBitsOutput_core_motifFreq.out";
		String motifConservationFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/evolutionConservation/coreTPD0.4_FeatureBitsOutput_core_phastCons30_motifCons.out";
		
		String outputFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/evolutionConservation/coreTPD0.4_evolutionConservationSignificance_core_phastCons30_sept21.tsv";

		//double prob = 0.218524153; // phastCons 20
		//double prob = 0.22523417; // phastCons 30
		double prob = 0.27063104254; 
		int motifs = 31;

		assessEvolutionSignificance(motifFreqFile, motifConservationFile, motifs, prob, outputFile);
		//testNormPvalue();
	}


	private static void assessEvolutionSignificance(String motifFreqFile, String motifConservationFile, int motifs, double prob, String outputFile) {

		List<String[]> motifSignificance = new ArrayList<>();

		/* load motif frequency */
		int[] motifFreqs = loadFeatureBitOutput(motifFreqFile, motifs);
		
		/* load motif conservation */
		int[] motifCons = loadFeatureBitOutput(motifConservationFile, motifs);
		
		for(int i=0; i< motifs ; i++) {
			if(motifFreqs[i] > 0 && motifCons[i] > 0) {
				
				double mean = motifFreqs[i] * prob;
				double sdev = mean* (1-prob); 

				//NormalDistribution nd = new NormalDistribution(mean, sdev);
				//double pval = nd.probability(motifCons[i], motifFreqs[i]);
				double pval = computeNormPvalue(mean, sdev, motifCons[i], motifFreqs[i]);
				
				String[] values = new String[] {Double.toString(motifCons[i]), Double.toString(motifFreqs[i]), Double.toString(motifCons[i]/(double) motifFreqs[i]), Double.toString(pval)};
				motifSignificance.add(values);
				
			} else {
				String[] values = new String[] {Double.toString(motifCons[i]), Double.toString(motifFreqs[i]), "NA", "0"};
				motifSignificance.add(values);
			}
		}
		printMotifSignificance(motifSignificance, outputFile);
	}

	private static int[] loadFeatureBitOutput(String inputFile, int motifs) {

		int[] values = new int[motifs];
		
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inputFile))));
			
			String line = in.readLine();
			int count = 0;
			while(line!=null) {
				
				String terms[] = line.split("\\s+");
				values[count] = Integer.parseInt(terms[1]);
				
				line = in.readLine();
				count++;
			}
			
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return values;
	}

	private static void printMotifSignificance(List<String[]> motifSignificance, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			out.write("motif#\tConservation\tFrequency\tFold-change\tp-val\n");

			for(int i=0; i<motifSignificance.size(); i++) {
				out.write((i+1) + "\t");

				String[] values = motifSignificance.get(i);
				for(int j=0; j<values.length; j++) {
					out.write(values[j] + "\t");
				}
				out.write("\n");
				out.flush();
			}

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	/* test that this works !! */
	public static double computeNormPvalue(double mean, double stdev, double value, double max){
		double prob = 0.0;
		
		for(double j = value; j <= max; j=j+0.0001){
			
			double s2Pi = Math.sqrt(2*Math.PI);
			double sd = Math.sqrt(stdev);
			double exp = (-((j-mean)*(j-mean)))/(2*stdev);
			double add = ((1/(sd*s2Pi))*Math.pow(Math.E,exp)*0.0001);
			prob = prob + add;
		}
		
		return prob;
		
	}

	
	public static void testNormPvalue(){
		double prob = 0.0;
		double mean = 20; 
		double stdev = 4;
		double value = 28;
		double max = 34;
		
		for(double j = value; j <= max; j=j+0.001){
			
			double s2Pi = Math.sqrt(2*Math.PI);
			double sd = Math.sqrt(stdev);
			double exp = (-((j-mean)*(j-mean)))/(2*stdev);
			double add = ((1/(sd*s2Pi))*Math.pow(Math.E,exp)*0.001);
			prob = prob + add;
		}
		
		System.out.println(prob);
		//return prob;
		
	}

	
	
	
//	private static void assessSignification (String inputFile, int motifs, double prob, String outputFile) {
//
//		List<String[]> motifSignificance = new ArrayList<>();
//
//		for(int i=1; i<=motifs; i++) {
//			System.out.println(i + ".");
//			Double phastCons = 0.0;
//			Double motifCons = 0.0;
//
//			try {
//				BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inputFile))));
//
//				String line = in.readLine();
//
//				outerloop:
//					while (line!= null) {
//
//						/* search for values corresponding to current motif */
//						if(line.length() <= 2 && line.startsWith(String.valueOf(i))) {
//
//
//							/* when at current motif; obtain values */
//							line = in.readLine();
//							int lineCount = 1;
//
//							while(line !=null && lineCount <= 4){ // line 1 = phastCons, 2 = value, 3 = UTR, 4 = value 
//
//								/* check that info is not incomplete (ie. switch to the next motif) */
//								if(line.length() > 2) {
//
//									switch(lineCount){
//									case 2: 
//										String[] terms = line.split("\\s+");
//										phastCons = Integer.parseInt(terms[0]) / (double) Long.parseLong(terms[3]);
//										break;
//									case 4:
//										String[] terms2 = line.split("\\s+");
//										motifCons = Integer.parseInt(terms2[0]) / (double) Long.parseLong(terms2[3]);
//										break;
//									}
//								} else {
//									break outerloop;
//								}
//
//								line = in.readLine();
//								lineCount++;
//							}
//						}
//						line = in.readLine();
//					}
//				in.close();
//			} catch (IOException e) {
//				e.printStackTrace();
//			}
//
//
//		}
//
//		printMotifSignificance(motifSignificance, outputFile);
//	}
}
