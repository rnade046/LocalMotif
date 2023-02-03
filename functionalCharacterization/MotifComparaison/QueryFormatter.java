package MotifComparaison;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Properties;

public class QueryFormatter {

	public static void main(String[] args) throws FileNotFoundException, IOException {
		
		System.out.println("**Loading parameters file** \n");
		Properties params = new Properties();
		params.load(new FileInputStream(args[0]));		

		String wd = params.getProperty("working_directory");
		String networkName = params.getProperty("network_name");

		int clusteringMeasure = Integer.parseInt(params.getProperty("clusteringMeasure", "0"));
		double percentThreshold = Double.parseDouble(params.getProperty("percentThreshold", "0.2"));

		String clusteringName = "";

		switch(clusteringMeasure) {
		case 0: clusteringName = "_TPD";
		break;
		case 1: clusteringName = "_TPPD_p" + percentThreshold;
		break;
		case 2: clusteringName = "_coreTPD_p" + percentThreshold;
		break;
		}
		
		String projectDirectory =  networkName + clusteringName + "_p" + params.getProperty("significantThreshold");
		
		String height = params.getProperty("height");
		String motifFamilyPrefix = wd + "/motifFamilies/" + projectDirectory + "/Groups_h"+ height + "/" + networkName + clusteringName + "_h" + height + "_ppm_motifFamilyGroup";
		String queryMotifTomTom = wd + "/motifFamilies/"+ projectDirectory + "/Groups_h" + height  + "/" + networkName + clusteringName + "_h"+ height + "_queryMotifsTomTom.txt";
		
		System.out.println("Formatting motif families ppm - all proteins");
		printFormattedQueryMotifsForTomTom(motifFamilyPrefix, Integer.parseInt(params.getProperty("motifFamilyGroups")), queryMotifTomTom);
		
		if(clusteringMeasure == 1 || clusteringMeasure == 2) {
			
			height = params.getProperty("coreHeight");
			motifFamilyPrefix = wd + "/motifFamilies/" + projectDirectory + "/Groups_CoreProteins_h"+ height + "/" + networkName + clusteringName + "_coreProteins_h" + height + "_ppm_motifFamilyGroup";
			queryMotifTomTom = wd + "/motifFamilies/"+ projectDirectory + "/Groups_CoreProteins_h" + height  + "/" + networkName + clusteringName + "_coreProteins_h"+ height + "_queryMotifsTomTom.txt";
			
			System.out.println("Formatting motif families ppm - core proteins");
			printFormattedQueryMotifsForTomTom(motifFamilyPrefix, Integer.parseInt(params.getProperty("coreFamilyGroups")), queryMotifTomTom);
			
		}
	}

	@SuppressWarnings("unused")
	private static void printFormatedQueryMotifsForMotifComp(String familyPrefix, int numFamilies, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			out.write("#INCLUSive Motif Model\n#\n"); // header

			for(int i=1; i<=numFamilies; i++) {

				String familyFile = familyPrefix + i + ".tsv";
				/* Load the ppm (position probability matrix) for every family*/
				double[][] ppm = loadPPM(familyFile);
				out.write("#ID = MotifFamily_" + i + "\n");
				out.write("#W = 8\n");
				
				for(int k=0; k<ppm[0].length; k++) {
					for(int j=0; j<ppm.length; j++) {
						out.write(ppm[j][k] + "\t");
					}
					out.write("\n");
					out.flush();
				}
				out.write("\n");
			}
			out.write("#END\n"); // end signal

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	private static void printFormattedQueryMotifsForTomTom(String familyPrefix, int numFamilies, String outputFile) {
		
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			out.write("MEME version 4\n\n" + "ALPHABET= ACGT\n\n" + 
			"Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n"); // header

			for(int i=1; i<=numFamilies; i++) {
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
