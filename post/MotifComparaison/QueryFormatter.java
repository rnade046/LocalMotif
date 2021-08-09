package MotifComparaison;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

public class QueryFormatter {

	public static void main(String[] args) {

		String motifFamilyPrefix = "C:\\Users\\Rachel\\Documents\\LESMoNlocal\\analysis\\motifFamilies\\ppm\\corrNetTop2_ppm_motifFamilyGroup";
		String queryMotifOutput = "C:\\Users\\Rachel\\Documents\\LESMoNlocal\\analysis\\motifFamilies\\ppm\\queryMotifs.txt";
		
		printFormatedQueryMotifs(motifFamilyPrefix, 10, queryMotifOutput);

	}

	private static void printFormatedQueryMotifs(String familyPrefix, int numFamilies, String outputFile) {

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
