package positionConservation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

public class ReverseComplement {

	public static void main(String[] args) {
		
		String inputFile = "C:\\Users\\rnade046\\Documents\\LESMoNlocal\\analysis\\corrNetTop2-400_coreTPD_p0.4_coreProteins_h0.7_motifFamiliesInfo.tsv";
		String outputFile = "C:\\Users\\rnade046\\Documents\\LESMoNlocal\\analysis\\MotifPosition\\corrNetTop2-400_coreTPD_p0.4_coreProteins_reverse-complement_motifs.txt";
	
		listReverseComplementMotifs(inputFile, outputFile);
	}

	private static void listReverseComplementMotifs(String inputFile, String outputFile) {

		InputStream in;
		try {
			in = new FileInputStream(new File(inputFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));
			
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			
			String line = input.readLine();
			line = input.readLine();
			int count = 1;
			out.write("Family\tReverse-Complement\n");
			
			while(line!=null) {
				
				String motif = line.split("\t")[1];
				out.write(count + "\t" + convertMotif(motif) + "\n");
				
				line = input.readLine();
				count++;
			}
			out.close();
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private static String convertMotif(String motif) {
		String finalMotif = "";
		
		String complement = "";
		for(int i=0; i<motif.length(); i++) {
			complement += mapComplement(motif.charAt(i));
		}
		
		for(int i=motif.length()-1; i>=0; i--) {
			finalMotif += complement.charAt(i);
		}
		
		return finalMotif;
	}

	private static char mapComplement(char c) {

		char n = 'X';

		switch(c) {
		case 'A' : n = 'T';
		break;
		case 'T' : n = 'A';
		break;
		case 'C' : n = 'G';
		break;
		case 'G' : n = 'C';
		break;
		case 'R' : n = 'Y';
		break;
		case 'Y' : n = 'R';
		break;
		case 'D' : n = 'H';
		break;
		case 'H' : n = 'D';
		break;
		case 'V' : n = 'B';
		break;
		case 'B' : n = 'V';
		break;
		case '*' : n = '*';
		break;
		}

		return n; 
	}
}
