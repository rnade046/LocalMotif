package MotifComparaison;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;

public class MiRRNAFormatter {
	
	public static void main(String[] args) {
		
		String mirDB = "C:\\Users\\rnade046\\Documents\\LESMoNlocal\\analysis\\motifComparison\\mature-mirRNA.fa";
		String outputFile = "C:\\Users\\rnade046\\Documents\\LESMoNlocal\\analysis\\motifComparison\\mirRNA-mature-TomTomFormatted.txt";
		String outputFile2 = "C:\\Users\\rnade046\\Documents\\LESMoNlocal\\analysis\\motifComparison\\mirRNA-mature-revComp-TomTomFormatted.txt";
		
		
		formatMirDBForTomTom(mirDB, outputFile);
		formatMirDBRevComplementforTomTom(mirDB, outputFile2);
		
	}
	
	private static void formatMirDBForTomTom(String mirDBFile, String outFile) {
		
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(mirDBFile))));
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outFile)));
			
			out.write("MEME version 4\n\n" + "ALPHABET= ACGU\n\n" + 
					"Background letter frequencies\nA 0.25 C 0.25 G 0.25 U 0.25\n\n"); // header
			
			String line = in.readLine();
			while(line!=null) {

				if(line.startsWith(">hsa")) {
					String motif = line.split("\\s+")[1];
							
					out.write("MOTIF " + motif +"\n");
					
					line = in.readLine(); // sequence
					out.write("letter-probability matrix: alength= 4 w= " + line.length() +"\n");
					
					ArrayList<ArrayList<Integer>> ppm = getPPM(line);
					for(int i=0; i<ppm.size(); i++) {
						for(int j=0; j<ppm.get(i).size(); j++) {
							out.write(ppm.get(i).get(j) + " ");
						}
						out.write("\n");
					}
					out.write("\n");
					out.flush();
					
				}
				line=in.readLine();
			}
			out.write("#END\n"); // end signal
			out.flush();
			
			in.close();
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	private static void formatMirDBRevComplementforTomTom(String mirDBFile, String outFile) {
		
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(mirDBFile))));
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outFile)));
			
			out.write("MEME version 4\n\n" + "ALPHABET= ACGU\n\n" + 
					"Background letter frequencies\nA 0.25 C 0.25 G 0.25 U 0.25\n\n"); // header
			
			String line = in.readLine();
			while(line!=null) {

				if(line.startsWith(">hsa")) {
					String motif = line.split("\\s+")[1];
							
					out.write("MOTIF " + motif +"\n");
					
					line = in.readLine(); // sequence
					out.write("letter-probability matrix: alength= 4 w= " + line.length() +"\n");
					
					String revComplement = getReverseComplement(line);
					
					ArrayList<ArrayList<Integer>> ppm = getPPM(revComplement);
					for(int i=0; i<ppm.size(); i++) {
						for(int j=0; j<ppm.get(i).size(); j++) {
							out.write(ppm.get(i).get(j) + " ");
						}
						out.write("\n");
					}
					out.write("\n");
					out.flush();
					
				}
				line=in.readLine();
			}
			out.write("#END\n"); // end signal
			out.flush();
			
			in.close();
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	private static ArrayList<ArrayList<Integer>> getPPM(String seq) {
		ArrayList<ArrayList<Integer>> ppm = new ArrayList<>();
		
		for(int i=0; i<seq.length(); i++) {
			
			ArrayList<Integer> position = new ArrayList<>();
			switch(seq.charAt(i)) {
			case 'A': Collections.addAll(position, 1, 0, 0, 0);
			break;
			case 'C': Collections.addAll(position, 0, 1, 0, 0);
			break;
			case 'G': Collections.addAll(position, 0, 0, 1, 0);
			break;
			case 'U': Collections.addAll(position, 0, 0, 0, 1);
			break;
			}
			
			ppm.add(position);
		}
		return ppm;
	}
	
	private static String getReverseComplement(String seq) {
	
		String complement = "";
		
		for(int i=seq.length()-1; i>=0; i--) {
			switch(seq.charAt(i)) {
			case 'A' : complement += 'U';
			break;
			case 'U' : complement += 'A';
			break;
			case 'C' : complement += 'G';
			break;
			case 'G' : complement += 'C';
			break;
			}
		}
		
		return complement;
	}
}
