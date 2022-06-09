package positionConservation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashSet;

public class ReverseComplement {

	public static void main(String[] args) {

		String inputFile ="C:\\Users\\rnade046\\Documents\\LESMoNlocal\\analysis\\MotifPosition\\corrNetTop2_longestSequence.tsv";
		String fastaFile = "C:\\Users\\rnade046\\Documents\\LESMoNlocal\\analysis\\MotifPosition\\corrNetTop2_3UTRlongestSequences.txt";
		String outputFile = "C:\\Users\\rnade046\\Documents\\LESMoNlocal\\analysis\\MotifPosition\\corrNetTop2_reverse-complement-sequences-3UTR.txt";

		generateReverseComplements(inputFile, fastaFile, outputFile);
		//listReverseComplementMotifs(inputFile, outputFile);
	}

	private static void generateReverseComplements(String refSeqIdsFile, String fastaFile, String outputFile) {

		/* obtain refseq of sequences to consider */
		HashSet<String> refSeqIdSet = getRefSeqIdSet(refSeqIdsFile);
		System.out.println("Sequences to find: " + refSeqIdSet.size());

		/* load sequences associated to ids in the fasta file */
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(fastaFile))));
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			
			String line = in.readLine(); // no header 
			
			boolean storeSeq = false;
			String seq = "";
			String id = "";
			
			int motifsFound = 0;
			
			while(line!=null) {
				
				if(storeSeq & !line.startsWith(">")) {
					seq += line;
				}
				
				if(line.startsWith(">")) { // new entry, check id
					
					/* if a sequence is loaded */
					if(!seq.isEmpty()) {
						
						/* generate reverse complement */
						String revC = convertSequence(seq);
						
						/* output reverse complement */
						out.write(">RevComplement_" + id + "\n");
						out.write(revC + "\n");
						out.flush();
						
						/* re-initialize variables */
						seq = "";
						id = "";
						storeSeq = false;
						
						if(motifsFound == refSeqIdSet.size()) {
							break;
						}
					}
					
					/* get current id */
					String[] col = line.split("_|\\.");
					//id = col[2] + "_" + col[3]; // col[2] = type of ID (eg. NM) ; col[3] = number ID (#####)
					id = col[1] + "_" + col[2];
					/* if id is in our list, store it's sequence to obtain it's reverse complement */
					if(refSeqIdSet.contains(id)) {
						storeSeq = true;
						motifsFound++;
						
						System.out.print(motifsFound + ".");
						
						if(motifsFound%50 == 0) {
							System.out.println();
						}
					}	
				}
				line = in.readLine();
			}
	
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private static HashSet<String> getRefSeqIdSet(String inputFile){

		HashSet<String> refSeqIdSet = new HashSet<>();
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(inputFile))));

			String line = in.readLine(); // no header
			while(line != null) {

				refSeqIdSet.add(line.split("\t")[1]); // [1] = refSeqId
				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return refSeqIdSet;
	}


	@SuppressWarnings("unused")
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
				out.write(count + "\t" + convertSequence(motif) + "\n");

				line = input.readLine();
				count++;
			}
			out.close();
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private static String convertSequence(String motif) {
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

		char n = ' ';

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
