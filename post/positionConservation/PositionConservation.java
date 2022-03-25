package positionConservation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

public class PositionConservation {

	private String fasta;
	private int motifLength;
	private HashMap<String, Integer> fastaIdx;
	private HashMap<String, String[]> proteinInfoMap;
	private HashSet<String> annotatedProteins;

	private Integer motifFoundCount = 0;
	// Constructor
	public PositionConservation(String fastaFile, String proteinInfoFile, String protFreqFile, int l) {

		this.fasta = fastaFile;
		this.motifLength = l;
		System.out.println("Generate fasta file index");
		this.fastaIdx = generateFastaFileIdx(fastaFile); // generate index for FASTA file
		System.out.println("Load protein info file");
		this.proteinInfoMap = loadProteinInfoFile(proteinInfoFile); // load motif = refSeqId1|2|..|n
		this.annotatedProteins = loadAnnotatedProteinsInNetwork(protFreqFile); // load {protein1 | prot2 |..| n} proteins in network
		// can remove annotatedProteins -- annotation subset files are generate with only the annotated proteins
	}

	public void getMotifPositions(String representativeMotifsFile, String extractedAnnotationsFile, String motifOutputPrefixFile) {

		/* Load representative motifs to test */ 
		HashMap<String, Integer> motifMap = loadRepresentativeMotifs(representativeMotifsFile);

		/* Load significant motif and its annotated proteins from extracted annotation 1 at a time */ 
		InputStream in;
		try {
			in = new FileInputStream(new File(extractedAnnotationsFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();
			int motifCount = 0;
			while(line != null) {

				String motif = line.split("\t")[0];

				if(motifMap.containsKey(motif)) {
					motifCount ++;
					System.out.println("Testing motif : " + motifCount);

					/* Determine possible instance motifs */
					HashSet<String> possibleInstances = getPossibleMotifInstances(motif);
					
					int[] motifPositions = new int[1000];
					this.motifFoundCount = 0;
					
					int protCount = 0;
					String[] proteins = line.split("\t")[2].split("\\|");
					System.out.println("proteins to check: " + proteins.length);

					for(String prot: proteins) {
						if(this.annotatedProteins.contains(prot) && this.proteinInfoMap.containsKey(prot)) {

							/* obtain longest 3'UTR sequence corresponding to protein */
							String longestSeq = getLongestSequence(prot);

							/* format sequence to fit 1000 nucleotides */
							String finalSeq = formatSequence(longestSeq);

							/* search for motif positions */
							motifPositions = getMotifPositions(possibleInstances, finalSeq, motifPositions);
						}
						protCount++;
						System.out.print(protCount + ".");
						if(protCount % 50 == 0) {
							System.out.println();
						}
					}
					System.out.println("Done");
					System.out.println("Sequences containing motifs in first/last 500 nucleotides: " + this.motifFoundCount);
					
					/* Print positions */
					String motifOutputFile = motifOutputPrefixFile+ "motif" + motifMap.get(motif);
					printMotifPosition(motifPositions, motifOutputFile);
				}
				line = input.readLine();
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	private static HashMap<String, Integer> loadRepresentativeMotifs(String inputFile){

		HashMap<String, Integer> motifSet = new HashMap<>();

		InputStream in;
		try {
			in = new FileInputStream(new File(inputFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // header
			line = input.readLine();

			while(line != null) {

				motifSet.put(line.split("\t")[1], Integer.parseInt(line.split("\t")[0])); // [1] = motif, [0] = motifCount
				line = input.readLine();
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return motifSet;
	}

	private String getLongestSequence(String prot) {

		String longestSeq = "";
		int seqLength = 0; 

		String[] refSeqIDs = this.proteinInfoMap.get(prot);

		/* For a given protein, check it's various sequences (corresponding to different refSeqIds) */
		for(String id: refSeqIDs) {
			if(this.fastaIdx.containsKey(id)) {

				/* Load sequence */
				int indexOfFasta = this.fastaIdx.get(id);
				String currentSeq = loadRNAsequence(indexOfFasta);

				/* Check if sequence is longer than current stored */
				if(currentSeq.length() > seqLength) {
					longestSeq = currentSeq;
					seqLength = currentSeq.length();
				}
			}
		}
		return longestSeq;
	}

	private static String formatSequence(String seq) {
		String finalSeq = "";

		int l = seq.length();

		if(l == 1000) {
			finalSeq = seq;

		} else if(l < 1000) {

			String fill = "";
			for(int i=0; i<(1000-l); i++) {
				fill += "X";
			}

			int half = (int) Math.floor(l/2);
			String start = seq.substring(0, half);
			String end = seq.substring(half, l);
					
			finalSeq = start + fill + end; 	
		} else { 
			
			String start = seq.substring(0, 500);
			String end = seq.substring(l-500, l);
			finalSeq = start + end;
		}

		return finalSeq;
	}

	private int[] getMotifPositions(HashSet<String> possibleMotifs, String sequence, int[] motifPositions){
		
		boolean motifFound = false;

		/* search for motif instances in first 500 nucleotides */
		for(int i=0; i<500-8; i++) {

			String substring = sequence.substring(i, i+8);

			/* if current substring is a motif instance, increase motif position count */
			if(possibleMotifs.contains(substring)) {
				motifFound = true;
				for(int j=i; j<i+8; j++) {
					motifPositions[j]++;
				}
			}
		}

		/* search for motif instances in last 500 nucleotides */
		for(int i=500; i<1000-8; i++) {

			String substring = sequence.substring(i, i+8);

			/* if current substring is a motif instance, increase motif position count */
			if(possibleMotifs.contains(substring)) {
				motifFound = true;
				for(int j=i; j<i+8; j++) {
					motifPositions[j]++;
				}
			}
		}
		
		if(motifFound) {
			this.motifFoundCount++;
		}
		return motifPositions;
	}


	private static void printMotifPosition(int[] motifPositions, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			for(int i=0; i<motifPositions.length; i++) {
				out.write((i+1) + "\t" + motifPositions[i] + "\n");
				out.flush();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Generate an index for the FASTA file; identifying the line number for the different refseq identifiers
	 * 
	 * @param fastaFile String - file path for the FASTA sequences
	 * @return indexOfFastaFile HashMap<String, Integer> - map of {refseqId : line count} 
	 */
	private static HashMap<String, Integer> generateFastaFileIdx(String fastaFile){

		HashMap<String, Integer> indexOfFastaFile = new HashMap<>();

		InputStream in;
		try {
			in = new FileInputStream(new File(fastaFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();
			int lineCount = 1;

			while(line != null) {
				/* store the line index of the start of a new sequence */
				if(line.startsWith(">")) {
					String[] col = line.split("_|\\.");
					String refSeqId = col[2] + "_" + col[3]; // col[2] = type of ID (eg. NM) ; col[3] = number ID (#####)

					indexOfFastaFile.put(refSeqId, lineCount);
				}

				line = input.readLine();
				lineCount++;
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return indexOfFastaFile;
	}

	private static HashMap<String, String[]> loadProteinInfoFile(String proteinInfoFile){

		HashMap<String, String[]> proteinToRefSeqIDMap = new HashMap<>();

		InputStream in;
		try {
			in = new FileInputStream(new File(proteinInfoFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();
			while(line != null) {

				String[] col = line.split("\t");
				if(col.length>1) { // some proteins in the network don't have an associated RefSeqID
					proteinToRefSeqIDMap.put(col[0], col[1].split("\\|")); // [0] = motif, [1] = list of refSeqIDs
				}
				line = input.readLine();
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return proteinToRefSeqIDMap;
	}


	private static HashSet<String> loadAnnotatedProteinsInNetwork(String proteinFreqFile){

		HashSet<String> annotatedProteinSet = new HashSet<>();

		InputStream in;
		try {
			in = new FileInputStream(new File(proteinFreqFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();
			while(line != null) {

				annotatedProteinSet.add(line.split("\t")[0]); // [0] = protein names

				line = input.readLine();
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return annotatedProteinSet;
	}

	/**
	 * Obtains full RNA seqeunce from the fasta file given the line number corresponding to the sequence header. 
	 * Will store all information following sequence header until the next sequence header 
	 * 
	 * @param fastaFile		String - file path to fasta file containing all rna sequences
	 * @param lineNumber	int - corresponds to the line number from which to start reading the file
	 * 
	 * @return sequence		String - full rna sequence 
	 */
	private String loadRNAsequence(int lineNumber) {
		String sequence = "";

		InputStream in;
		try {
			in = new FileInputStream(new File(this.fasta));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // first line

			/* read all lines until the line number specified */
			for(int i=1; i<=lineNumber; i++) {
				line = input.readLine();
			}

			/* store line info until we reach the end of the file or the next sequence header */
			while(line !=null && !line.startsWith(">")) {
				sequence += line;
				line = input.readLine();
			}
			input.close();

		}catch (IOException e) {
			e.printStackTrace();
		}
		return sequence;
	}

	@SuppressWarnings("unused")
	private ArrayList<Integer> searchForMotifPositionsInSequence(String motif, String sequence){

		ArrayList<Integer> motifPositions = new ArrayList<>();
		HashSet<String> possibleInstances = getPossibleMotifInstances(motif);


		int currentPosition = 0;

		for(int i=(sequence.length() - this.motifLength); i>=0; i--) {
			String subSeq = sequence.substring(i, i+8);
			if(possibleInstances.contains(subSeq)) {
				motifPositions.add(currentPosition);
			}

			currentPosition++;
		}

		return motifPositions;
	}

	private HashSet<String> getPossibleMotifInstances(String motif){

		HashSet<String> motifs = new HashSet<>();

		HashMap<Character, Character[]> degenCharacterMap = defineDegenCharacterMap();
		int solutions = calculateSolutions(motif, degenCharacterMap);

		/* generate all solutions */ 
		for(int i = 0; i < solutions; i++) {
			int j = 1;
			String seq = "";
			/* generate a given motif */
			for(int k=0; k < motif.length() ; k++) {
				Character[] set = degenCharacterMap.get(motif.charAt(k));

				char charToAppend = set[(i/j)%set.length]; // magic 
				seq += charToAppend;
				j *= set.length;
			}

			motifs.add(seq);

		}	
		return motifs;

	}

	/**
	 * Determine the number of possible solutions that will be generated for the given motif length
	 * This function works simply because there's the same amount of possible letters for each nucleotide
	 * Therefore the number of solutions will be the same. 
	 * 
	 * @param motifLength	int - number of nucleotides considered for a motif (motif size)
	 * @param charMap		HashMap<Character, Character[]> - map of {nucleotides : list of substitutions}
	 * 
	 * @return solutions	int - number of solutions for a given motif of length 8
	 */
	private int calculateSolutions(String motif, HashMap<Character, Character[]> degenCharacterMap) {
		int solutions = 1;
		for(int i = 0; i < motif.length(); i++) {
			solutions *= degenCharacterMap.get(motif.charAt(i)).length;
		}
		return solutions;
	} 


	/**
	 * Define the RNA nucleotide mapping to their possible substitutions (ie. degenerate characters)
	 * Note: for now this is fix - this could be passed as an input parameter in the future to enable 
	 * flexible motif representation for users
	 * Maps RNA nucleotides {A, U, C, G} to possible substitutions (themselves + degenerate characters) 
	 * 
	 * @return charMap	 HashMap<Character, Character[]> - map {nucleotide : list of possible characters}
	 */
	private HashMap<Character, Character[]> defineDegenCharacterMap(){
		HashMap<Character, Character[]> charMap = new HashMap<>();

		Character[] a = {'A'};
		Character[] t = {'T'};
		Character[] g = {'G'};
		Character[] c = {'C'};
		Character[] r = {'A', 'G'};
		Character[] y = {'T', 'C'};
		Character[] d = {'A', 'T', 'G'};
		Character[] b = {'T', 'G', 'C'};
		Character[] h = {'A', 'T', 'C'};
		Character[] v = {'A', 'G', 'C'};
		Character[] x = {'A', 'T', 'G', 'C'};

		charMap.put('A', a);
		charMap.put('T', t);
		charMap.put('G', g);
		charMap.put('C', c);
		charMap.put('R', r);
		charMap.put('Y', y);
		charMap.put('D', d);
		charMap.put('B', b);
		charMap.put('H', h);
		charMap.put('V', v);
		charMap.put('*', x);
		return charMap;
	}

	@SuppressWarnings("unused")
	private void printMotifPositions(String proteinName, String[] refSeqIds, ArrayList<ArrayList<Integer>> motifPositionByRefSeqID, String outputFile) {

		/* Decide on motif position */
		int minPosition = Integer.MAX_VALUE; // min

		/* Map of positions and their frequencies */
		HashMap<Integer, Integer> positionMode = new HashMap<>(); // mode 
		int modeMax = 1;

		HashMap<String, Integer> utrModeMap = new HashMap<>(); // mode of entire 3'utr positions
		int utrModeMax = 1;

		ArrayList<Integer> allPositions = new ArrayList<>();

		for(ArrayList<Integer> motifPositions: motifPositionByRefSeqID) {
			if(!motifPositions.isEmpty()) {

				String utrPositions = ""; // utr mode initializaion

				for(int position: motifPositions) {

					// minimum
					if(position < minPosition) {
						minPosition = position;
					}

					// mode 
					if(positionMode.containsKey(position)) {
						int currentPosition =  positionMode.get(position);
						if( currentPosition + 1 > modeMax) { 
							modeMax++;
						}
						positionMode.put(position, currentPosition +1 );
					} else { 
						positionMode.put(position, 1);
					}

					// median
					allPositions.add(position);

					utrPositions += position + ",";
				}

				// update utr mode 
				if(utrModeMap.containsKey(utrPositions)) {
					int currentUtrPositionFrq = utrModeMap.get(utrPositions);
					if(currentUtrPositionFrq + 1 > utrModeMax) {
						utrModeMax++;
					}
					utrModeMap.put(utrPositions, currentUtrPositionFrq +1);
				} else {
					utrModeMap.put(utrPositions, 1);
				}
			}

		}

		// mode
		ArrayList<Integer> modes = new ArrayList<>();
		for(int pos : positionMode.keySet()) {
			if(positionMode.get(pos) == modeMax) {
				modes.add(pos);
			}
		}

		int mode = Integer.MAX_VALUE;
		if(!modes.isEmpty()) {
			if(modes.size() > 1) {
				Collections.sort(modes);	
			}

			mode = modes.get(0);
		} 

		// utr mode 
		String utrMode = "";
		for(String utrPositions : utrModeMap.keySet()) {
			if(utrModeMap.get(utrPositions) == utrModeMax) {
				utrMode = utrPositions;
			}
		}

		// median
		Collections.sort(allPositions);
		int median = Integer.MAX_VALUE;
		if(!allPositions.isEmpty()) {
			median = allPositions.get(((int) Math.ceil(allPositions.size()/ (double) 2)) -1);	
		} 


		try {

			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile), true));

			File f = new File(outputFile);
			if (f.length() == 0) {
				out.write("ProteinName\tPositions\tPositionMedian\tPositionMode\tPositionMin\tDetails\n");
			}
			out.write(proteinName + "\t" + utrMode + "\t" + median + "\t" + mode +"\t" + minPosition + "\t");


			for(int i=0; i<refSeqIds.length; i++) {
				ArrayList<Integer> positions = motifPositionByRefSeqID.get(i);
				if(!positions.isEmpty()) {

					out.write(refSeqIds[i] + "=[");
					Collections.sort(positions);
					for(int j=0; j<positions.size(); j++){
						out.write(positions.get(j) + ",");
					}
					out.write("],");
				}
			}
			out.write("\n");

			out.flush();
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}


}
