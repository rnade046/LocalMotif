package ClusteredMotifs;

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

	// Constructor
	public PositionConservation(String fastaFile, String proteinInfoFile, int l) {

		this.fasta = fastaFile;
		this.motifLength = l;
		System.out.println("Generate fasta file index");
		this.fastaIdx = generateFastaFileIdx(fastaFile); // generate index for FASTA file
		System.out.println("Load protein info file");
		this.proteinInfoMap = loadProteinInfoFile(proteinInfoFile); // load protein info file
	}

	public void getMotifPositions(String extractedAnnotationsFile, String motifOutputPrefixFile) {

		/* Load significant motif and its annotated proteins from extracted annotation 1 at a time */ 
		InputStream in;
		try {
			in = new FileInputStream(new File(extractedAnnotationsFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();

			while(line != null) {

				String[] col = line.split("\t"); // [0] = motif, [1] = proteinList
				int motifCount = 1;
				/* For every protein get list of refseqIDs*/
				for(String prot: col[1].split("\\|")) {

					if(this.proteinInfoMap.containsKey(prot)) {

						String[] refSeqIDs = this.proteinInfoMap.get(prot);
						ArrayList<ArrayList<Integer>> listMotifPositionsByRefSeqId = new ArrayList<>();

						for(String id: refSeqIDs) {
							if(this.fastaIdx.containsKey(id)) {
								int indexOfFasta = this.fastaIdx.get(id);

								/* Load sequence */
								String sequence = loadRNAsequence(indexOfFasta);

								/* Search for motif in sequence */
								listMotifPositionsByRefSeqId.add(searchForMotifPositionsInSequence(col[0], sequence)); // col[0] = motifs
							} else {
								listMotifPositionsByRefSeqId.add(new ArrayList<Integer>());
							}
						}

						/* Print results 1 file per motif; protein \t position \t RefSeqId1=[Position,Position,etc.], RefSeqID2=[Position, Position]  */
						/* Protein position is decided as the median, mode or smallest possible or random? */ 
						String motifOutputFile = motifOutputPrefixFile + col[0] + "_motif" + motifCount;
						printMotifPositions(prot, refSeqIDs, listMotifPositionsByRefSeqId, motifOutputFile);
					}
				}

				line = input.readLine();
				motifCount++;
			}
			input.close();
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
				proteinToRefSeqIDMap.put(col[0], col[1].split("\\|")); // [0] = motif, [1] = list of refSeqIDs
				line = input.readLine();
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return proteinToRefSeqIDMap;
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

	private void printMotifPositions(String proteinName, String[] refSeqIds, ArrayList<ArrayList<Integer>> motifPositionByRefSeqID, String outputFile) {

		/* Decide on motif position */
		int minPosition = Integer.MAX_VALUE; // min

		HashMap<Integer, Integer> positionMode = new HashMap<>(); // mode 
		int modeMax = 0;

		ArrayList<Integer> allPositions = new ArrayList<>();

		for(ArrayList<Integer> motifPositions: motifPositionByRefSeqID) {
			if(!motifPositions.isEmpty()) {
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
		if(modes.size() > 1) {
			Collections.sort(modes);	
		}
		mode = modes.get(0);

		// median
		Collections.sort(allPositions);
		int median = allPositions.get((int) Math.ceil(allPositions.size()/ (double) 2));

		try {
			File f = new File(outputFile);
			BufferedWriter out = new BufferedWriter(new FileWriter(f, true));
			
			if (f.length() == 0) {
				out.write("ProteinName\tPositionMedian\tPositionMode\tPositionMin\tDetails");
			}
			out.write(proteinName + "\t" + median + "\t" + mode +"\t" + minPosition + "\t");
			
			
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

			out.flush();
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

}
