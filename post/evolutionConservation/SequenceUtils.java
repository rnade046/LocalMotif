package evolutionConservation;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.HashSet;

public class SequenceUtils {

	public static String formatSequence(String seq, int testLength) {
		String finalSeq = "";

		int l = seq.length();

		if(l == testLength) {
			finalSeq = seq;

		} else if(l < testLength) {

			String fill = "";
			for(int i=0; i<(testLength-l); i++) {
				fill += "X";
			}

			int half = (int) Math.floor(l/2);
			String start = seq.substring(0, half);
			String end = seq.substring(half, l);

			finalSeq = start + fill + end;
		} else {

			int half = (int) Math.floor(testLength/2);
			String start = seq.substring(0, half);
			String end = seq.substring(l-half, l);
			finalSeq = start + end;
		}

		return finalSeq;
	}


	public static HashSet<String> getPossibleMotifInstances(String motif){

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
	private static int calculateSolutions(String motif, HashMap<Character, Character[]> degenCharacterMap) {
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
	private static HashMap<Character, Character[]> defineDegenCharacterMap(){
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

	public static HashSet<String> loadRepresentativeMotifs(String inputFile){

		HashSet<String> motifs = new HashSet<>();

		InputStream in;
		try {
			in = new FileInputStream(new File(inputFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // no header

			while(line != null) {
				motifs.add(line.split("\t")[0]);

				line = input.readLine();
			}

			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return motifs;
	}

	public static String getFastaForId(String refSeqId, String fastaFile) {

		String sequence = "";

		/* search for FASTA */
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader (new FileInputStream(new File(fastaFile))));
			String line = in.readLine();

			boolean storeSeq = false;
			String seq = "";
			String id = "";

			while(line!=null) {

				/* store sequence (add multiple lines) */
				if(storeSeq) {
					seq += line;
				}

				if(line.startsWith(">")) {

					/* store length */
					if(!seq.isEmpty()) {
						sequence = seq;

						break;
					}

					storeSeq = false;
					seq = "";

					String[] col = line.split("_|\\.");
					id = col[2] + "_" + col[3]; // col[2] = type of ID (eg. NM) ; col[3] = number ID (#####)

					if(refSeqId.equals(id)) {
						storeSeq = true;
					}
				}
				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return sequence;
	}

	public static String[] formatSequenceInBins(String sequence, int bins) {

		String[] sequenceBins = new String[bins];

		int[] binSizes = new int[bins];
		int[] binSizesReodered = new int[bins];

		int minimumBinSize = (int) Math.floor(sequence.length() / (double) bins); 
		int binRemainder = sequence.length() % bins; 

		/* determine size of bins (adding remainder to minimum bin sizes) */
		for(int i=0; i<bins; i++) {

			// determining indexes for unequal bins : https://stackoverflow.com/questions/2130016/splitting-a-list-into-n-parts-of-approximately-equal-length
			int binStart = i*minimumBinSize + Math.min(i, binRemainder);
			int binEnd = (i+1) * minimumBinSize + Math.min(i+1, binRemainder);
			int binRange = binEnd - binStart;

			binSizes[i] = binRange;
		}

		/* re-ordering bin sizes to pad edges */
		int countEven = 0;
		int countOdd = bins-1;
		for(int i=0; i<bins; i++) {
			if(i%2 ==0) {
				binSizesReodered[countEven] = binSizes[i];
				countEven++;
			} else { 
				binSizesReodered[countOdd] = binSizes[i];
				countOdd--;
			}
		}

		/* get sequence substring */
		int currentIdx = 0;
		for(int i=0; i<bins; i++) {

			int start = currentIdx;
			int end = start + binSizesReodered[i];
			currentIdx += binSizesReodered[i];

			sequenceBins[i] = sequence.substring(start, end);
		}
		return sequenceBins;

	}
}
