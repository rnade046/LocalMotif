package positionConservation;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

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

	public static List<String> loadRepresentativeMotifs(String inputFile){

		List<String> motifs = new ArrayList<>();

		InputStream in;
		try {
			in = new FileInputStream(new File(inputFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // header
			line = input.readLine();

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
	
}
