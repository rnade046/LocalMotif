package motifInfo;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class RegexSearch {

	public static HashMap<String, List<String>> searchForMotifPostions(String formattedMotif, String fastaFile, HashSet<String> refSeqIds){

		HashMap<String, List<String>> motifPositions = new HashMap<>();
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(fastaFile))));			

			String line = in.readLine();

			boolean readSeq = false; 
			String seq = "";
			String header = "";
			String refSeqId = "";

			while(line!=null) {

				if(readSeq) {
					seq += line;
				}

				/* check if sequence is associated to a protein in network */
				if(line.startsWith(">")) {
					readSeq = false;

					/* if sequence has been loaded - check for motif in the sequence  */
					if(!seq.isEmpty()) {
						List<String> positions = searchSeqForMotif(formattedMotif, seq, header);
						motifPositions.put(refSeqId, positions);
					}

					//  >hg38_ncbiRefSeq_NM_001276352.2 range=chr1:67092165-67093579 5'pad=0 3'pad=0 strand=- repeatMasking=none					
					String id = line.split("[\\_\\s++\\.]")[2] + "_"+ line.split("[\\_\\s++\\.]")[3];

					if(refSeqIds.contains(id)) { // do not consider alternate chromosomes
						if(!line.contains("alt") && !line.contains("_fix")) {
							refSeqId = id;
							readSeq = true;
							header = line;
							seq = "";
						}
					}
					
				}

				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return motifPositions;
	}

	
	public static List<String> searchSeqForMotif(String motif, String sequence, String header) {

		List<String> positions = new ArrayList<>();

		Pattern pattern = Pattern.compile(motif);		// compile motif as REGEX
		Matcher matcher = pattern.matcher(sequence);	// match pattern to sequence 

		/* check for all instances of motif in sequence */
		while(matcher.find()) {

			/* obtain start position of motif */
			int pos =matcher.start();
			positions.add(Integer.toString(pos));
		}
		return positions;
	}

	/**
	 * Set character map - hard coded with alphabet used throughout analysis 
	 * 
	 * Values follow known IUPAC regular expression e.g. R = "[AG]" 
	 * @return characterMap		HashMap<Character, String> 
	 */
	public static HashMap<Character, String> setCharacterMapForRegularExpression(){

		HashMap<Character, String> characterMap = new HashMap<>();

		characterMap.put('A', "A");
		characterMap.put('C', "C");
		characterMap.put('G', "G");
		characterMap.put('T', "T");
		characterMap.put('R', "[AG]");
		characterMap.put('Y', "[CT]");
		characterMap.put('D', "[ATG]");
		characterMap.put('B', "[TGC]");
		characterMap.put('H', "[AUC]");
		characterMap.put('V', "[AGC]");
		characterMap.put('*', ".");
		return characterMap;
	}

	/**
	 * Take input motif and translate to corresponding regular expression. 
	 * 
	 * @param motif				String - initial motif
	 * @param characterMap		HashMap<Character, String> - conversion map with regular expressions
	 * @return formattedMotif	String - motif formatted with regular expression
	 */
	public static String formatMotifWithRegularExpression(String motif, HashMap<Character, String> characterMap) {

		String formattedMotif = "";

		for(int i = 0; i < motif.length(); i++){
			formattedMotif += characterMap.get(motif.charAt(i));
		}

		return formattedMotif;
	}
}
