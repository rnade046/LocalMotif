package evolutionConservation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class DetermineGenomePositionsOfMotifsTestedInNetwork {

	public static void main(String[] args) {
	
		String proteinInfoFile = "corrNetTop2-400_proteinsInNetwork_info.tsv";
		String companionFile = "/home/rnade046/scratch/annotationFiles_Full_FWD/CompanionFiles_corrNetTop2-400_n20_2000/corrNetTop2-400_n20_2000_motifAnnotationsCompanionFile_"+ args[0];
		String annotationFile = "/home/rnade046/scratch/annotationFiles_Full_FWD/corrNet_degenAnnotations_FWD_"+ args[0] +".tsv";
		String fastaFile = "human_3UTRsequences.txt";
		String outputFile = "positions/genomePositions_" + args[0];
		
	
		/* search for motifs in given annotation/companion file */
		assessMotifs(proteinInfoFile, companionFile, annotationFile, fastaFile, outputFile);

	}



	public static void assessMotifs(String proteinInfo, String companionFile, String annotationFile, String fastaFile, String outputFile) {

		/* load refseqIds associated to motif */
		HashMap<String, HashSet<String>> refSeqIds = loadProteinRefSeqIdMap(proteinInfo);
		
		/* load motifs to test */
		HashSet<String> testedMotifs = SequenceUtils.loadRepresentativeMotifs(companionFile);

		/* set character map */
		HashMap<Character, String> characterMap = setCharacterMapForRegularExpression();

		HashSet<String> allMotifPositions = new HashSet<>();

		/* search through annotation file; for motifs */
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(annotationFile))));

			String line = in.readLine();
			int lineCount = 1;
			while(line!=null) {
				System.out.print(lineCount + ".");
				if(lineCount%50 ==0) {
					System.out.println();
				}
				lineCount++;
				
				String motif = line.split("\t")[0];
				if(testedMotifs.contains(motif)) { // [0] = motifs

					/* obtain refSeqIds known to contain motif */
					HashSet<String> allRefSeqIds = new HashSet<>();

					String[] proteins = line.split("\t")[2].split("\\|"); // [2] = list of proteins
					for(String p: proteins) {
						if(refSeqIds.containsKey(p)) {
							allRefSeqIds.addAll(refSeqIds.get(p));
						}
					}

					/* determine possible instances of motif */
					String formattedMotif = formatMotifWithRegularExpression(motif, characterMap);

					/* search for motifs in known refSeqIds */
					allMotifPositions.addAll(searchForMotifPostions(formattedMotif, fastaFile, allRefSeqIds));
					
					//printMotifPositions(allMotifPositions, outputFile);
				}
				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		/* combine motif positions and print */ 
		printMotifPositions(allMotifPositions, outputFile);
	}

	/**
	 * Load RefSeq IDs associated to motif
	 * @param proteinInfo
	 * @return proteinMap	HashMap<String, HashSet<String>> - protein = id1|2|3..|4
	 */
	public static HashMap<String, HashSet<String>> loadProteinRefSeqIdMap(String proteinInfo){

		HashMap<String, HashSet<String>> proteinMap = new HashMap<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(proteinInfo))));

			String line = in.readLine();
			while(line!=null) {

				String[] col = line.split("\t"); // [0] = protein, [1] = refSeqIds

				if(col.length > 1) {
					proteinMap.put(col[0], new HashSet<String>(Arrays.asList(col[1].split("\\|"))));
				}
				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return proteinMap;
	}
	private static HashSet<String> searchForMotifPostions(String formattedMotif, String fastaFile, HashSet<String> refSeqIds){

		HashSet<String> motifPositions = new HashSet<>();
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(fastaFile))));			

			String line = in.readLine();

			boolean readSeq = false; 
			String seq = "";
			String header = "";

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
						motifPositions.addAll(positions);
					}

					//  >hg38_ncbiRefSeq_NM_001276352.2 range=chr1:67092165-67093579 5'pad=0 3'pad=0 strand=- repeatMasking=none					
					String id = line.split("[\\_\\s++\\.]")[2] + "_"+ line.split("[\\_\\s++\\.]")[3];

					if(refSeqIds.contains(id)) { // do not consider alternate chromosomes

						if(!line.contains("alt") && !line.contains("_fix")) {
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
	
	

	private static List<String> searchSeqForMotif(String motif, String sequence, String header) {

		List<String> positions = new ArrayList<>();

		Pattern pattern = Pattern.compile(motif);		// compile motif as REGEX
		Matcher matcher = pattern.matcher(sequence);	// match pattern to sequence 

		String chromosome = header.split(" ")[1].split("=")[1].split(":")[0];
		chromosome = chromosome+":";
		int startPos =  Integer.parseInt(header.split(" ")[1].split("=")[1].split(":")[1].split("-")[0]);
		int endPos =  Integer.parseInt(header.split(" ")[1].split("=")[1].split(":")[1].split("-")[1]);

		/* check for all instances of motif in sequence */
		while(matcher.find()) {

			/* obtain start position of motif */
			int pos =matcher.start();

			String genomicPosition = "";

			if(!header.contains("+")){

				int motifEnd = endPos-pos;
				int motifStart = motifEnd - 7;

				genomicPosition = chromosome+motifStart+"-"+motifEnd;
			}
			else{
				int motifStart = startPos+pos;
				int motifEnd = motifStart+7;
				genomicPosition = chromosome+motifStart+"-"+motifEnd;
			}

			positions.add(genomicPosition); 
		}
		return positions;
	}

	/**
	 * Set character map - hard coded with alphabet used throughout analysis 
	 * 
	 * Values follow known IUPAC regular expression e.g. R = "[AG]" 
	 * @return characterMap		HashMap<Character, String> 
	 */
	private static HashMap<Character, String> setCharacterMapForRegularExpression(){

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
	private static String formatMotifWithRegularExpression(String motif, HashMap<Character, String> characterMap) {

		String formattedMotif = "";

		for(int i = 0; i < motif.length(); i++){
			formattedMotif += characterMap.get(motif.charAt(i));
		}

		return formattedMotif;
	}

	public static void printMotifPositions(HashSet<String> motifPositions, String outputFile) {

		/* combine motif positions */
		HashMap<String, HashSet<Integer>> combinedPositions = new HashMap<>(); // CHR = Position 1, 2, 3

		for(String pos : motifPositions) {

			String[] p = pos.split("\\:"); // [0] = chr; [1] = position
			if(combinedPositions.containsKey(p[0])){

				HashSet<Integer> currentPositions = combinedPositions.get(p[0]);
				String[] indexes = p[1].split("\\-");

				for(int i=Integer.parseInt(indexes[0]); i<=Integer.parseInt(indexes[1]); i++) {
					currentPositions.add(i);
				}

				combinedPositions.put(p[0], currentPositions);
			} else {

				HashSet<Integer> positions = new HashSet<>();
				String[] indexes = p[1].split("\\-");

				for(int i=Integer.parseInt(indexes[0]); i<=Integer.parseInt(indexes[1]); i++) {
					positions.add(i);
				}

				combinedPositions.put(p[0], positions);
			}
		}

		/* print and sort positions; 1 line per chromosomes*/
		BufferedWriter out;
		try {
			out = new BufferedWriter(new FileWriter(new File(outputFile)));

			for(Entry<String, HashSet<Integer>> e: combinedPositions.entrySet()) {

				out.write(e.getKey() + "\t");

				List<Integer> order = new ArrayList<Integer>(e.getValue());
				Collections.sort(order);

				int startIndex = order.get(0); // start index 
				int idx1 = startIndex;
				int idx2 = 0; 
				for(int i=1; i<order.size(); i++) {

					idx2 = order.get(i);

					if(idx2 != (idx1+1)) {
						out.write(startIndex + "-" + idx1 + ",");
						out.flush();
						
						startIndex = idx2;
					}
					
					idx1 = idx2;
				}

				out.write("\n");
				out.flush();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

}
