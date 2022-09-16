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
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class DetermineMotifGenomePostions {

	public static void main(String[] args) {

		String wd = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/";
		
		String motifFile = wd + "corrNetTop2-400_coreTPD_p0.4_coreProteins_h0.7_motifFamiliesInfo.tsv";
		String fastaFile = wd + "input_files/human_3UTRsequences.txt";
		String refSeqIdFile = wd + "corrNetTop2_proteinsInNetwork_info.tsv";
		
		String positionFilePrefix = wd + "evolutionConservation/coreTPD/corrNet2-400_coreTPD_p0.4_h0.7_genomicPositions_motif";
		String bedFilePrefix = wd + "evolutionConservation/coreTPD/corrNet2-400_coreTPD_p0.4_h0.7_genomicPositions_motif";
		
		String coreProteinAnnotations = wd + "corrNetTop2-400_coreTPD_p0.4_coreProteinsByMotif.tsv";
		
		String positionCoreFilePrefix = wd + "evolutionConservation/coreTPD/corrNet2-400_coreTPD_p0.4_h0.7_genomicPositions_core_motif";
		String bedCoreFilePrefix = wd + "evolutionConservation/coreTPD/corrNet2-400_coreTPD_p0.4_h0.7_genomicPositions_core_motif";

		//getGenomePositionsForMotif(motifFile, fastaFile, refSeqIdFile, positionFilePrefix, bedFilePrefix);
		
		getGenomePositionsForMotifFromCoreProteins(motifFile, fastaFile, refSeqIdFile, coreProteinAnnotations, positionCoreFilePrefix, bedCoreFilePrefix);
	}

	private static void getGenomePositionsForMotif(String motifFile, String fastaFile, String refSeqIdFile, String positionFilePrefix, String bedFilePrefix) {

		/* get motifs to test  */
		List<String> motifsToAssess = loadMotifsToAssess(motifFile);

		/* set character map for regular expression */ 
		HashMap<Character, String> characterMap = setCharacterMapForRegularExpression();

		/* assess motifs individually */
		for(int i=0; i<motifsToAssess.size(); i++) {
			int motifNumber = i+1;
			System.out.println("Motif: " + motifNumber);
			
			String motif = motifsToAssess.get(i);

			/* convert motif to regular expression */ 
			String formattedMotif = formatMotifWithRegularExpression(motif, characterMap);
			
			String motifInfo = motif + "|" + formattedMotif;

			/* search FASTA for motif - get GENOME positions */ 
			HashSet<String> motifPositions = searchForMotifPostions(formattedMotif, fastaFile, refSeqIdFile);

			/* print all positions */
			printGenomePositionsOfMotif(motifInfo, motifNumber, motifPositions, positionFilePrefix + motifNumber + ".txt");

			/* print BED file */
			printBedFiles(motifInfo, motifNumber, motifPositions, bedFilePrefix + motifNumber + ".bed");
		}
	}
	
	private static void getGenomePositionsForMotifFromCoreProteins(String motifFile, String fastaFile, String refSeqIdFile, String coreProteinAnnotations, String positionFilePrefix, String bedFilePrefix) {

		/* get motifs to test  */
		List<String> motifsToAssess = loadMotifsToAssess(motifFile);
		
		/* set character map for regular expression */ 
		HashMap<Character, String> characterMap = setCharacterMapForRegularExpression();

		/* assess motifs individually */
		for(int i=0; i<motifsToAssess.size(); i++) {
			int motifNumber = i+1;
			System.out.println("Motif: " + motifNumber);
			
			String motif = motifsToAssess.get(i);
			
			/* load core proteins */ 
			HashSet<String> coreProteinSet = loadCoreProteins(motif, coreProteinAnnotations);

			/* convert motif to regular expression */ 
			String formattedMotif = formatMotifWithRegularExpression(motif, characterMap);
			
			String motifInfo = motif + "|" + formattedMotif;

			/* search FASTA for motif - get GENOME positions */ 
			HashSet<String> motifPositions = searchForMotifPostionsFromCore(formattedMotif, fastaFile, refSeqIdFile, coreProteinSet);

			/* print all positions */
			printGenomePositionsOfMotif(motifInfo, motifNumber, motifPositions, positionFilePrefix + motifNumber + ".txt");

			/* print BED file */
			printBedFiles(motifInfo, motifNumber, motifPositions, bedFilePrefix + motifNumber + ".bed");
		}
	}

	/**
	 * Load motifs from a .TSV containing the motif name in the 1st column
	 * 
	 * @param motifFile			String - File.tsv
	 * @return motifsToAssess	List<String> - list of motifs in order 
	 */
	private static List<String> loadMotifsToAssess(String motifFile){

		List<String> motifsToAssess = new ArrayList<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(motifFile))));			

			String line = in.readLine(); // header
			line = in.readLine(); 
			
			while(line!=null) {
				motifsToAssess.add(line.split("\t")[0]); // [0] motif 
				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return motifsToAssess;
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
	
	
	private static HashSet<String> searchForMotifPostions(String formattedMotif, String fastaFile, String refSeqIdFile){

		HashSet<String> motifPositions = new HashSet<>();
		
		/* for each gene select (1) RefSeqId to test */
		HashSet<String> refSeqIds = getRefSeqIds(refSeqIdFile);

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

	private static HashSet<String> getRefSeqIds(String refSeqIdFile){

		HashSet<String> refSeqIds = new HashSet<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(refSeqIdFile))));			

			String line = in.readLine();

			while(line!=null) {
				
				if(line.split("\t").length > 1) {
					String[] ids = line.split("\t")[1].split("\\|"); // [0] protein name, [1] list of refSeqIds
					
					int numMotifs = ids.length;
					Random rand = new Random();
					
					int randomIdx = rand.nextInt(numMotifs); // random number between 0 and numMotifs (excluding upper bound)
					refSeqIds.add(ids[randomIdx]);	
				}

				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return refSeqIds;	
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

	private static void printGenomePositionsOfMotif(String motif, int number, HashSet<String> motifPositions, String positionsFile) {
		
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(positionsFile)));
			
			out.write("Motif " + number + " : " + motif + " | Number of Occurences : " + motifPositions.size() + "\n");
			
			for(String pos: motifPositions) {
				out.write(pos + "\n");
				out.flush();
			}
			
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	private static void printBedFiles(String motif, int number, HashSet<String> motifPositions, String bedFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(bedFile)));
			
			out.write("Track - Motif " + number + " : " + motif + "\n");
			
			for(String pos: motifPositions) {
				
				String[] split = pos.split(":");
				out.write(split[0] +"\t"+ (Integer.parseInt(split[1].split("-")[0])+3) +"\t"+(Integer.parseInt(split[1].split("-")[0])+4) + "\n");
				out.flush();
			}
			
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	
	private static HashSet<String> loadCoreProteins(String motif, String coreProteinAnnotationFile){
		
		HashSet<String> coreProteins = new HashSet<>();
		
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(coreProteinAnnotationFile))));			

			String line = in.readLine();

			while(line!=null) {

				String currentMotif = line.split("\t")[0];
				if(currentMotif.equals(motif)) {
					coreProteins.addAll(Arrays.asList(line.split("\t")[2].split("\\|")));
					break;
				}

				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return coreProteins;
	}

	private static HashSet<String> searchForMotifPostionsFromCore(String formattedMotif, String fastaFile, String refSeqIdFile, HashSet<String> coreProteins){

		HashSet<String> motifPositions = new HashSet<>();
		
		/* for each gene select (1) RefSeqId to test */
		HashSet<String> refSeqIds = getCoreRefSeqIds(refSeqIdFile, coreProteins);

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
	
	private static HashSet<String> getCoreRefSeqIds(String refSeqIdFile, HashSet<String> coreProteins){

		HashSet<String> refSeqIds = new HashSet<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(refSeqIdFile))));			

			String line = in.readLine();

			while(line!=null) {
				
				String[] col = line.split("\t");
				
				if(col.length > 1) {
					
					if(coreProteins.contains(col[0])){ // [0] protein name
						
						String[] ids = col[1].split("\\|"); // [0] protein name, [1] list of refSeqIds
						
						int numMotifs = ids.length;
						Random rand = new Random();
						
						int randomIdx = rand.nextInt(numMotifs); // random number between 0 and numMotifs (excluding upper bound)
						refSeqIds.add(ids[randomIdx]);	
					}
				}

				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return refSeqIds;	
	}
}
