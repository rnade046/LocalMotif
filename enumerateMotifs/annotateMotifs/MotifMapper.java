package annotateMotifs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class MotifMapper {

	private HashMap<Character, String> characterMap;
	private HashMap<String, String> idToProteinMap;
	private String annotationFilePrefix;
	private File motifDirectory;
	private String seqFasta;
	private String motifFilePrefix;

	public MotifMapper(String ids, String a, String fasta, String dir, String motifs) {

		this.characterMap = setCharacterMapForRegularExpression();
		this.idToProteinMap = loadProteinRefSeqIdMap(ids);

		this.annotationFilePrefix = a;
		this.motifDirectory = new File(dir);
		this.seqFasta = fasta;
		this.motifFilePrefix = motifs;
	}

	public void mapMotifsToTheirAssociatedProteins() {

		/* iterate over list of motifs */
		int fileCount = this.motifDirectory.list().length;

		for(int i=0; i<fileCount; i++) {

			System.out.println("searching for motifs in file: " + i);

			/* load motifs to assess and determine their regular expression */
			ArrayList<Motif> motifs = initializeMotifs(motifFilePrefix + i);

			Set<String> refSeqIds = this.idToProteinMap.keySet();

			try {
				BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(this.seqFasta))));			

				String line = in.readLine();

				boolean readSeq = false; 
				String seq = "";
				String id = "";
				int seqCount = 1;

				while(line!=null) {

					if(readSeq) {
						seq += line;
					}

					/* check if sequence is associated to a protein in network */
					if(line.startsWith(">")) {
						readSeq = false;
						
						/* if sequence has been loaded - check for motif in the sequence  */
						if(!seq.isEmpty() && refSeqIds.contains(id)) {
							
							//System.out.println("sequence : " + id);

							if(seqCount%10 == 0) {
								System.out.print(seqCount + ".");
							}
							
							if(seqCount%100 == 0) {
								System.out.println();
							}
							
							for(Motif m : motifs) {
								
								boolean motifFound = searchSeqForMotif(m.getRegexMotif(), seq);

								if(motifFound) {
									m.addProtein(this.idToProteinMap.get(id));
								}
							}
						}
						//  >hg38_ncbiRefSeq_NM_001276352.2 range=chr1:67092165-67093579 5'pad=0 3'pad=0 strand=- repeatMasking=none					
						id = line.split("[\\_\\s++\\.]")[2] + "_"+ line.split("[\\_\\s++\\.]")[3];

						/* reinitialize parameters for next sequence */
						if(refSeqIds.contains(id)) { 

							if(!line.contains("alt") && !line.contains("_fix")) { // do not consider alternate chromosomes
								readSeq = true;
								seq = "";
								seqCount++;
							}
						}
					}
					line = in.readLine();
				}
				in.close();
			} catch (IOException e) {
				e.printStackTrace();
			}

			printAnnotationFile(this.annotationFilePrefix + i, motifs);
		}
	}

	private ArrayList<Motif> initializeMotifs(String motifFile){

		ArrayList<Motif> motifList = new ArrayList<>();

		try {
			FileInputStream in = new FileInputStream(new File(motifFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String motif = input.readLine();
			while(motif !=null) {

				/* determine possible instances of motif */
				String regexMotif = formatMotifWithRegularExpression(motif);
				motifList.add(new Motif(motif, regexMotif));

				/* read next motif in file */
				motif = input.readLine();
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return motifList;
	}

	/**
	 * Print annotation file { motif -- #proteins -- protein1|protein2|..| } 
	 * 
	 * @param annotationFile	String - annotation file path
	 * @param motif				String - UPAC motif	
	 * @param proteins			HashSet<String> - proteins associated to motif
	 */
	private void printAnnotationFile(String annotationFile, ArrayList<Motif> motifs) {

		BufferedWriter out;
		try {
			out = new BufferedWriter(new FileWriter(new File(annotationFile), true));

			for(Motif m: motifs) {
				out.write(m + "\t" + m.getProteins().size() + "\t");

				for(String p : m.getProteins()) {
					out.write(p + "|");
					out.flush();
				}
				out.write("\n");
				out.flush();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Take input motif and translate to corresponding regular expression. 
	 * 
	 * @param motif				String - initial motif
	 * @param characterMap		HashMap<Character, String> - conversion map with regular expressions
	 * @return formattedMotif	String - motif formatted with regular expression
	 */
	private String formatMotifWithRegularExpression(String motif) {

		String formattedMotif = "";

		for(int i = 0; i < motif.length(); i++){
			formattedMotif += this.characterMap.get(motif.charAt(i));
		}

		return formattedMotif;
	}

	/**
	 * Search for regex motif in the provided sequence 
	 * 
	 * @param motif			String - motif formatted with regular expression
	 * @param sequence		String - RNA sequence to search
	 * @return motifFound	Boolean - true if motif is found in sequence
	 */
	private static boolean searchSeqForMotif(String motif, String sequence) {

		boolean motifFound = false; 

		Pattern pattern = Pattern.compile(motif);		// compile motif as REGEX
		Matcher matcher = pattern.matcher(sequence);	// match pattern to sequence 

		/* check if motif is contained in the sequence */
		if(matcher.find()) {
			motifFound = true;
		}
		return motifFound;
	}

	/**
	 * Set character map - hard coded with alphabet used throughout analysis 
	 * 
	 * Values follow known IUPAC regular expression e.g. R = "[AG]" 
	 * @return characterMap		HashMap<Character, String> 
	 */
	private HashMap<Character, String> setCharacterMapForRegularExpression(){

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

	private HashMap<String, String> loadProteinRefSeqIdMap(String proteinInfo){

		HashMap<String, String> proteinIdMap = new HashMap<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(proteinInfo))));

			String line = in.readLine();
			while(line!=null) {

				String[] col = line.split("\t"); // [0] = protein, [1] = refSeqIds

				if(col.length > 1) {
					for(String id: col[1].split("\\|")) {
						proteinIdMap.put(id, col[0]);
					}
				}
				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return proteinIdMap;
	}
}
