package sequenceLength;

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


public class DetermineSequenceLength {

	public static void main(String[] args) {
		
		String wd = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/";
		String fastaFile = wd + "/input_files/human_3UTRsequences.txt";
		
		String refSeqIdFile = wd + "corrNetTop2_proteinsInNetwork_info.tsv";
		String outptutFile = wd + "nucleotide-composition/corrNetTop2_sequenceLengths_allSeqsConsidered.tsv";
		
		System.out.println("* Assess - all sequences *");
		determine3UTRsequenceLength(refSeqIdFile, fastaFile, outptutFile);
		
		String refSeqIdLongestSeqsFile = wd + "/MotifPosition/corrNetTop2_longestSequence.tsv";
		String outputLongestSeqsFile = wd + "nucleotide-composition/corrNetTop2_sequenceLengths_LongestSeqs.tsv";
		
		System.out.println("* Assess - longest sequences *");
		determine3UTRsequenceLength(refSeqIdLongestSeqsFile, fastaFile, outputLongestSeqsFile);
		
		String annotationFile = wd + "corrNetTop2-400_coreTPD_p0.4_coreProteinsByMotif.tsv";
		String outputFamilyPrefix = wd + "nucleotide-composition/families/sequenceLenght_allSeqs_coreProteins_motifFamily";
		String outputAll = wd + "nucleotide-composition/families/sequenceLenght_motifFamilies.tsv";
		String motifFamily = wd + "corrNetTop2-400_coreTPD_p0.4_coreProteins_h0.7_motifFamiliesInfo.tsv";
		
		System.out.println("* Assess - sequence lengths of motif families *");
		determineSequenceLengthsForMotifFamilies(motifFamily, annotationFile, refSeqIdFile, fastaFile, outputFamilyPrefix, outputAll);
	}

	private static void determine3UTRsequenceLength(String refSeqIdsFile, String fastaFile, String outputFile) {

		/* load RefseqIds */
		HashSet<String> idSet = loadRefSeqIds(refSeqIdsFile);

		/* determine nucleotide composition from FASTA file */
		List<Integer> seqLengths = countSequenceLength(fastaFile, idSet);
		
		/* output composition to file */
		printCounts(seqLengths, outputFile);
		
	}
	
	private static void determineSequenceLengthsForMotifFamilies(String motifFamilies, String annotationFile, String refSeqIdsFile, String fastaFile, String outputFilePrefix, String outputAll) {
		
		/* Load motif families */
		List<Motif> motifList = loadMotifFamilies(motifFamilies);
		
		
		/*load all refseqIDs */ 
		HashMap<String, String[]> idSets = loadRefSeqIdsPerProtein(refSeqIdsFile);
		
		for(Motif m : motifList) {
			
			System.out.println(m.getMotifName());
			/* For each family determine list of proteins */ 
			m.setProteins(loadProteinsPerMotif(m.getMotifName(), annotationFile));
			
			/* For each protein determine refseqIDs */
			m.setRefSeqIds(idSets);
			
			/* count sequence lengths and print */
			m.setSequenceLenghts(countSequenceLength(fastaFile, m.getIds()));
			
			printCounts(m.getLengths(), outputFilePrefix + m.getMotifOrder() + ".tsv");
		}		
		combineCount(motifList, outputAll);
	}

	private static HashSet<String> loadRefSeqIds(String refSeqIdFile){

		HashSet<String> idSet = new HashSet<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(refSeqIdFile))));

			String line = in.readLine();

			while(line != null) {
				
				if(line.split("\t").length > 1 ) {
					String[] ids = line.split("\t")[1].split("\\|"); // [1] = list of IDs {id1|id2|..|idx} 
					idSet.addAll(Arrays.asList(ids));
				}

				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return idSet;
	}
	
	
	private static List<Integer> countSequenceLength(String fastaFile, HashSet<String> ids){
		
		List<Integer> seqLengths = new ArrayList<>();
		
		/* read through FASTA file */
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(fastaFile))));			

			String line = in.readLine();

			boolean readSeq = false; 
			String seq = "";

			while(line!=null) {

				if(readSeq & !line.startsWith(">")) {
					seq += line;
				}

				/* check if sequence is associated to a protein in network */
				if(line.startsWith(">")) {
					readSeq = false;

					/* if sequence has been loaded - assess nucleotide composition  */
					if(!seq.isEmpty()) {
						seqLengths.add(seq.length());
						seq = "";
					}

					//  >hg38_ncbiRefSeq_NM_001276352.2 range=chr1:67092165-67093579 5'pad=0 3'pad=0 strand=- repeatMasking=none					
					String id = line.split("[\\_\\s++\\.]")[2] + "_"+ line.split("[\\_\\s++\\.]")[3];

					if(ids.contains(id)) {
						readSeq = true;
						seq = "";
					}
				}
				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return seqLengths;
	}
	
	private static void printCounts(List<Integer> seqCounts, String outputFile) {
		
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			
			for(int i=0; i<seqCounts.size(); i++) {
				out.write(seqCounts.get(i) + "\n");
			}
			
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private static ArrayList<Motif> loadMotifFamilies(String motifFamilies){
		
		ArrayList<Motif> motifList = new ArrayList<>();
		
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(motifFamilies))));

			String line = in.readLine(); // header
			line = in.readLine();

			while(line != null) {
				
				String[] col = line.split("\t");
				motifList.add(new Motif(col[0], Integer.parseInt(col[1])));

				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return motifList;
	}
	

	private static HashMap<String, String[]> loadRefSeqIdsPerProtein(String refSeqIdFile){

		HashMap<String, String[]> idSet = new HashMap<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(refSeqIdFile))));

			String line = in.readLine();

			while(line != null) {
				
				if(line.split("\t").length > 1 ) {
					String[] ids = line.split("\t")[1].split("\\|"); // [1] = list of IDs {id1|id2|..|idx} 
					idSet.put(line.split("\t")[0], ids);
				}

				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return idSet;
	}
	
	private static HashSet<String> loadProteinsPerMotif(String motif, String annotationFile){
		HashSet<String> proteinSet = new HashSet<>();
		
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(annotationFile))));

			String line = in.readLine();

			while(line != null) {
				
				if(line.split("\t")[0].equals(motif)) {
					proteinSet.addAll(Arrays.asList(line.split("\t")[2].split("\\|")));
					break;
				}
				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return proteinSet;
	}
	
	private static void combineCount(List<Motif> motifList, String outputFile) {
		
		/* determine max count */
		int max = 0;
		for(Motif m: motifList) {
			if(m.getLengths().size() > max) {
				max = m.getLengths().size();
			}
		}
		
		/* print combined data */
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			
			/* header */
			for(int i=0; i<motifList.size(); i++) {
				out.write("F" + (i+1) + "\t");
			}
			out.write("\n");
			
			/* body */
			for(int j=0; j<max; j++) {
				for(int i=0; i<motifList.size(); i++) {
					
					if(motifList.get(i).getLengths().size() > j) {
						out.write(motifList.get(i).getLengths().get(j) + "\t");
					} else { 
						out.write("NA\t");
					}
				}
				out.write("\n");
			}
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
}
