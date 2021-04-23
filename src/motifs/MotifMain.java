package motifs;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.HashSet;
import java.util.Properties;

public class MotifMain {

	public static void main(String[] args) throws IOException{

		/**
		 * This program takes a list of protein names as HGNC symbols output from tested network from LESMONlocal, 
		 * generates the list of non degenerate and degenerate motifs from the FASTA sequence of these proteins 3'UTR mRNAs,
		 * and formats the output as a mapping of motifs to list of proteins that will be used by LESMoNlocal to perform
		 * the motif enrichment analysis. 
		 * 
		 * Note to self: eventually should figure out how to run all of these in one run (e.g. have bash script call other scripts)
		 * 
		 * Steps: 
		 * 1- In R, use biomaRt package to obtain the mapping of HGNC symbols to RefSeqIds		>> To integrate
		 * 2- Enumerate motifs associated to all RefSeqIds from proteins associated 3'UTR mRNA sequence 
		 * 3- Enumerate possible degenerate motifs from list of motifs 			>> Run in parallel on CC
		 * 4- Map motifs and degenerate motifs to their associated proteins 	>> Run in parallel on CC
		 *
		 * Required input: 
		 * 	- List of proteins in the network (as HGNC symbols)
		 * 	- 3'UTR mRNA sequences of human genes (FASTA format)
		 */ 
		
		System.out.println("Loading parameters file \n");
		Properties params = new Properties();
		params.load(new FileInputStream(args[0]));
		
		boolean enumerateMotifs = Boolean.parseBoolean(params.getProperty("enumerateMotifs"));
		boolean enumerateDegenMotifs = Boolean.parseBoolean(params.getProperty("enumerateDegenMotifs"));
		boolean enumeratePossibleDegenMotifs = Boolean.parseBoolean(params.getProperty("enumeratePossibleDegenMotifs"));
		boolean mapMotifs = Boolean.parseBoolean(params.getProperty("mapMotifs"));
		
		int motifLength = Integer.parseInt(params.getProperty("motifLength"));
		int maxDegenThreshold = Integer.parseInt(params.getProperty("maxDegenThreshold"));
		
		String wd = params.getProperty("working_directory");
		
		/* Command line arguments */ 
		String motifsToTestFile = args[1];
		String degenMotifsToRefSeqFile = args[2];

		/* Local computer - file paths */
		String projectName = params.getProperty("project_name");
		
		/* Input Files */
		String fastaFile = wd + params.getProperty("fastaFile").replaceAll("\\s+", ""); // output from UCSC genome browser
		String mapProteinToRefSeqFile = wd + params.getProperty("mapGeneSymbolsToRefSeqIds").replaceAll("\\s+", "");// output from BiomaRt
		
		/* Output Files */
		String mapMotifsToRefSeqIdsFile = wd + "MapMotifsToRefSeqIds"+ projectName +".tsv"; // output from motif enumeration
		String listOfUniqueMotifsFile = wd + "ListOfUniqueMotifs_"+ projectName + ".txt"; // output from motif enumeration
		String motifMapFile = wd + "motifMappedToProteinsInNetwork"+ projectName +".tsv"; // output from motif enumeration
		
		String mapMotifsToDegenMotifsFilePath = wd + "MapOfMotifsToDegenMotifs_"+ projectName +".tsv"; // output from motif degeneration
		
		String motifSet = wd + motifsToTestFile; // if enumerating degen motifs = list of non degen motifs, if mapping degen motifs = list of possible degen motifs
		String mapOfMotifs = wd + degenMotifsToRefSeqFile;
		
		/* Generate mapping of protein HGNC symbols to mRNA RefSeqIds >> To call R */
		
		// MOTIF ENUMERATION CAN BE RUN LOCALLY // 
		if(enumerateMotifs) {
			System.out.println("**Enumerating motifs**");
			
			/* Load list of RefSeq IDs in network */ 
			HashSet<String> refSeqSet = SeqLoader.loadRefSeqIDsToTest(mapProteinToRefSeqFile);
			System.out.println("Number of loaded ref seq Ids: " + refSeqSet.size() + "\n");

			/* Enumerate all unique motifs within the network */ 
			MotifEnumerator e = new MotifEnumerator(fastaFile, motifLength, refSeqSet);
			e.generateMotifList(mapMotifsToRefSeqIdsFile, listOfUniqueMotifsFile);
			
			System.out.println("**Generating annotation file for non degenerate motifs**");
			
			/* Map motifs to the proteins in the network */
			MapMotifs.mapMotifsToRefSeqIds(mapMotifsToRefSeqIdsFile, mapProteinToRefSeqFile, motifMapFile);
		}	
		
		// MISSING STEP : Split # of motif in listOfUniqueMotifsToTest into a max of 1000 files ??? 
		
		if(enumeratePossibleDegenMotifs) {
			System.out.println("**Generating all possible degen motifs**");
			MotifDegeneration d1 = new MotifDegeneration(motifLength, maxDegenThreshold);
			d1.generateAllPossibleMotifs(motifSet);
		}
		
		// MOTIF DEGENERATION SHOULD BE RUN REMOTELY //
		if(enumerateDegenMotifs) {
			System.out.println("**Enumerating degenerate motifs from non degenerate motifs**");
			MotifDegeneration d = new MotifDegeneration(motifLength, maxDegenThreshold);
			d.enumerateDegenerateMotifs(motifSet, mapOfMotifs);
		}
		
		// Update to run with multiple files
		if(mapMotifs) {
			System.out.println("**Mapping degen motifs to proteins in network**");
			MapMotifs.mapDegenMotifsToRefSeqIds(mapMotifsToRefSeqIdsFile, motifSet, mapMotifsToDegenMotifsFilePath, mapOfMotifs, mapProteinToRefSeqFile);
			
		}
		
	}

}
