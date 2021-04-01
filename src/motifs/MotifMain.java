package motifs;

import java.util.HashSet;

public class MotifMain {

	public static void main(String[] args) {
	
		String wd = "C:\\Users\\Rachel\\Documents\\LESMoNlocal\\";
		
		String fastaFile = wd + "input_files\\human_3UTRsequences.txt";
		String rnaIdListFile = wd + "IO_files\\refSeqRNAids-HumanCellMap.tsv";
		
		String mapMotifsToRefSeqIdsFile = wd + "IO_files\\mapMotifsToRefSeqIdsFile.tsv";
		String motifsToTestFile = wd + "IO_files\\motifsToTestFile.tsv";
		
		// Load list of RefSeq IDs in network
		HashSet<String> refSeqSet = SeqLoader.loadRefSeqIDsToTest(rnaIdListFile);
		System.out.println("Number of loaded ref seq Ids: " + refSeqSet.size());
		
		MotifEnumerator e = new MotifEnumerator(fastaFile, 8, refSeqSet);
		//HashSet<String> motifSet = e.generateMotifList(mapMotifsToRefSeqIdsFile, motifsToTestFile);
		e.generateMotifList(mapMotifsToRefSeqIdsFile, motifsToTestFile);
		// Get an index of the different Sequences; RefSeqId = lineCount 
		//HashMap<String, Integer> idxOfSeqInFastaMap = SeqLoader.generateIndexOfFastaFile(fastaFile);
		
		// For a given ID in network; load its formated sequence, enumerate motif nmers
		// HashMap <RefSeq ID, List<motif>> && print to file
		
		// Make Set of unique motifs to degenerate
		
		
		
		
	}
	
}
