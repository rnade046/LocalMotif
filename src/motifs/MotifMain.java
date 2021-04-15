package motifs;

import java.io.File;
import java.util.HashSet;

public class MotifMain {

	public static void main(String[] args) {

		boolean runRemote = true;
		/* Command line arguments */ 
		String degenMotifsToTestFile = args[0];
		String degenMotifsToRefSeqFile = args[1];

		/* Local computer - file paths */
		String wd = "C:\\Users\\Rachel\\Documents\\LESMoNlocal\\";

		String fastaFile = wd + "input_files\\human_3UTRsequences.txt";
		String rnaIdListFile = wd + "IO_files\\refSeqRNAids-HumanCellMap.tsv";

		String mapMotifsToRefSeqIdsFile = wd + "IO_files\\mapMotifsToRefSeqIdsFile.tsv";
		String refSeqToProtMapFile = wd + "IO_files\\protToRefSeqRNAids-HumanCellMap.tsv";

		String motifsToTestFilePath = wd + "IO_files\\test-motifs.txt";
		String mapOfMotifsToDegenMotifsFilePath = wd + "IO_files\\degenMotifs_sample.tsv";
		
		String degenMotifSet = wd + "IO_files\\degenMotifsToMap_000";
		String degenMotifToRedSeqFile = wd + "IO_files\\motifAnnotationFile.tsv";
		/* Remote computer - file paths */ 
		if(runRemote) {
			wd = "/home/rnade046/projects/rrg-mlaval/rnade046/motifDegen_Full_FWD/";

			fastaFile = wd + "input_files\\human_3UTRsequences.txt";
			rnaIdListFile = wd + "IO_files\\refSeqRNAids-HumanCellMap.tsv";

			mapMotifsToRefSeqIdsFile = wd + "mapMotifsToRefSeqIdsFile.tsv";
			refSeqToProtMapFile = wd + "protToRefSeqRNAids-HumanCellMap.tsv";
			
			motifsToTestFilePath = wd + "motifsToTestFile.tsv";
			mapOfMotifsToDegenMotifsFilePath = wd + "degenerateMotifs_FWD.tsv";
			degenMotifSet = wd + degenMotifsToTestFile;
			degenMotifToRedSeqFile = wd + degenMotifsToRefSeqFile;
		}

		int motifLength = 8;
		int maxDegenThreshold = 7;

		File f = new File(motifsToTestFilePath);
		if(!f.exists() && !f.isDirectory()) {
			/* Load list of RefSeq IDs in network */ 
			HashSet<String> refSeqSet = SeqLoader.loadRefSeqIDsToTest(rnaIdListFile);
			System.out.println("Number of loaded ref seq Ids: " + refSeqSet.size());

			/* Enumerate all unique motifs within the network */ 
			MotifEnumerator e = new MotifEnumerator(fastaFile, motifLength, refSeqSet);
			e.generateMotifList(mapMotifsToRefSeqIdsFile, motifsToTestFilePath);
		}	
		
		File f1 = new File(mapOfMotifsToDegenMotifsFilePath);
		if(!f1.exists() && !f1.isDirectory()) {
			MotifDegeneration d = new MotifDegeneration(motifLength, maxDegenThreshold);
			d.enumerateDegenerateMotifs(motifsToTestFilePath, mapOfMotifsToDegenMotifsFilePath);
			d.generateAllPossibleMotifs(degenMotifSet);
		}
		
		MapMotifs.mapMotifsToRefSeqIds(mapMotifsToRefSeqIdsFile, degenMotifSet, mapOfMotifsToDegenMotifsFilePath, degenMotifToRedSeqFile, refSeqToProtMapFile);
	}

}
