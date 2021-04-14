package motifs;

import java.io.File;
import java.util.HashSet;

public class MotifMain {

	public static void main(String[] args) {

		boolean runRemote = false;
		/* Command line arguments */ 
		String motifsToTestFileName = args[0];
		String degenMotifsFileName = args[1];

		/* Local computer - file paths */
		String wd = "C:\\Users\\Rachel\\Documents\\LESMoNlocal\\";

		String fastaFile = wd + "input_files\\human_3UTRsequences.txt";
		String rnaIdListFile = wd + "IO_files\\refSeqRNAids-HumanCellMap.tsv";

		String mapMotifsToRefSeqIdsFile = wd + "IO_files\\mapMotifsToRefSeqIdsFile.tsv";

		String motifsToTestFilePath = wd + motifsToTestFileName;
		String motifsTodegenMotifsFilePath = wd + degenMotifsFileName;
		
		String degenMotifSet = wd + "IO_files\\degenMotifSet_sample.tsv";
		String degenMotifToRedSeqFile = wd + "IO_files\\motifAnnotationFile.tsv";
		/* Remote computer - file paths */ 
		if(runRemote) {
			wd = "/home/rnade046/projects/rrg-mlaval/rnade046/motifDegen_Full_FWD/";

			fastaFile = wd + "input_files\\human_3UTRsequences.txt";
			rnaIdListFile = wd + "IO_files\\refSeqRNAids-HumanCellMap.tsv";

			mapMotifsToRefSeqIdsFile = wd + "IO_files\\mapMotifsToRefSeqIdsFile.tsv";

			motifsToTestFilePath = wd + motifsToTestFileName;
			motifsTodegenMotifsFilePath = wd + degenMotifsFileName;
			degenMotifSet = wd + "degenMotifSet.tsv";
			degenMotifToRedSeqFile = wd + "motifAnnotationFile.tsv";
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
		
		File f1 = new File(motifsTodegenMotifsFilePath);
		if(!f1.exists() && !f1.isDirectory()) {
			MotifDegeneration d = new MotifDegeneration(motifLength, maxDegenThreshold);
			d.enumerateDegenerateMotifs(motifsToTestFilePath, motifsTodegenMotifsFilePath);
			d.generateAllPossibleMotifs(degenMotifSet);
		}
		
		MapMotifs.mapMotifsToRefSeqIds(mapMotifsToRefSeqIdsFile, degenMotifSet, motifsTodegenMotifsFilePath, degenMotifToRedSeqFile);
	}

}
