package annotateMotifs;

import java.io.File;

public class MapMotifsToProteins {

	public static void main(String[] args) {
		
		int motifLength = 8;
		int maxDegenThreshold = 7;

		String wd = args[0];

		File dir = new File(wd + "/motif_enumeration/");
		if(! dir.exists()) {
			dir.mkdir();
		}

		File dir3 = new File(wd + "/motif_enumeration/degenMotifSet");
		if (! dir3.exists()){
			dir3.mkdir();
		}

		String degenMotifSetPrefix = wd + "/motif_enumeration/degenMotifSet/degenMotifs_";
	
		if(dir3.list().length < 1000) {
			System.out.println("**Generating all possible degen motifs**");
			MotifDegeneration d1 = new MotifDegeneration(motifLength, maxDegenThreshold);
			d1.generateAllPossibleMotifs(degenMotifSetPrefix);
		}
		
		String fastaFile = wd + "/input_files/human_3UTRsequences.txt";
		String motifDir = wd + "/motif_enumeration/degenMotifSet/";
		String ids = wd + "/input_files/corrNetTop2-400_proteinsInNetwork_info.tsv";
		
		String annotationFile = wd + "/motif_enumeration/annotationFile/corrNet_degenAnnotations_FWD_";
		
		File dir2 = new File(wd + "/motif_enumeration/annotationFile/");
		if (! dir2.exists()){
			dir2.mkdir();
		}
		
		MotifMapper m = new MotifMapper(ids, annotationFile, fastaFile, motifDir, degenMotifSetPrefix);
		m.mapMotifsToTheirAssociatedProteins();
	}

}
