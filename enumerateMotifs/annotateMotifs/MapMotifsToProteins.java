package annotateMotifs;

import java.io.File;

public class MapMotifsToProteins {

	public static void main(String[] args) {
		
		String wd = args[0];
		String fastaFile = wd + "/input_files/human_3UTRsequences.txt";
		String degenMotifSetPrefix = wd + "/motif_enumeration/degenMotifSet/degenMotifs_";
		String motifDir = wd + "/motif_enumeration/degenMotifSet/";
		String ids = wd + "/input_files/corrNetTop2-400_proteinsInNetwork_info.tsv";
		
		String annotationFile = wd + "/motif_enumeration/annotationFile/corrNet_degenAnnotations_FWD_";
		
		File dir = new File(wd + "/motif_enumeration/annotationFile/");
		if (! dir.exists()){
			dir.mkdir();
		}
		
		MotifMapper m = new MotifMapper(ids, annotationFile, fastaFile, motifDir, degenMotifSetPrefix);
		m.mapMotifsToTheirAssociatedProteins();
	}

}
