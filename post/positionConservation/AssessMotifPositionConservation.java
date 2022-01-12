package positionConservation;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Properties;

public class AssessMotifPositionConservation {

	public static void main(String[] args) throws FileNotFoundException, IOException {

		System.out.println("**Loading parameters file** \n");
		Properties params = new Properties();
		params.load(new FileInputStream(args[0]));		

		String wd = params.getProperty("working_directory");
		String networkName = params.getProperty("network_name");

		int clusteringMeasure = Integer.parseInt(params.getProperty("clusteringMeasure", "0"));
		double percentThreshold = Double.parseDouble(params.getProperty("percentThreshold", "0.2"));

		String clusteringName = "";

		switch(clusteringMeasure) {
		case 0: clusteringName = "_TPD";
		break;
		case 1: clusteringName = "_TPPD_p" + percentThreshold;
		break;
		case 2: clusteringName = "_coreTPD_p" + percentThreshold;
		break;
		}
		
		File directory = new File(wd + "/MotifPosition/"); 
		if (! directory.exists()){
			System.out.println("creating directory: MotifPosition/");
			directory.mkdir();
		}

		String protAnnotationFreqFile = wd + networkName + "_protFreqAnnotation.tsv";
		String proteinToRefSeqIdFile = wd +  networkName + "_proteinsInNetwork_info.tsv";

		String extractedAnnotationsFile = wd + networkName + clusteringName + "_annotationSubset.tsv";
		String corePorteinsFile = wd + networkName + clusteringName +"_coreProteinsByMotif.tsv";

		String fastaFile = wd + params.getProperty("fastaFile");
		String motifOutputPrefixFile = wd + "MotifPosition/motifPositionConservation_"; 
		PositionConservation p = new PositionConservation(fastaFile, proteinToRefSeqIdFile, protAnnotationFreqFile, 8);
		p.getMotifPositions(extractedAnnotationsFile, motifOutputPrefixFile, Integer.parseInt(args[1]), Integer.parseInt(args[2]));
	
		/* motif positions of core proteins */
		if(clusteringMeasure == 1 || clusteringMeasure == 2) {
			motifOutputPrefixFile = wd + "MotifPositionConservation_coreProteins_";
			p.getMotifPositions(corePorteinsFile, motifOutputPrefixFile, Integer.parseInt(args[1]), Integer.parseInt(args[2]));				
			
		}
		
	}

}
