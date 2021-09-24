import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.Properties;

import ClusteredMotifs.FunctionalEnrichment;
import ClusteredMotifs.IdentifyMotifs;
import ClusteredMotifs.MotifFamily;
import ClusteredMotifs.Similarity;

public class MainResults {

	/**
	 * Note: making changes for local testing
	 */
	
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
		
		String motifClusteringPrefix = wd + "motifClustering/" + networkName + clusteringName + "_testedDegenMotifClustering_";
		int numOfFiles = Integer.parseInt(params.getProperty("numDegenMotifFiles"));

		double pvalThreshold = Double.parseDouble(params.getProperty("significantThreshold"));
		
		String significantMotifsFile = wd +  networkName + clusteringName +"_signficantMotifs_p" + pvalThreshold +".tsv";

		String annotationPrefixFile = params.getProperty("degenAnnotationPrefix") + "corrNetTop2_degenMotifMappedToProteinsInNetwork_";
		String protAnnotationFreqFile = wd + networkName + "_protFreqAnnotation.tsv";
		
		String extractedAnnotationsFile = wd + networkName + clusteringName + "_annotationSubset.tsv";
		
		String motifsInMatrixFile = wd +  "motifFamilies/" + networkName + clusteringName + "_motifsMatrix_p" + pvalThreshold + ".tsv";
		String similarityMatrix = wd + "motifFamilies/" +  networkName + clusteringName + "_similarity_DistanceMatrix_p" + pvalThreshold + ".tsv" ;
		
		String motifFamilyFilePrefix = wd +  "motifFamilies/" + "motifFamily_ward_group";
		int numberOfFamilies = Integer.parseInt(params.getProperty("motifFamilyGroups", "10"));
		
		String enumeratedMotifs = wd +  "motifFamilies/" + networkName + "_enumeratedMotifsPerRefSeqId.tsv";
		String proteinToRefSeqIdFile = wd + "motifFamilies/" +  networkName + "_proteinsInNetwork_info.tsv";
		
		String motifInstancesPrefix = wd +  "motifFamilies/" +  networkName + "_ppm_motifFamilyGroup";
		String motifInfoFile = wd +  "motifFamilies/" + networkName + "_motifFamiliesInfo.tsv";
		
		
		/* Identify motifs that pass significant threshold: (1) print details to separate file, (2) store motif and file # in map */ 
		File f = new File(significantMotifsFile);
		if(!f.exists() && !f.isDirectory()) { 
			System.out.println("**Identifying significant motifs**");
			IdentifyMotifs.getSignificantMotifs(motifClusteringPrefix, numOfFiles, pvalThreshold, significantMotifsFile);
		}
		
		System.out.println("**Loading significant motifs**");
		HashMap<String, Integer> motifMapOfFileIdx = IdentifyMotifs.loadSignificantMotifs(significantMotifsFile);
		System.out.println("Number of significant motifs: " + motifMapOfFileIdx.size() + "\n");
		
		File f2 = new File(extractedAnnotationsFile);
		if(!f2.exists() && !f2.isDirectory()) { 
			System.out.println("**Identifying annotated proteins**");
			IdentifyMotifs.getAnnotatedProteinInfo(motifMapOfFileIdx, protAnnotationFreqFile, extractedAnnotationsFile,	annotationPrefixFile);
		}
		
		if(Boolean.parseBoolean(params.getProperty("goEnrich"))) {
			
			File directory = new File(wd + "/Ontologizer/"); 
			if (! directory.exists()){
				System.out.println("creating directory: Ontologizer/");
				directory.mkdir();
			}
			
			String proteinsInNetworkFile = wd + "Ontologizer/" + networkName + "_proteinsInNetwork.txt";
			String annotatedProteinsPrefix = wd + "Ontologizer/" + networkName + clusteringName + "_annotatedProteinsByMotif_";
			System.out.println("**Formatting files for ontologizer analysis**");
			FunctionalEnrichment.formatFilesForOntologizer(protAnnotationFreqFile, extractedAnnotationsFile, proteinsInNetworkFile, annotatedProteinsPrefix);
		}
		
		
		if(Boolean.parseBoolean(params.getProperty("computeSimilarity"))) {
			
			File directory2 = new File(wd + "/motifFamilies/"); 
			if (! directory2.exists()){
				System.out.println("creating directory: motifFamilies/");
				directory2.mkdir();
			}
			
			System.out.println("**Loading annotation info**");
			/* Search through annotation Files to get proteins annotated by significant motifs: (1) store in map for similarity measuring, (2) print to file for local testing */
			//HashMap<String, String[]> motifMapOfAnnotatedProteins = IdentifyMotifs.getAnnotatedProteinInfo(motifMapOfFileIdx, annotationPrefixFile);
			HashMap<String, String[]> motifMapOfAnnotatedProteins = IdentifyMotifs.loadAnnotatedProteinInfo(extractedAnnotationsFile);
			System.out.println("Found motif info: " + motifMapOfAnnotatedProteins.size() + "\n");
			System.out.println("**Computing similarity**");
			Similarity.computeMotifSimilary(motifMapOfAnnotatedProteins, motifsInMatrixFile, similarityMatrix);
		}
		
		/* output to R : perform hierarchical clustering*/
		
		/* Assess motif families */ 
		if(Boolean.parseBoolean(params.getProperty("assessMotifFamilies"))) {
			System.out.println("**Assessing motif families**");
			MotifFamily.assessMotifFamilies(motifFamilyFilePrefix, numberOfFamilies, significantMotifsFile, 1, enumeratedMotifs, proteinToRefSeqIdFile, motifInstancesPrefix, motifInfoFile);

		}
		
	}

	
}
