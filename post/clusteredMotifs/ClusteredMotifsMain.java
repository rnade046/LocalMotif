package clusteredMotifs;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.Properties;

public class ClusteredMotifsMain {

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
		String corePorteinsFile = wd + networkName + clusteringName +"_coreProteinsByMotif.tsv";

		String proteinToRefSeqIdFile = wd +  networkName + "_proteinsInNetwork_info.tsv";
		String enumeratedMotifs = wd +  "motifFamilies/" + networkName + "_enumeratedMotifsPerRefSeqId.tsv";

		
		/* Identify motifs that pass significant threshold: (1) print details to separate file, (2) store motif and file # in map */ 
		File f = new File(significantMotifsFile);
		if(!f.exists() && !f.isDirectory()) { 
			System.out.println("**Identifying significant motifs**");
			IdentifyMotifs.getSignificantMotifs(motifClusteringPrefix, numOfFiles, pvalThreshold, significantMotifsFile);
		}

		System.out.println("**Loading significant motifs**");
		HashMap<String, Integer> motifMapOfFileIdx = IdentifyMotifs.loadSignificantMotifs(significantMotifsFile);
		System.out.println("Number of significant motifs: " + motifMapOfFileIdx.size() + "\n");
		
		/* Get annotations (motif = protein1|protein2|...|proteinN) of significant motifs for easy access in future analysis */ 
		File f2 = new File(extractedAnnotationsFile);
		if(!f2.exists() && !f2.isDirectory()) { 
			System.out.println("**Identifying annotated proteins**");
			IdentifyMotifs.getAnnotatedProteinInfo(motifMapOfFileIdx, protAnnotationFreqFile, extractedAnnotationsFile,	annotationPrefixFile);
		}

		File f3 = new File(corePorteinsFile);
		if(!f3.exists() && !f3.isDirectory()) {
			
			/* For CoreTPD and TPPD; load annotation subset and determine core proteins; print to seperate file */
			if(clusteringMeasure == 1 || clusteringMeasure == 2) {
				
				System.out.println("**Identifying Core Proteins**");
				
				String distanceMatrixFile = wd + networkName + "_removedOverConnectedProteins_" + params.getProperty("maxInteractions") + "_distanceMatrix2.txt";

				IdentifyCoreProteins.getCoreProteins(extractedAnnotationsFile, Double.parseDouble(params.getProperty("percentThreshold")), 
						proteinToRefSeqIdFile, distanceMatrixFile, corePorteinsFile);
			}
		}
		
		/* Compute similarity between significant proteins for determination of motif family*/
		if(Boolean.parseBoolean(params.getProperty("computeSimilarity"))) {

			File directory2 = new File(wd + "/motifFamilies/"); 
			if (! directory2.exists()){
				System.out.println("creating directory: motifFamilies/");
				directory2.mkdir();
			}
			
			String motifsInMatrixFile = wd +  "motifFamilies/" + networkName + clusteringName + "_motifsMatrix_p" + pvalThreshold + ".tsv";
			String similarityMatrix = wd + "motifFamilies/" +  networkName + clusteringName + "_similarity_DistanceMatrix_p" + pvalThreshold + ".tsv" ;
			
			System.out.println("**Loading annotation info**");
			/* Search through annotation Files to get proteins annotated by significant motifs: (1) store in map for similarity measuring, (2) print to file for local testing */
			//HashMap<String, String[]> motifMapOfAnnotatedProteins = IdentifyMotifs.getAnnotatedProteinInfo(motifMapOfFileIdx, annotationPrefixFile);
			HashMap<String, String[]> motifMapOfAnnotatedProteins = IdentifyMotifs.loadAnnotatedProteinInfo(extractedAnnotationsFile);
			System.out.println("Found motif info: " + motifMapOfAnnotatedProteins.size() + "\n");
			System.out.println("**Computing similarity**");
			Similarity.computeMotifSimilary(motifMapOfAnnotatedProteins, motifsInMatrixFile, similarityMatrix);
			
			/* For CoreTPD and TPPD; load annotation subset and determine core proteins; print */
			if(clusteringMeasure == 1 || clusteringMeasure == 2) {
				
				System.out.println("**Similarity calculation for Core Proteins**");
				
				motifsInMatrixFile = wd +  "motifFamilies/" + networkName + clusteringName + "_CorePorteins_motifsMatrix_p" + pvalThreshold + ".tsv";
				similarityMatrix = wd + "motifFamilies/" +  networkName + clusteringName + "_CoreProteins_similarity_DistanceMatrix_p" + pvalThreshold + ".tsv" ;
				
				motifMapOfAnnotatedProteins = IdentifyMotifs.loadAnnotatedProteinInfo(corePorteinsFile);
				System.out.println("Found motif info: " + motifMapOfAnnotatedProteins.size() + "\n");
				System.out.println("**Computing similarity**");
				Similarity.computeMotifSimilary(motifMapOfAnnotatedProteins, motifsInMatrixFile, similarityMatrix);

			}
		}

		/* output to R : perform hierarchical clustering*/

		/* Assess motif families */ 
		if(Boolean.parseBoolean(params.getProperty("assessMotifFamilies"))) {
			
			String motifInstancesPrefix = wd + "motifFamilies/" + networkName + clusteringName + "_motifInstances_motifFamily";
			String motifPPMPrefix = wd +  "motifFamilies/" +  networkName + clusteringName + "_ppm_motifFamilyGroup";
			String motifInfoFile = wd +  "motifFamilies/" + networkName + clusteringName + "_motifFamiliesInfo.tsv";

			String motifFamilyFilePrefix = wd +  "motifFamilies/" + "motifFamily_CoreTPD40_ward2_group";
			int numberOfFamilies = Integer.parseInt(params.getProperty("motifFamilyGroups", "10"));
			
			System.out.println("**Assessing motif families**");
			MotifFamily.assessMotifFamilies(motifFamilyFilePrefix, numberOfFamilies, significantMotifsFile, enumeratedMotifs, proteinToRefSeqIdFile, motifInstancesPrefix, motifPPMPrefix, motifInfoFile, extractedAnnotationsFile);

			if(clusteringMeasure == 1 || clusteringMeasure == 2) {
				
				motifInstancesPrefix = wd + "motifFamilies/" + networkName + clusteringName + "_coreProteins_motifInstances_motifFamily";
				motifPPMPrefix = wd +  "motifFamilies/" +  networkName + clusteringName + "_coreProteins_ppm_motifFamilyGroup";
				motifInfoFile = wd +  "motifFamilies/" + networkName + clusteringName + "_coreProteins_motifFamiliesInfo.tsv";
				
				motifFamilyFilePrefix = wd +  "motifFamilies/" + "motifFamily_CoreTPD40_CoreProteins_ward2_group";
				int numberOfCoreFamilies = Integer.parseInt(params.getProperty("coreFamilyGroups", "10"));
				
				System.out.println("**Assessing motif instances for Core Proteins**");
				MotifFamily.assessMotifFamilies(motifFamilyFilePrefix, numberOfCoreFamilies, significantMotifsFile, enumeratedMotifs, proteinToRefSeqIdFile, motifInstancesPrefix, motifPPMPrefix, motifInfoFile, corePorteinsFile);

			}
			
		}

	}

}
