import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.Properties;

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
		String projectName = params.getProperty("project_name");

		int clusteringMeasure = Integer.parseInt(params.getProperty("clusteringMeasure", "0"));
		double percentThreshold = Double.parseDouble(params.getProperty("percentThreshold", "0.2"));

		String clusteringName = "";

		switch(clusteringMeasure) {
		case 0: clusteringName = "";
		break;
		case 1: clusteringName = "_TPPD" + percentThreshold;
		break;
		case 2: clusteringName = "_coreTPD_p" + percentThreshold;
		break;
		}
		
		String motifClusteringPrefix = wd + "motifClustering/" + projectName + clusteringName + "_testedDegenMotifClustering_";
		int numOfFiles = Integer.parseInt(params.getProperty("numDegenMotifFiles"));

		double pvalThreshold = Double.parseDouble(params.getProperty("significantThreshold"));
		
		String significantMotifsFile = wd + "motifFamilies/" +  projectName + clusteringName +"_signficantMotifs_p" + pvalThreshold +".tsv";

		//String annotationPrefixFile = params.getProperty("degenAnnotationPrefix");
		String extractedAnnotationsFile = wd +  "motifFamilies/" + projectName + clusteringName + "_annotationSubset.tsv";
		
		String motifsInMatrixFile = wd +  "motifFamilies/" + projectName + clusteringName + "_motifsMatrix_p" + pvalThreshold + ".tsv";
		String similarityMatrix = wd + "motifFamilies/" +  projectName + clusteringName + "_similarity_DistanceMatrix_p" + pvalThreshold + ".tsv" ;
		
		String motifFamilyFilePrefix = wd +  "motifFamilies/" + "motifFamily_ward_group";
		int numberOfFamilies = Integer.parseInt(params.getProperty("motifFamilyGroups"));
		
		String enumeratedMotifs = wd +  "motifFamilies/" + projectName + "_enumeratedMotifsPerRefSeqId.tsv";
		String proteinToRefSeqIdFile = wd + "motifFamilies/" +  projectName + "_proteinsInNetwork_info.tsv";
		
		String motifInstancesPrefix = wd +  "motifFamilies/" +  projectName + "_ppm_motifFamilyGroup";
		String motifInfoFile = wd +  "motifFamilies/" + projectName + "_motifFamiliesInfo.tsv";
		
		//System.out.println("**Identifying significant motifs**");
		/* Identify motifs that pass significant threshold: (1) print details to separate file, (2) store motif and file # in map */ 
		//HashMap<String, Integer> motifMapOfFileIdx = IdentifyMotifs.getSignificantMotifs(motifClusteringPrefix, numOfFiles, pvalThreshold, significantMotifsFile);
		//HashMap<String, Integer> motifMapOfFileIdx = IdentifyMotifs.loadSignificantMotifs(significantMotifsFile);
		//System.out.println("Number of significant motifs: " + motifMapOfFileIdx.size() + "\n");
		
		
		if(Boolean.parseBoolean(params.getProperty("computeSimilarity"))) {
			
			System.out.println("**Loading annotation info**");
			/* Search through annotation Files to get proteins annotated by significant motifs: (1) store in map for similarity measuring, (2) print to file for local testing */
			//HashMap<String, String[]> motifMapOfAnnotatedProteins = IdentifyMotifs.getAnnotatedProteinInfo(motifMapOfFileIdx, annotationPrefixFile);
			HashMap<String, String[]> motifMapOfAnnotatedProteins = IdentifyMotifs.getAnnotatedProteinInfoForTesting(extractedAnnotationsFile);
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
