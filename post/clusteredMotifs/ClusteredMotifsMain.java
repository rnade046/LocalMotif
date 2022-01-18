package clusteredMotifs;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
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

			String familyFolder = networkName + clusteringName + "_p" + pvalThreshold + "/";
			File dir3 = new File(wd + "/motifFamilies/" + familyFolder);
			if (! dir3.exists()){
				System.out.println("creating directory: " + familyFolder + "\n");
				dir3.mkdir();
			}

			String motifsInMatrixFile = wd +  "motifFamilies/" + familyFolder + networkName + clusteringName + "_p" + pvalThreshold + "_MotifsInMatrix.tsv";
			String similarityMatrix = wd + "motifFamilies/"+ familyFolder +  networkName + clusteringName + "_p" + pvalThreshold + "_DistanceMatrix.tsv" ;
			f3 = new File(similarityMatrix);
			if(!f3.exists() && !f3.isDirectory()) { 
				System.out.println("**Loading annotation info**");
				/* Search through annotation Files to get proteins annotated by significant motifs: (1) store in map for similarity measuring, (2) print to file for local testing */
				//HashMap<String, String[]> motifMapOfAnnotatedProteins = IdentifyMotifs.getAnnotatedProteinInfo(motifMapOfFileIdx, annotationPrefixFile);
				HashMap<String, String[]> motifMapOfAnnotatedProteins = IdentifyMotifs.loadAnnotatedProteinInfo(extractedAnnotationsFile);
				System.out.println("Found motif info: " + motifMapOfAnnotatedProteins.size() + "\n");


				System.out.println("**Computing similarity**");
				Similarity.computeMotifSimilary(motifMapOfAnnotatedProteins, motifsInMatrixFile, similarityMatrix);
			}
			/* For CoreTPD and TPPD; load annotation subset and determine core proteins; print */
			if(clusteringMeasure == 1 || clusteringMeasure == 2) {



				motifsInMatrixFile = wd +  "motifFamilies/" + familyFolder + networkName + clusteringName + "_CoreProteins_p" + pvalThreshold + "_MotifsInMatrix.tsv";
				similarityMatrix = wd + "motifFamilies/" + familyFolder +  networkName + clusteringName + "_CoreProteins_p" + pvalThreshold + "_DistanceMatrix.tsv" ;

				f3 = new File(similarityMatrix);
				if(!f3.exists() && !f3.isDirectory()) { 
					System.out.println("**Similarity calculation for Core Proteins**");
					HashMap<String, String[]> motifMapOfAnnotatedProteins = IdentifyMotifs.loadAnnotatedProteinInfo(corePorteinsFile);
					System.out.println("Found motif info: " + motifMapOfAnnotatedProteins.size() + "\n");

					System.out.println("**Computing similarity**");
					Similarity.computeMotifSimilary(motifMapOfAnnotatedProteins, motifsInMatrixFile, similarityMatrix);
				}
			}
			/* output to R : perform hierarchical clustering*/
			String wd4r = wd + "motifFamilies/" + familyFolder;
			String projectName = networkName + clusteringName;
			Runtime.getRuntime().exec("Rscript HierarchicalClustering_1.R " + wd4r + " " + projectName + " " + pvalThreshold + " ward.D2"); 
			//Runtime.getRuntime().exec("Rscript HierarchicalClustering_1.R C:\\Users\\Rachel\\Documents\\LESMoNlocal\\analysis\\motifFamilies/ corrNetTop2-400_coreTPD_p0.4 5.41545109270352E-7 ward.D2");

			// run second R script when ready (e.g. testing h values [min max interval])
			double[] hTest = Arrays.stream(params.getProperty("heightToTest", "0").split("\\s+")).mapToDouble(Double::parseDouble).toArray();
			if(hTest[0] != 0) {
				Runtime.getRuntime().exec("Rscript HierarchicalClustering_2.R " + wd4r + " " + projectName + " " + pvalThreshold + " ward.D2 " + hTest[0] + " " + hTest[1] + " " + hTest[2]); 
			}

			double height = Double.parseDouble(params.getProperty("height", "0"));
			if(height != 0) {

				File dir4 = new File(wd + "/motifFamilies/" + familyFolder + "/Groups_h" + height + "/");
				if (! dir4.exists()){
					System.out.println("creating directory: " + familyFolder + "\n");
					dir4.mkdir();
				}
				String condition = "Groups";
				Runtime.getRuntime().exec("Rscript HierarchicalClustering_3.R " + wd4r + " " + projectName + " " + pvalThreshold + " " + condition + " "+ height  + " ward.D2"); 
			}

			if(clusteringMeasure == 1 || clusteringMeasure == 2) {

				/* output to R : perform hierarchical clustering for CoreProteins */
				projectName = networkName + clusteringName + "_CoreProteins";
				Runtime.getRuntime().exec("Rscript HierarchicalClustering_1.R " + wd4r + " " + projectName + " " + pvalThreshold + " ward.D2"); 

				// run second R script when ready (e.g. testing h values [min max interval])
				double[] hTest2 = Arrays.stream(params.getProperty("coreHeightToTest", "0").split("\\s+")).mapToDouble(Double::parseDouble).toArray();
				if(hTest[0] != 0) {
					Runtime.getRuntime().exec("Rscript HierarchicalClustering_2.R " + wd4r + " " + projectName + " " + pvalThreshold + " ward.D2 " + hTest2[0] + " " + hTest2[1] + " " + hTest2[2]); 

				}

				// run second R script when ready (e.g. h is set)
				double height2 = Double.parseDouble(params.getProperty("coreHeight", "0"));
				if(height2 != 0) {

					String condition = "Groups_CoreProteins";
					File dir4 = new File(wd + "/motifFamilies/" + familyFolder + condition+ "_h" + height2 + "/");
					if (! dir4.exists()){
						System.out.println("creating directory: " + familyFolder + "\n");
						dir4.mkdir();
					}

					Runtime.getRuntime().exec("Rscript HierarchicalClustering_3.R " + wd4r + " " + projectName + " " + pvalThreshold  + " " + condition + " " + height2  + " ward.D2"); 
				}
			}
		}
		/* Assess motif families */ 
		if(Boolean.parseBoolean(params.getProperty("assessMotifFamilies"))) {
			
			String familyFolder = networkName + clusteringName + "_p" + pvalThreshold + "/";

			double height = Double.parseDouble(params.getProperty("height", "0"));
			String condition = "Groups_h" + height + "/";

			String motifInstancesPrefix = wd + "motifFamilies/" + familyFolder + condition + networkName + clusteringName + "_h" + height + "_motifInstances_motifFamily";
			String motifPPMPrefix = wd +  "motifFamilies/" + familyFolder + condition + networkName + clusteringName+ "_h" + height + "_ppm_motifFamilyGroup";
			String motifInfoFile = wd +  "motifFamilies/" + familyFolder + condition + networkName + clusteringName + "_h" + height+ "_motifFamiliesInfo.tsv";


			String motifFamilyFilePrefix = wd +  "motifFamilies/" + familyFolder + condition + networkName + clusteringName + "_p" + pvalThreshold + "_ward.D2_h" + height +"_group";

			int numberOfFamilies = Integer.parseInt(params.getProperty("motifFamilyGroups", "10"));

			System.out.println("**Assessing motif families**");
			MotifFamily.assessMotifFamilies(motifFamilyFilePrefix, numberOfFamilies, significantMotifsFile, enumeratedMotifs, proteinToRefSeqIdFile, motifInstancesPrefix, motifPPMPrefix, motifInfoFile, extractedAnnotationsFile);
			
			String wd4r = wd + "motifFamilies/" + familyFolder + condition;
			String projectName = networkName + clusteringName;
			
			Runtime.getRuntime().exec("Rscript seqLogo.R " + wd4r + " " + projectName + " " + numberOfFamilies + " " + height); 
			
			
			if(clusteringMeasure == 1 || clusteringMeasure == 2) {

				double height2 = Double.parseDouble(params.getProperty("coreHeight", "0"));
				condition = "Groups_CoreProteins_h" + height2 + "/";

				motifInstancesPrefix = wd + "motifFamilies/" + familyFolder + condition + networkName + clusteringName + "coreProteins_h" + height2 + "_motifInstances_motifFamily";
				motifPPMPrefix = wd +  "motifFamilies/" + familyFolder + condition + networkName + clusteringName + "coreProteins_h" + height2 + "_ppm_motifFamilyGroup";
				motifInfoFile = wd +  "motifFamilies/" + familyFolder + condition + networkName + clusteringName + "coreProteins_h" + height2 + "_motifFamiliesInfo.tsv";

				motifFamilyFilePrefix = wd +  "motifFamilies/" + familyFolder + condition + networkName + clusteringName + "_CoreProteins_p" + pvalThreshold + "_ward.D2_h"+ height2 +"_group";
				int numberOfCoreFamilies = Integer.parseInt(params.getProperty("coreFamilyGroups", "10"));

				System.out.println("**Assessing motif instances for Core Proteins**");
				MotifFamily.assessMotifFamilies(motifFamilyFilePrefix, numberOfCoreFamilies, significantMotifsFile, enumeratedMotifs, proteinToRefSeqIdFile, motifInstancesPrefix, motifPPMPrefix, motifInfoFile, corePorteinsFile);
				
				wd4r = wd + "motifFamilies/" + familyFolder + condition;
				projectName = networkName + clusteringName + "_coreProteins";
				Runtime.getRuntime().exec("Rscript seqLogo.R " + wd4r + " " + projectName + " " + numberOfFamilies + " " + height); 

			}

		}

	}


}
