package benchmark;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

public class CompareProteinSetSimilarity {

	public static void main(String[] args) {
		
		String wd = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/";

		String mclClustersFile = wd + "benchmark/mclOutput_i2.txt";
		String proteinInfoFile = wd + "corrNetTop2-400_proteinsInNetwork_info.tsv";

		String motifFamilySummaryFile = wd + "corrNetTop2-400_coreTPD_p0.4_coreProteins_h0.7_motifFamiliesInfo.tsv";
		String annotationSubsetFile = wd + "corrNetTop2-400_coreTPD_p0.4_coreProteinsByMotif.tsv";

		String outputFile = wd + "benchmark/corrNetTop2-400_coreTPD_p0.4_MCL_i2_similarity_coreProts.tsv";
		//String protDistribution = "C:\\Users\\Rachel\\Documents\\LESMoNlocal\\analysis\\benchmark\\corrNetTop2-400_coreTPD_p0.4_MCL_i2_proteinDistribution_coreProts.tsv";

		/* Load MCL clusters (protein sets) as list<Set<Proteins>> */ 
		System.out.println("**Loading mcl clusters**");
		List<HashSet<String>> mclClusters = loadMCLclusters(mclClustersFile, proteinInfoFile);
		System.out.println("Loaded mcl clusters: " + mclClusters.size());

		/* Load motif family subset */ 
		System.out.println("**Loading motif family subset**");
		List<HashSet<String>> motifFamilySubset = loadMotifFamilyProteinSet(annotationSubsetFile, motifFamilySummaryFile);

		/* Iterate through annotation subset - calculate percent similarity - output to file */ 
		System.out.println("**Computing similarity**");
		computeProteinSetSimilarity(motifFamilySubset, mclClusters, outputFile);

//		System.out.println("**Assess Protein distribution**");
//		assessDistributionOfProteinsInMCLclusters(motifFamilySubset, mclClusters, protDistribution);
	}

	/**
	 * Load the protein set corresponding to each MCL cluster, line order is important
	 * 
	 * @param mclClustersFile	String - input file, each line corresponds to 1 MCL cluster 
	 * @return mclClusters		List<HashSet<String>> - List<Set {Protein 1|2|..|n} >
	 */
	private static List<HashSet<String>> loadMCLclusters(String mclClustersFile, String proteinInfoFile){
		List<HashSet<String>> mclClusters = new ArrayList<HashSet<String>>();

		HashSet<String> proteinsWithSequence = loadProteinsWithSequence(proteinInfoFile);

		InputStream in;
		try {
			in = new FileInputStream(new File(mclClustersFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // no header
			while(line!=null) {
				/* Each line corresponds to an mcl cluster (or set of proteins) */
				HashSet<String> proteinSet = new HashSet<>(Arrays.asList(line.split("\t")));
				HashSet<String> proteinSet2 = new HashSet<>();

				for(String prot: proteinSet) {
					if(proteinsWithSequence.contains(prot)) {
						proteinSet2.add(prot);
					}
				}

				mclClusters.add(proteinSet2);
				line = input.readLine();
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return mclClusters;
	}

	private static HashSet<String> loadProteinsWithSequence(String proteinInfoFile){

		HashSet<String> proteinsWithSequence = new HashSet<>();
		try {
			InputStream in = new FileInputStream(new File(proteinInfoFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();
			while(line !=null) {

				if(line.split("\t").length > 1) {  // check if protein has corresponding refSeqIds (possible there are none)
					proteinsWithSequence.add(line.split("\t")[0]);
				}
				line = input.readLine();
			}

			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return proteinsWithSequence;
	}

	/**
	 * Load representative motifs for each motif family contained in motif info file
	 * Note: order of motifs is important
	 * 
	 * @param motifInfoFile		String - input file
	 * @return motifFamilies 	Map<String, int> - Map<Motif1=0 > 0 = index
	 */
	private static HashMap<String, Integer> getMotifFamilies(String motifInfoFile){

		HashMap<String, Integer> motifFamilies = new HashMap<>();

		InputStream in;
		try {
			in = new FileInputStream(new File(motifInfoFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // header 
			line = input.readLine();
			int motifCount = 0;
			while(line!=null) {
				motifFamilies.put(line.split("\t")[1], motifCount); // [1] = motif representation XXXXXXXX
				line = input.readLine();
				motifCount++;
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return motifFamilies;
	}

	/**
	 * Load representative motifs from motif info file and obtain their corresponding proteins from the annotation subset file
	 * 
	 * @param annotationSubset	String - file motif1[XXXXX] \t protein1|protein2|...|n
	 * @param motifInfoFile 	String - file listing representative motifs
	 * 
	 * @return motifFamilyProteinSets	List<HashSet<String>> - list of protein sets corresponding to motifs families
	 */
	private static List<HashSet<String>> loadMotifFamilyProteinSet(String annotationSubset, String motifInfoFile){

		/* Get list of representative motifs */
		HashMap<String, Integer> motifFamiliesMap = getMotifFamilies(motifInfoFile);
		System.out.println("Loaded representative motifs: " + motifFamiliesMap.size() + "\n");

		/* Initialize list*/
		List<HashSet<String>> motifFamilyProteinSets = new ArrayList<>(motifFamiliesMap.size());
		HashSet<String> set = new HashSet<>();
		for(int i=0; i<motifFamiliesMap.size(); i++) {
			motifFamilyProteinSets.add(set);
		}
		InputStream in;
		try {
			in = new FileInputStream(new File(annotationSubset));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // no header
			int motifCount=0;  // keep count of motif families identified in annotationSubset

			/* iterate through annotation subset file to identify protein sets corresponding to motif families*/
			while(line!=null || motifCount < motifFamiliesMap.size()) {

				String motif = line.split("\t")[0]; 
				if(motifFamiliesMap.containsKey(motif)) {
					motifCount++;
					System.out.print(motifCount + ".");

					int idxOfMotif = motifFamiliesMap.get(motif);				
					motifFamilyProteinSets.set(idxOfMotif, new HashSet<>(Arrays.asList(line.split("\t")[2].split("\\|"))));
				}
				line=input.readLine();
			}

			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		System.out.println("Done");
		return motifFamilyProteinSets;
	}

	/**
	 * Compare protein sets of mcl clusters (size n) and motif families (size m). Print matrix of size (n x m). 
	 * Note only mcl clusters with greater than 3 proteins with a corresponding sequence will be output (therefore actual size may be smaller than n)
	 * 
	 * @param motifFamilySubset		List<HashSet<String>> - set of proteins corresponding to motif families
	 * @param mclClusters			List<HashSet<String>> - set of proteins corresponding to MCL clusters
	 * @param outputFile			String - File - similarity matrix between MCL clusters and motif families
	 */
	private static void computeProteinSetSimilarity(List<HashSet<String>> motifFamilySubset, List<HashSet<String>> mclClusters, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			/* header */
			out.write("mcl-cluster\t");
			for(int i=1; i<= motifFamilySubset.size(); i++) {
				out.write("motif" + i + "\t");
			}
			out.write("Max-Family\tMax-value\n");
			out.flush();

			/* For each MCL cluster and motifFamily, compute similarity overlap of protein set and print to file */
			for(int i=0; i<mclClusters.size(); i++) {
				HashSet<String> cluster = mclClusters.get(i);
				System.out.print(i + ".");
				
				/* Ensure MCL cluster is associated to at least 3 sequences */
				if(!cluster.isEmpty() && cluster.size() >= 3) {
					out.write((i+1) + "\t");

					/* compare MCL cluster protein set similarity against each motif family */
					List<Double> similarityList = new ArrayList<>();
					for(int j=0; j<motifFamilySubset.size(); j++) {
						HashSet<String> familySet = motifFamilySubset.get(j);

						/* count number of shared proteins */ 
						int countSharedProteins = 0;
						for(String protein: cluster) {
							if(familySet.contains(protein)) {
								countSharedProteins++;
							}
						}
						/* compute similarity */ 
						double similarity = countSharedProteins / (double) Math.min(cluster.size(), familySet.size());

						similarityList.add(similarity);
						out.write(similarity + "\t");
					}
					double max = Collections.max(similarityList);
					int idx = similarityList.indexOf(max) + 1;
					out.write(idx + "\t" + max + "\n");
					out.flush();
				}
			}
			out.close();
			System.out.println("Done");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}


	@SuppressWarnings("unused")
	private static void assessDistributionOfProteinsInMCLclusters(List<HashSet<String>> motifFamilySubset, List<HashSet<String>> mclClusters, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			/* header */
			out.write("mcl-cluster\t");
			for(int i=1; i<= motifFamilySubset.size(); i++) {
				out.write("motif" + i + "\t");
			}
			out.write("\n");
			out.flush();

			/* For each MCL cluster and motifFamily, compute similarity overlap of protein set and print to file */
			for(int i=0; i<mclClusters.size(); i++) {
				HashSet<String> cluster = mclClusters.get(i);
				out.write((i+1) + "\t");
				System.out.print(i + ".");

				List<Double> similarityList = new ArrayList<>();
				for(int j=0; j<motifFamilySubset.size(); j++) {
					HashSet<String> familySet = motifFamilySubset.get(j);

					/* count number of shared proteins */ 
					int countSharedProteins = 0;
					for(String protein: cluster) {
						if(familySet.contains(protein)) {
							countSharedProteins++;
						}
					}
					/* compute similarity */ 
					double similarity = countSharedProteins / (double) familySet.size();

					similarityList.add(similarity);
					out.write(similarity + "\t");
				}
				out.write("\n");
				out.flush();
			}
			out.close();
			System.out.println("Done");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
