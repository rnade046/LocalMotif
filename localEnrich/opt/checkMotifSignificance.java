package opt;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.HashSet;

public class checkMotifSignificance {

	public static void main(String[] args) {

		String motifsToCheckFile = "/home/rnade046/projects/rrg-mlaval/rnade046/motifsToCheck_v1.tsv";

		String motifsPrefix = "/home/rnade046/projects/rrg-mlaval/rnade046/motifs_Full_FWD/motifClustering/corrNetTop2_testedDegenMotifClustering_";
		String nullPrefix = "/home/rnade046/projects/rrg-mlaval/rnade046/motifDegen_Full_REV/motifClustering/corrNetTop2_nullModel_testedDegenMotifClustering_";
		
		String outputResults = "/home/rnade046/projects/rrg-mlaval/rnade046/Clustering_motifsToCheck_v1.tsv";
		int numFiles = 999;
		
		/* Load known motifs - Set<String> */
		System.out.println("**Loading motifs to test**");
		HashSet<String> motifsToCheck = loadMotifsToCheck(motifsToCheckFile);
		System.out.println("Loaded motifs: " + motifsToCheck.size() + "\n");
		
		/* Check motifs - Map<String, Double>*/
		System.out.println("**Searching motif clustering**");
		HashMap<String, Double[]> motifsClusteringMap = searchForMotifClustering(motifsToCheck, motifsPrefix, numFiles);
		System.out.println("Found motifs: " + motifsClusteringMap.size() + "\n");
		
		/* Check degen - Map<String, Double>*/ 
		System.out.println("**Searching null clustering**");
		HashMap<String, Double[]> nullClusteringMap = searchForMotifClustering(motifsToCheck, nullPrefix, numFiles);
		System.out.println("Found motifs: " + nullClusteringMap.size() + "\n");
		
		/* print results */
		System.out.println("**Printing results**");
		printClusteringResults(motifsToCheck, motifsClusteringMap, nullClusteringMap, outputResults);
	
	}

	private static HashSet<String> loadMotifsToCheck(String motifsToCheckFile) {

		HashSet<String> motifsToCheck = new HashSet<>();

		try {
			InputStream in = new FileInputStream(new File(motifsToCheckFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();

			while(line!=null) {

				motifsToCheck.add(line);
				line = input.readLine();
			}

			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return motifsToCheck;
	}

	private static HashMap<String, Double[]> searchForMotifClustering(HashSet<String> motifsToCheck, String motifPrefix, int numFiles) {

		HashMap<String, Double[]> mapOfMotifClustering = new HashMap<>();

		/* iterate over all degen motif annotation files */
		for(int i=0; i < numFiles; i++) {

			if(i%10==0) {
				System.out.print(i + ".");
			}

			if(i%100==0) {
				System.out.println();
			}

			String motifFile = motifPrefix + i; // current annotation file (ex. mapDegenMotifsToProteins_200) 

			InputStream in;
			try {
				in = new FileInputStream(new File(motifFile));
				BufferedReader input = new BufferedReader(new InputStreamReader(in));

				String line = input.readLine(); // no header

				while(line!=null) {

					String motif = line.split("\t")[0];

					if(motifsToCheck.contains(motif)) {
						double nprot = Double.parseDouble(line.split("\t")[1]);
						double tpd = Double.parseDouble(line.split("\t")[2]);
						double pval = Double.parseDouble(line.split("\t")[3]);
						Double[] array = {nprot, tpd, pval};
						
						mapOfMotifClustering.put(motif, array);

						if(mapOfMotifClustering.size() ==  motifsToCheck.size()) {
							break;
						}
					}
					line = input.readLine();
				}

				input.close();
			} catch(IOException e) {
				e.printStackTrace();
			}
		}
		System.out.println("Done\n");
		return mapOfMotifClustering;
	}
	
	private static void printClusteringResults(HashSet<String> motifsToCheck, HashMap<String, Double[]> motifsClusteringMap, 
			HashMap<String, Double[]> nullClusteringMap, String outputFile) {
	
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			
			out.write("Motif\tnProt\tTPD\tpval\tNullnProt\tNullTPD\tNullpval\n");
			
			for(String motif: motifsToCheck) {
				
				if(motifsClusteringMap.containsKey(motif) && nullClusteringMap.containsKey(motif)) {

					out.write(motif + "\t" + motifsClusteringMap.get(motif)[0].intValue() + "\t" + motifsClusteringMap.get(motif)[1] + "\t" + motifsClusteringMap.get(motif)[2] 
					+ "\t" + nullClusteringMap.get(motif)[0].intValue() + "\t"+ nullClusteringMap.get(motif)[1] + "\t"+ nullClusteringMap.get(motif)[2]+ "\n");
				}
				
			}
			
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}


}
