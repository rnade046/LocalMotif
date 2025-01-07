

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;

public class SummarizeTomTomResults {

	public static void main(String[] args) {
		String tomtomPrefix = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/benchmark/fasta-coreTPD/coreTPDBenchmark2/tomtom-anr-cluster";
		String outputFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/benchmark/fasta-coreTPD/coreTPD-anr-detailed-summary.tsv";
		String summaryFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/benchmark/fasta-coreTPD/coreTPD-anr-summary.tsv";
		
		summarizeResults(tomtomPrefix, outputFile, summaryFile);
	}

	public static void summarizeResults(String tomtomPrefix, String outputFile, String summaryFile) {
		
		HashMap<Integer, int[]> summarizeMap = new HashMap<>();
		HashSet<String> uniqueMotifs = new HashSet<>();
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			out.write("MCLcluster\tMEMEmotif\tMotifFamily\tFDR\n");
			
			for(int i=1; i<107; i++) {

				/* check if file exists */
				File f = new File(tomtomPrefix + i + ".tsv");
				if(f.exists() && !f.isDirectory()) { 

					System.out.println(i);

					/* store results */
					HashMap<String, String> motifMap = loadResultsForCluster(f);
					
					int significantMotifs =0;
					for(Entry<String, String> e : motifMap.entrySet()) {
						
						/* print results */
						String[] v =  e.getValue().split("\\-");
						out.write(i + "\t" + e.getKey() + "\t" + v[0] + "\t" + v[1] + "\n");
						out.flush();
					
						if(Double.parseDouble(v[1]) < 0.1) {
							significantMotifs++;
							uniqueMotifs.add(v[0]);
						}
					
						
					}
					summarizeMap.put(i, new int[] {motifMap.size(), significantMotifs});
				}
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		printSummary(summaryFile, summarizeMap, uniqueMotifs);
	}

	private static HashMap<String, String> loadResultsForCluster(File file){

		HashMap<String, String> motifMap = new HashMap<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(file)));

			String line = in.readLine(); // header 
			line = in.readLine();

			while(line != null) {

				if(!line.isEmpty() && !line.startsWith("#")) {
					String[] col = line.split("\t");
					/* if smallest FDR for MEME motif; store updated value and associated family motif */
					if(motifMap.containsKey(col[1])) {
						String[] values = motifMap.get(col[1]).split("\\-"); // [0] = family motif , [1] = FDR
						if(Double.parseDouble(values[1]) > Double.parseDouble(col[5])) {
							motifMap.replace(col[1], col[0] + "-" + col[5]);
						}
					} else {
						motifMap.put(col[1], col[0] + "-" + col[5]);
					}
				}
				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return motifMap;
	}
	
	private static void printSummary(String summaryFile, HashMap<Integer, int[]> summaryMap, HashSet<String> uniqueMotifs) {
		 
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(summaryFile)));
			out.write("MCLcluster\ttotalMotifs\tsignificantMotifs\n");
			
			for(Entry<Integer, int[]> s: summaryMap.entrySet()) {
				out.write(s.getKey() + "\t" + s.getValue()[0] + "\t" + s.getValue()[1] + "\n");
				out.flush();
			}
			
			out.write("\nUnique Motifs: " + uniqueMotifs.size() + "\n");
			for(String m: uniqueMotifs) {
				out.write(m + "\n");
				out.flush();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
}
