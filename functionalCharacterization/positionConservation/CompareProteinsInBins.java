package positionConservation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class CompareProteinsInBins {


	public static void main(String[] args) {

		String proteinPerBinFilePrefix = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/MotifPosition/coreTPD0.4/prots-per-bin/corrNet2-400_coreTPD_p0.4_3UTR_FWD_proteinsPerBin__motif";
		int numberOfMotifs = 31;
		double threshold = 0.02;

		String proteinOverlapOutputFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/MotifPosition/coreTPD0.4/prots-per-bin/corrNet2-400_coreTPD_p0.4_3UTR_FWD_proteinOverlapPerBin_";
		String summaryProteinOverlapFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/MotifPosition/coreTPD0.4/prots-per-bin/corrNet2-400_coreTPD_p0.4_3UTR_FWD_proteinOverlapSummary_";
		String summaryFreqsOutputFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/MotifPosition/coreTPD0.4/prots-per-bin/corrNet2-400_coreTPD_p0.4_3UTR_FWD_summaryOfFrequencies.tsv";

		System.out.println("*Compare protein overlap*");
		compareProteinOverlap(proteinPerBinFilePrefix, numberOfMotifs, threshold, proteinOverlapOutputFile);
		
		System.out.println("*\nSummarize protein overlap*");
		summarizeOverlapOfProteins(proteinOverlapOutputFile, summaryProteinOverlapFile);

		System.out.println("\n*Summarize protein frequencies*");
		summarizeBinFrequencies(proteinPerBinFilePrefix, numberOfMotifs, summaryFreqsOutputFile);
	}

	/**
	 * Compare protein overlap within a position bin (e.g. 0-125) that contribute to a given frequency threshold (e.g. 0.02)
	 */
	private static void compareProteinOverlap(String protsPerBinFilePrefix, int numberOfMotifs, double threshold, String outputFile) {

		/* Make a list (1 entry per bin); where each element is a map {protein = number of occurrences} */
		List<HashMap<String, Integer>> proteinsPerBin = initializeMap();

		for(int i=0; i<numberOfMotifs; i++) {
			System.out.print(i + ".");
			/* Load proteins per bin file */
			try {
				BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(protsPerBinFilePrefix + (i+1) + ".tsv"))));
				
				String line = in.readLine(); // header
				line = in.readLine();
				
				int lineCount = 0; 
				
				while(line!=null) {
					
					String[] col = line.split("\t"); // [0] = motif#, [1] = maxFreq, [2] = numberProteins, [3] = list of proteins
					
					/* if current bin passes threshold, count proteins within that bin */
					if(Double.parseDouble(col[1]) >= threshold) {
						
						for(String prot: col[3].split("\\|")) {
							if(proteinsPerBin.get(lineCount).containsKey(prot)) {
								int currentCount = proteinsPerBin.get(lineCount).get(prot);
								proteinsPerBin.get(lineCount).put(prot, currentCount+1); 
							} else { 
								proteinsPerBin.get(lineCount).put(prot, 1);
							}
						}
					}
					
					line = in.readLine();
					lineCount++;
				}
				
				in.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
			
			printProteinOverlapPerBin(proteinsPerBin, outputFile);
		}

	}

	private static List<HashMap<String, Integer>> initializeMap(){

		List<HashMap<String, Integer>> proteinsPerBin = new ArrayList<>();

		for(int i=0; i<1000; i=i+125) {
			HashMap<String, Integer> proteinCount = new HashMap<>();
			proteinsPerBin.add(proteinCount);
		}

		return proteinsPerBin;
	}
	
	private static void printProteinOverlapPerBin(List<HashMap<String, Integer>> proteinsPerBin, String outputFile){
	
		for(int i=0; i<proteinsPerBin.size(); i++) {
			
			try {
				int binStart = i*125;
				int binEnd = binStart + 125;
				
				BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile + binStart + "_" + binEnd + ".tsv")));
				
				HashMap<String, Integer> proteinCount = proteinsPerBin.get(i);
				for(Map.Entry<String, Integer> e : proteinCount.entrySet()) {
					out.write(e.getKey() + "\t" + e.getValue() + "\n");
					out.flush();
				}
				out.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	
	private static void summarizeOverlapOfProteins(String proteinOverlapFile, String summaryOverlapFile) {
		
		/* load file as hash map {#contributing motifs = list of protein} */
		for(int i=0; i<1000; i=i+125) {
			int binStart = i;
			int binEnd = binStart + 125;
			
			HashMap<Integer, List<String>> proteinOverlapMap = new HashMap<>();
			
			try {
				BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(proteinOverlapFile + binStart + "_" + binEnd + ".tsv"))));
				
				String line = in.readLine(); // no header
				while(line != null) {
					
					String protein = line.split("\t")[0];
					Integer involvement = Integer.parseInt(line.split("\t")[1]); // [0] = protein name, [1] = number of motifs protein is involved in
					
					if(proteinOverlapMap.containsKey(involvement)) {
						List<String> currentList = proteinOverlapMap.get(involvement);
						currentList.add(protein);
						proteinOverlapMap.put(involvement, currentList);
					} else {
						List<String> currentList = new ArrayList<>();
						currentList.add(protein);
						proteinOverlapMap.put(involvement, currentList);
					}
					
					line = in.readLine();
				}
				
				in.close();
				
				BufferedWriter out = new BufferedWriter(new FileWriter(new File(summaryOverlapFile + binStart + "_" + binEnd + ".tsv")));
				
				out.write("#motifsWithProteinOverlap\t#Proteins\tListOfProteins\n");
				for(Map.Entry<Integer, List<String>> entry : proteinOverlapMap.entrySet()) {
					out.write(entry.getKey() + "\t" + entry.getValue().size() + "\t");
					
					for(String prot : entry.getValue()) {
						out.write(prot + "|");
					}
					
					out.write("\n");
					out.flush();
				}
				
				out.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		
	}
	
	private static void summarizeBinFrequencies(String protsPerBinPrefixFile, int numberOfMotifs, String summaryOfBinOutputFile) {
		
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(summaryOfBinOutputFile)));
			
			/* header */
			out.write("motifs\t");
			for(int i=0; i<1000; i=i+125) {
				out.write(i + "-" + (i+125) + "\t");
			}
			out.write("\n");
			out.flush();
			
			/* table summary */ 
			for(int i=0; i<numberOfMotifs; i++) {
				BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(protsPerBinPrefixFile + (i+1) + ".tsv"))));
				
				String line = in.readLine(); // header
				line = in.readLine();
				
				out.write(i+1 + "\t");
				
				while(line != null) {
					
					out.write(line.split("\t")[1] + "\t");
					
					line = in.readLine();
				}
				
				out.write("\n");
				in.close();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
