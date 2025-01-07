package overconnectedProteins;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;

public class DetermineLocalizationOfOverConnectedProteins {

	public static void main(String[] args) {

		String originalNetworkFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/corrNetTop2_protFreqAnnotation.tsv";
		String finalNetworkFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/corrNetTop2-400_protFreqAnnotation.tsv";
		
		String localizationFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/localization/preys-latest.txt";
		
		String outputFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/localization/corrNet2-400_overConnectedProteins_localizations.tsv";
		String freqOutFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/localization/corrNet2-400_overConnectedProteins_FrequenciesOfLocalizations.tsv";
		
		/* get list of removed proteins */
		HashSet<String> overConnectedProteins = determineListOfRemovedProteins(originalNetworkFile, finalNetworkFile);
		
		/* get their localizations */
		List<Protein> proteinList = loadPreyInfo(localizationFile, overConnectedProteins);
		
		/* output results */
		printProteinLocalization(outputFile, proteinList);
		locationFrequencyMap(freqOutFile, proteinList);
	}

	
	/***
	 * Compare original list of proteins in network (corrNet2) and final list (corrNet2-400) 
	 * to obtain list of proteins that were deemed overly connected.
	 * 
	 * @param originalNetworkProteinFile	String - file 
	 * @param finalNetworkProteinFile		String - file
	 * @return overConnected Proteins 		List<Protein> 
	 */
	private static HashSet<String> determineListOfRemovedProteins(String originalNetworkProteinFile, String finalNetworkProteinFile){

		HashSet<String> overConnectedProteins = new HashSet<>();

		/* get list of proteins in disconnected network */
		HashSet<String> finalProteinSet = new HashSet<>();
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(finalNetworkProteinFile))));

			String line = in.readLine();
			while(line!=null) {
				
				finalProteinSet.add(line.split("\t")[0]);
				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		/* while inputting original network - add proteins to list that are not in final network */ 
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(originalNetworkProteinFile))));

			String line = in.readLine();
			while(line!=null) {
				
				String prot = line.split("\t")[0];
				if(!finalProteinSet.contains(prot)) {
					overConnectedProteins.add(prot);
				}
				line = in.readLine();
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return overConnectedProteins;
	}

	private static List<Protein> loadPreyInfo(String preyFile, HashSet<String> proteinSet){ 
		
		List<Protein> proteinList = new ArrayList<>();
		
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(preyFile))));
			
			String line = in.readLine(); // header 
			line = in.readLine();
			
			while(line!=null) {
				
				String[] col = line.split("\t");
				
				if(proteinSet.contains(col[0])) {
					
					int safe = 25;
					int nmf = 19;
					
					if(!col[3].equals("-")) {
						safe = Integer.parseInt(col[3]);
					}
					
					if(!col[1].equals("19")) {
						nmf = Integer.parseInt(col[1]);
					}
					
					proteinList.add(new Protein(col[0], safe, nmf));
				}
				
				line = in.readLine();
			}
			
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	
		return proteinList;
	}
	
	private static void printProteinLocalization(String outputFile, List<Protein> proteinList) {
		
		BufferedWriter out;
		try {
			out = new BufferedWriter(new FileWriter(new File(outputFile)));
			
			out.write("Protein\tNMF\tSafe\n");
			
			for(Protein p :proteinList) {
				
				String[] info = p.getProteinInfo();
				out.write(info[0] + "\t" + info[1] + "\t" + info[2] + "\n");
				out.flush();
				
			}
			
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	private static void locationFrequencyMap(String outputFile, List<Protein> proteinList) {
		
		HashMap<String, Integer> nmfFreqMap = new HashMap<>();
		HashMap<String, Integer> safeFreqMap = new HashMap<>();
		
		/* count occurrences of localizations */
		for(Protein p: proteinList) {
			String[] info = p.getProteinInfo();
			
			// NMF
			if(nmfFreqMap.containsKey(info[1])) {
				nmfFreqMap.put(info[1], nmfFreqMap.get(info[1])+1);
			} else {
				nmfFreqMap.put(info[1], 1);
			}
			
			// SAFE
			if(safeFreqMap.containsKey(info[2])) {
				safeFreqMap.put(info[2], safeFreqMap.get(info[2])+1);
			} else {
				safeFreqMap.put(info[2], 1);
			}
		
		}
		
		/* print frequency maps */
		BufferedWriter out;
		try {
			out = new BufferedWriter(new FileWriter(new File(outputFile)));
			
			out.write("NMF\tCount\n");
			
			for(Entry<String, Integer> e: nmfFreqMap.entrySet()) {
				out.write(e.getKey() + "\t" + e.getValue() + "\n");
				out.flush();
			}
			
			out.write("\n\n\n");
			
			out.write("SAFE\tCount\n");
			
			for(Entry<String, Integer> e: safeFreqMap.entrySet()) {
				out.write(e.getKey() + "\t" + e.getValue() + "\n");
				out.flush();
			}
			
			
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
