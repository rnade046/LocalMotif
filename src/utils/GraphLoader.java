package utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import graph.Interaction;

public class GraphLoader {

	public static ArrayList<Interaction> loadInteractionRepository(String repositoryFile, String baitMapFile, String preyMapFile, double fdr){

		HashMap<String, Double> interactionToSpectralCountMap = loadSaintReport(repositoryFile, fdr);
		HashMap<String, ArrayList<ArrayList<String>>> baitMap = loadBaitMapping(baitMapFile);
		HashMap<String, ArrayList<ArrayList<String>>> preyMap = loadPreyMapping(preyMapFile);
		
	
		//printProteinsInNetwork(interactionToSpectralCountMap);
		
		ArrayList<Interaction> interactionList = formatInteractionList(interactionToSpectralCountMap, baitMap, preyMap);
		
		return interactionList;
	}
	
	/**
	 * Loads the significant PPIs from the saint express report from cell-map.org as a mapping of the 
	 * interaction and the average spectral counts. PPIs are deemed significant if their BFDR is <= to
	 * the FDR cut-off outlined by the user. 
	 * 
	 * @param inputRepositoryFile	String - file path to SaintExpress report
	 * @param fdr					double - FDR significant cut off
	 * 
	 * @return interactionToSpectralCountMap	HashMap<String, Double> {bait \t prey} : {average spectral count} 
	 */
	private static HashMap<String,Double> loadSaintReport(String inputRepositoryFile, double fdr) {

		HashMap<String, Double> interactionToSpectralCountMap = new HashMap<String, Double>();// To contain interaction name and ID, and number of occurence
		try {
			/* Read BioGrid ; get all possible human protein-protein interactions */

			InputStream in = new FileInputStream(new File(inputRepositoryFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // header
			line = input.readLine();
			
			while (line != null) { // stops when there are no more lines in file (in)

				String[] col = line.split("\t"); // split line by tabs; obtain individual columns

				if(Double.parseDouble(col[15]) <= fdr) { // col[15] = BFDR
					/* Keep interactions that pass FDR thresholds */ 
					String interactor1 = col[0]; // get name of bait 
					String interactor2 = col[1].split("\\.")[0]; // get RefSeq Protein Accession of prey & remove .version
					
					interactionToSpectralCountMap.put(interactor1+"\t"+interactor2, Double.parseDouble(col[5])); //col[5] = Average Spectral count
				} 
				line = input.readLine(); // read next line
			}
			input.close(); // close BufferedReader
		} catch (IOException e) {
			e.printStackTrace();
		}
		return interactionToSpectralCountMap;
	}
	
	/**
	 * Load bait mapping file which contains various accession IDs for the baits used to generate de human cell map. 
	 * The loaded file is an accessory file on cell-map.org. The information is stored in a HashMap and maps the 
	 * bait name as decided by the Gingras lab to its RefSeq mRNA accession and its conventional gene name (HGNC symbol).
	 * 
	 * Note: allow multiple RefSeq mRNA accessions but will only have one gene symbol
	 * 
	 * @param inputFile	String - file path to bait map
	 * 
	 * @return baitMap	HashMap<String, ArrayList<String>> {bait name} : {RefSeq mRNA accession List | gene symbol List}
	 */
	private static HashMap<String, ArrayList<ArrayList<String>>> loadBaitMapping(String inputFile){
		HashMap<String, ArrayList<ArrayList<String>>> baitMap = new HashMap<>();
		
		InputStream in;
		try {
			in = new FileInputStream(new File(inputFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));
			
			String line = input.readLine(); //header
			line = input.readLine();
			
			while(line != null) {
				String[] col = line.split("\t");
				
				/* update RefSeqIDs if bait has already been seen*/ 
				if(baitMap.containsKey(col[1])) {
					ArrayList<ArrayList<String>> map = baitMap.get(col[1]);
					ArrayList<String> refSeq = map.get(0);
					refSeq.add(col[9]);
				/* create new mapping if bait hasn't been seen */	
				} else { 
					ArrayList<ArrayList<String>> map = new ArrayList<ArrayList<String>>(); // map: RefSeq mRNA accession | gene symbol
					ArrayList<String> refSeq = new ArrayList<>();
					ArrayList<String> symbol = new ArrayList<>();
					
					if(col.length > 9 && col[9] != null) {
						refSeq.add(col[9]); //RefSeq mRNA accession
					}
					symbol.add(col[3]); //gene symbol
					
					map.add(0, refSeq);
					map.add(1, symbol);
					
					baitMap.put(col[1], map); // col[1] = bait name as seen in 1st column of Saint Express file
				} 
				line = input.readLine();
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return baitMap;
	}
	
	/**
	 * Load prey mapping file which contains accessions corresponding to preys in the Saint Express report. This mapping file was
	 * generated using the biomaRt package in R. The information is stored in a HashMap and maps the RefSeq protein accession to 
	 * the RefSeq mRNA accession and it's conventional gene name (HGNC symbol).
	 * 
	 * Note: allow multiple RefSeq mRNA accessions but will only have one gene symbol
	 * 
	 * @param inputFile	String - file path to prey map
	 * 
	 * @return preyMap	HashMap<String, String[]> - map {RefSeq protein Accession} : [RefSeq mRNA Accession, HGNC symbol]
	 */
	private static HashMap<String,ArrayList<ArrayList<String>>> loadPreyMapping(String inputFile){
		HashMap<String, ArrayList<ArrayList<String>>> preyMap = new HashMap<>();
		
		InputStream in;
		try {
			in = new FileInputStream(new File(inputFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));
			
			String line = input.readLine(); //header
			line = input.readLine();
			
			while(line != null) {
				String[] col = line.split("\t");
				
				/* update RefSeqIDs if bait has already been seen*/ 
				if(preyMap.containsKey(col[0])) {
					ArrayList<ArrayList<String>> map = preyMap.get(col[0]);
					ArrayList<String> refSeq = map.get(0);
					refSeq.add(col[1]);
				/* create new mapping if bait hasn't been seen */	
				} else { 
					ArrayList<ArrayList<String>> map = new ArrayList<ArrayList<String>>(); // map: RefSeq mRNA accession | gene symbol
					ArrayList<String> refSeq = new ArrayList<>();
					ArrayList<String> symbol = new ArrayList<>();
					
					refSeq.add(col[1]); //RefSeq mRNA accession
					symbol.add(col[2]); //gene symbol
					
					map.add(0, refSeq);
					map.add(1, symbol);
					
					preyMap.put(col[0], map); // col[0] = RefSeq protein accession which maps to info in SaintExpress report 
				} 
				
				line = input.readLine();
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
		return preyMap;
	}
	
	private static ArrayList<Interaction> formatInteractionList(HashMap<String, Double> interactionToSpectralCountMap, HashMap<String, ArrayList<ArrayList<String>>> baitMap, HashMap<String, ArrayList<ArrayList<String>>> preyMap){
		
		ArrayList<Interaction> interactionList = new ArrayList<>(); // list to contain formated interactions
		
		/* sets to contain missing mappings */ 
		HashSet<String> missingBaitMapping = new HashSet<>();
		HashSet<String> missingPreyMapping = new HashSet<>();
		
		for(String interaction: interactionToSpectralCountMap.keySet()) {
			
			/* bait and prey identifiers for significant saint PPI */
			String bait_name = interaction.split("\t")[0];
			String prey_RefSeqProtAcc = interaction.split("\t")[1];
			
			/* get accessions for bait [RefSeq mRNA Acc, HGNC symbol] */ 
			ArrayList<ArrayList<String>> baitAcc = null;
			if(baitMap.containsKey(bait_name)) {
				baitAcc = baitMap.get(bait_name);
			} else { 
				missingBaitMapping.add(bait_name);
			}
			
			/* get accessions for prey [RefSeq mRNA Acc, HGNC symbol] */ 
			ArrayList<ArrayList<String>> preyAcc = null;
			if(preyMap.containsKey(prey_RefSeqProtAcc)) {
				preyAcc = preyMap.get(prey_RefSeqProtAcc);
			} else { 
				missingPreyMapping.add(prey_RefSeqProtAcc);
			}
			
			/* Interaction constructor: 
			 * Interaction(String _Protein1, String _Protein2, String _ID1, String _ID2, double _w)
			 * Interaction(Bait-gene-symbol, Prey-gene-symbol, Bait-RefSeq-mRNA-Acc, Prey-RefSeq-mRNA-Acc, weight)*/
			Interaction ppi = null;
			if (preyAcc == null) {
				ppi = new Interaction(baitAcc.get(1).get(0), prey_RefSeqProtAcc, baitAcc.get(0), null, 1/interactionToSpectralCountMap.get(interaction));
			} else { 
				ppi = new Interaction(baitAcc.get(1).get(0), preyAcc.get(1).get(0), baitAcc.get(0), preyAcc.get(0), 1/interactionToSpectralCountMap.get(interaction));
			}
			interactionList.add(ppi);
		}
		
		if(missingBaitMapping.size() != 0) {
			System.out.print("Missing baits mappings (" + missingBaitMapping.size() +") : ");
			for(String bait: missingBaitMapping) {
				System.out.print(bait + "|");
			}
			System.out.println();
		}
		
		if(missingPreyMapping.size() != 0) {
			System.out.print("Missing prey mappings (" + missingPreyMapping.size() +") : ");
			for(String prey: missingPreyMapping) {
				System.out.print(prey + "|");
			}
			System.out.println("\n");
			
		}
		return interactionList;
	}
	
	/**
	 * Function to print information related to the cell map network 
	 * @param interactionToSpectralCountMap
	 */
	@SuppressWarnings("unused")
	private static void printProteinsInNetwork(HashMap<String, Double> interactionToSpectralCountMap) {
	
		HashSet<String> proteinSet = new HashSet<>();
		
		for(String interaction : interactionToSpectralCountMap.keySet()) {
			
			String[] interactors = interaction.split("\t");
			proteinSet.add(interactors[0]);
			//proteinSet.add(interactors[1]);
		}
		
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File("C:\\Users\\Rachel\\Documents\\LESMoNlocal\\baits_refseqIDs_FDR0.01.tsv")));
		
			for(String protein: proteinSet) {
				//Write the RefSeq Protein accession
				//String[] prot = protein.split("\\.");
				//out.write(prot[0] + "\n");
				
				//Wrtie protein name
				out.write(protein + "\n");
				
				out.flush();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
}
