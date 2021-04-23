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

public class SaintGraphLoader {

	public static ArrayList<Interaction> loadInteractionRepositoryFromSaintExpressReport(String repositoryFile, String cellMapBaitMapFile, String bioMartBaitMapFile, String bioMartPreyMapFile, double fdr){

		HashMap<String, Double> interactionToSpectralCountMap = loadSaintReport(repositoryFile, fdr);
		HashMap<String,String> cellMapBaitMap = loadCellMapBaitMapping(cellMapBaitMapFile);
		HashMap<String,String> cellMapPreyMap = loadCellMapPreytMapping(repositoryFile, fdr);

		HashMap<String, ArrayList<String>> bioMartBaitMap = loadBioMartMapping(bioMartBaitMapFile);
		HashMap<String, ArrayList<String>> bioMartPreyMap = loadBioMartMapping(bioMartPreyMapFile);


		//printProteinsInNetwork(cellMapPreyMap);

		ArrayList<Interaction> interactionList = formatInteractionList(interactionToSpectralCountMap, cellMapBaitMap, bioMartBaitMap, cellMapPreyMap, bioMartPreyMap);

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
					String interactor1 = col[0]; // get name of bait (SAINT report bait name)
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
	 * Load prey mapping file which contains various accession IDs for the baits used to generate the human cell map. 
	 * The loaded file is an accessory file on cell-map.org. The information is stored in a HashMap and maps the 
	 * bait name as decided by the Gingras lab to its conventional gene name (HGNC symbol).
	 * 
	 * @param inputFile	String - file path to bait map
	 * 
	 * @return baitMap	HashMap<String, String> {proteinRefSeq} : {gene symbol List}
	 */
	private static HashMap<String, String> loadCellMapPreytMapping(String inputFile, double fdr){
		HashMap<String, String> preyMap = new HashMap<>();

		InputStream in;
		try {
			in = new FileInputStream(new File(inputFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); //header
			line = input.readLine();

			while(line != null) {
				String[] col = line.split("\t");
				if(Double.parseDouble(col[15]) <= fdr) { // col[15] = BFDR
					preyMap.put(col[1].split("\\.")[0], col[2]); // col[1] = prey RefSeqProtein ID
				}												 // col[2] = prey HGNC symbol
				line = input.readLine();
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return preyMap;
	}

	/**
	 * Load bait mapping file which contains various accession IDs for the baits used to generate the human cell map. 
	 * The loaded file is an accessory file on cell-map.org. The information is stored in a HashMap and maps the 
	 * bait name as decided by the Gingras lab to its conventional gene name (HGNC symbol).
	 * 
	 * @param inputFile	String - file path to bait map
	 * 
	 * @return baitMap	HashMap<String, String> {bait name} : {gene symbol List}
	 */
	private static HashMap<String, String> loadCellMapBaitMapping(String inputFile){
		HashMap<String, String> baitMap = new HashMap<>();

		InputStream in;
		try {
			in = new FileInputStream(new File(inputFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); //header
			line = input.readLine();

			while(line != null) {
				String[] col = line.split("\t");

				baitMap.put(col[1], col[3]); // col[1] = bait name as seen in 1st column of Saint Express file
				// col[3] = baits HGNC symbol
				line = input.readLine();
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return baitMap;
	}

	/**
	 * Load bait mapping file generated using the R package BioMart. It contains the RefSeq mRNA ID for a given HGNC symbol
	 * 
	 * @param inputFile 	String - path to BioMart Bait mapping file
	 * @return baitMap		HashMap<String, ArrayList<String>> - {HGNC symbol : List of RefSeq mRNA IDs} 
	 */
	private static HashMap<String, ArrayList<String>> loadBioMartMapping(String inputFile){
		HashMap<String, ArrayList<String>> baitMap = new HashMap<>();

		InputStream in;
		try {
			in = new FileInputStream(new File(inputFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); //no header
			line = input.readLine();

			while(line != null) {
				String[] col = line.split("\t");

				/* if map contains HGNC symbol, append new MRNA Seq ID to exisitng list */
				if(baitMap.containsKey(col[0])) { // col[0] = baits HGNC symbol
					ArrayList<String> refSeqList = baitMap.get(col[0]);
					refSeqList.add(col[1]);
					/* otherwise create new mapping */ 
				}else { 
					ArrayList<String> refSeqList = new ArrayList<>();
					refSeqList.add(col[1]);		 //col[1] = RefSeq mRNA ID
					baitMap.put(col[0], refSeqList); 
				}

				line = input.readLine();
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return baitMap;
	}

	private static ArrayList<Interaction> formatInteractionList(HashMap<String, Double> interactionToSpectralCountMap, 
			HashMap<String, String> cellMapBaitMap, HashMap<String, ArrayList<String>> bioMartBaitMap, HashMap<String, String> cellMapPreyMap,
			HashMap<String, ArrayList<String>> bioMartPreyMap){

		ArrayList<Interaction> interactionList = new ArrayList<>(); // list to contain formated interactions

		/* sets to contain missing mappings */ 
		HashSet<String> missingBaitMapping = new HashSet<>();
		HashSet<String> missingPreyMapping = new HashSet<>();

		for(String interaction: interactionToSpectralCountMap.keySet()) {

			/* bait and prey identifiers for significant saint PPI */
			String bait_name = interaction.split("\t")[0];
			String prey_RefSeqProtAcc = interaction.split("\t")[1];

			/* Identify bait HGNC symbol and it's RefSeq mRNA IDs  */ 
			String baitHGNCsymbol = null;
			ArrayList<String> baitRefSeqIds = null;
			if(cellMapBaitMap.containsKey(bait_name)) {
				baitHGNCsymbol = cellMapBaitMap.get(bait_name);
				if(bioMartBaitMap.containsKey(baitHGNCsymbol)) {
					baitRefSeqIds = bioMartBaitMap.get(baitHGNCsymbol);	
				} else {
					missingBaitMapping.add(baitHGNCsymbol);
				}
			} else { 
				missingBaitMapping.add(bait_name);
			}

			/* get RefSeq mRNA ids */ 
			ArrayList<String> preyRefSeqRNAids = null;
			String preyHGNCsymbol = cellMapPreyMap.get(prey_RefSeqProtAcc);
			if(bioMartPreyMap.containsKey(prey_RefSeqProtAcc)) {
				preyRefSeqRNAids = bioMartPreyMap.get(prey_RefSeqProtAcc);
			} else { 
				missingPreyMapping.add(prey_RefSeqProtAcc);
			}

			/* Interaction constructor: 
			 * Interaction(String _Protein1, String _Protein2, List<String> _ID1, List<String> _ID2, double _w)
			 * Interaction(Bait-gene-symbol, Prey-gene-symbol, Bait-RefSeq-mRNA-Acc, Prey-RefSeq-mRNA-Acc, weight)*/
			Interaction ppi = null;
			if (preyRefSeqRNAids == null) {
				ppi = new Interaction(baitHGNCsymbol, preyHGNCsymbol, baitRefSeqIds, null, 1/interactionToSpectralCountMap.get(interaction));
			} else { 
				ppi = new Interaction(baitHGNCsymbol, preyHGNCsymbol, baitRefSeqIds, preyRefSeqRNAids, 1/interactionToSpectralCountMap.get(interaction));
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
	private static void printProteinsInNetwork(HashMap<String,String> interactionToSpectralCountMap) {

		HashSet<String> proteinSet = new HashSet<>();

		/*for(String interaction : interactionToSpectralCountMap.keySet()) {

			String[] interactors = interaction.split("\t");
			proteinSet.add(interactors[0]);
			//proteinSet.add(interactors[1]);
		}*/

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File("C:\\Users\\Rachel\\Documents\\LESMoNlocal\\baits_HGNCsymbols.tsv")));

			for(String protein: interactionToSpectralCountMap.keySet()) {
				//Write the RefSeq Protein accession
				//String[] prot = protein.split("\\.");
				//out.write(prot[0] + "\n");

				//Wrtie protein name
				out.write(interactionToSpectralCountMap.get(protein) + "\n");

				out.flush();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
}
