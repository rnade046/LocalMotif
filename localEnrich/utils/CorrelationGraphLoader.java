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
import java.util.Map.Entry;

import graph.Interaction;


public class CorrelationGraphLoader {

	public static ArrayList<Interaction> loadGraphFromCorrelationNetwork(String inputRepository, String fastaFile, String mapProtToRefSeqIdsFile,
			String proteinsInNetworkFile, double threshold, boolean removeOverlyConnectedProteins, int maxInteractions) {
		
		/* Load interactions with correlation score greater or equal to correlation threshold */
		HashMap<String, Double> confidentInteractionsMap = getConfidentInteractions(inputRepository, threshold);
		//HashSet<String> confidentProteinSet = getConfidentProteins(confidentInteractionsMap);
		//printProteinsInNetwork(confidentProteinSet, proteinsInNetworkFile);
		
		/* Load RefSeq IDs for which we have a correspond FASTA sequence */
		HashSet<String> possibleRefSeqIds = generateRefSeqSet(fastaFile);
		/* Load list of Proteins and their possible RefSeq IDs keeping only the IDs for which we have a sequence */
		HashMap<String, ArrayList<String>> mapProtToRefSeqIds = getRefSeqIdsInNetwork(mapProtToRefSeqIdsFile, possibleRefSeqIds);
		
		ArrayList<Interaction> confidentInteractions;
		
		if(removeOverlyConnectedProteins) {
			/* Check number of interactions each protein is involved in; return Set of proteins to remove */
			HashSet<String> proteinsToRemove = getOverlyConnectedProteins(confidentInteractionsMap, maxInteractions);
			System.out.println("Number of proteins with >= " + maxInteractions + ": " + proteinsToRemove.size());
			confidentInteractions = formatInteractionListWithoutOverelyConnectedProteins(confidentInteractionsMap, mapProtToRefSeqIds, proteinsToRemove);
		} else {
			confidentInteractions = formatInteractionList(confidentInteractionsMap, mapProtToRefSeqIds);
		}
		
		return confidentInteractions;
	}
	
	
	/**
	 * Load cell map correlation network, keeping interactions who's correlation scores is greater or equal to the designated threshold
	 * 
	 * @param inputRepository	String - file path to correlation network 
	 * @param threshold			double - correlation threshold required to be a "confident interaction" 
	 * 
	 * @return condifentInteractions	HashMap<String, Double> - protein1_protein2 = correlation score 
	 */
	private static HashMap<String, Double> getConfidentInteractions(String inputRepository, double threshold) {
		HashMap<String, Double> confidentInteractions = new HashMap<>();
		
		try {
			InputStream in = new FileInputStream(new File(inputRepository));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));
			
			String line = input.readLine(); // header
			line = input.readLine();
			while(line!= null) {
				
				String[] col = line.split("\t");
				
				if(Double.parseDouble(col[2]) >= threshold) {
					
					/* format of possible interactions col[0] = protein1, col[1] = protein2*/
					String interactionOpt1 = col[0] + "_" + col[1];
					String interactionOpt2 = col[1] + "_" + col[0];
					
					/* if this pair of interactors hasn't been seen previously, they are added as a new entery in the map */
					if(!confidentInteractions.containsKey(interactionOpt1) || !confidentInteractions.containsKey(interactionOpt2)) {
						confidentInteractions.put(interactionOpt1, Double.parseDouble(col[2])); // col[2] = correlation score
					}
				}
			line = input.readLine();
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return confidentInteractions;
	}
	/**
	 * Generate an index for the FASTA file; identifying the line number for the different refseq identifiers
	 * 
	 * @param fastaFile String - file path for the FASTA sequences
	 * @return indexOfFastaFile HashMap<String, Integer> - map of {refseqId : line count} 
	 */
	private static HashSet<String> generateRefSeqSet(String fastaFile){
		HashSet<String> refSeqSet = new HashSet<>();

		InputStream in;
		try {
			in = new FileInputStream(new File(fastaFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine();

			while(line != null) {
				/* store the line index of the start of a new sequence */
				if(line.startsWith(">")) {
					String[] col = line.split("_|\\.");
					String refSeqId = col[2] + "_" + col[3]; // col[2] = type of ID (eg. NM) ; col[3] = number ID (#####)

					refSeqSet.add(refSeqId);
				}

				line = input.readLine();
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return refSeqSet;
	}
	
	/**
	 * Get the list of refSeqIds associated to a given protein. RefSeqIds need to have a corresponding sequence in the fasta file
	 * 
	 * @param mapProtToRefSeqIdsFile	String - file path output from BiomaRt outlining map of HGNC symbol to RefSeqId
	 * @param refSeqSet					HashSet<String> - refSeqIds with a corresponding sequence in the fasta file
	 * @return mapRefSeqIds				HashMap<String, ArrayList<String>> - map of protein as HGNC symbol = list of RefSeqIds
	 */
	private static HashMap<String, ArrayList<String>> getRefSeqIdsInNetwork(String mapProtToRefSeqIdsFile, HashSet<String> refSeqSet){
		
		HashMap<String, ArrayList<String>> mapRefSeqIds = new HashMap<>();
		
		try {
			InputStream in = new FileInputStream(new File(mapProtToRefSeqIdsFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));
			
			String line = input.readLine(); // header
			line = input.readLine();
			while(line!=null) {
				
				String protein = line.split("\t")[0];
				String refSeqID = line.split("\t")[1];
				
				/* Add RefSeqId only if it has corresponding sequence in fasta file */
				if(refSeqSet.contains(refSeqID)) {
					
					/* Add refSeqId to existing list if protein has been seen b4, 
					 * or create new listing if protein hasn't been seen b4 */ 
					if(mapRefSeqIds.containsKey(protein)) {
						mapRefSeqIds.get(protein).add(refSeqID);
					} else {
						ArrayList<String> refSeqIdList = new ArrayList<>();
						refSeqIdList.add(refSeqID);
						mapRefSeqIds.put(protein, refSeqIdList);
					}
				}
				line = input.readLine();
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return mapRefSeqIds;
	}
	
	private static HashSet<String> getOverlyConnectedProteins(HashMap<String, Double> confidentInteractionsMap, int maxInteractions){
		
		HashSet<String> proteinsToRemove = new HashSet<>();
		
		/* Count number of interactions each protein is involved in */
		HashMap<String, Integer> interactionCountByProtein = new HashMap<>();
		for(String interaction : confidentInteractionsMap.keySet()) {
			String interactor1 = interaction.split("\\_")[0];
			String interactor2 = interaction.split("\\_")[1];
			
			/* update count for first interactor */
			if(interactionCountByProtein.containsKey(interactor1)){
				interactionCountByProtein.put(interactor1, interactionCountByProtein.get(interactor1) + 1);
			} else {
				interactionCountByProtein.put(interactor1, 1);
			}
			
			/* update count for second interactor */
			if(interactionCountByProtein.containsKey(interactor2)){
				interactionCountByProtein.put(interactor2, interactionCountByProtein.get(interactor2) + 1);
			} else {
				interactionCountByProtein.put(interactor2, 1);
			}
		}

		/* Check if protein is involved in more than the accepted number of interactions */
		for(Entry<String, Integer> entry : interactionCountByProtein.entrySet()) {
			if(entry.getValue() >= maxInteractions) {
				proteinsToRemove.add(entry.getKey());
			}
		}
		
		return proteinsToRemove;
	}
	
	private static ArrayList<Interaction> formatInteractionList(HashMap<String, Double> confidentInteractions, HashMap<String, ArrayList<String>> mapProtToRefSeqIds ){

		ArrayList<Interaction> interactionList = new ArrayList<>(); // list to contain formated interactions

		for(String interaction: confidentInteractions.keySet()) {

			/* bait and prey identifiers for significant saint PPI */
			String prot1 = interaction.split("\\_")[0];
			String prot2 = interaction.split("\\_")[1];


			/* Interaction constructor: 
			 * Interaction(String _Protein1, String _Protein2, List<String> _ID1, List<String> _ID2, double _w)
			 * Interaction(Bait-gene-symbol, Prey-gene-symbol, Bait-RefSeq-mRNA-Acc, Prey-RefSeq-mRNA-Acc, weight)*/
			Interaction ppi = new Interaction(prot1, prot2, mapProtToRefSeqIds.get(prot1), mapProtToRefSeqIds.get(prot2), (1-confidentInteractions.get(interaction)));
			
			interactionList.add(ppi);
		}

		return interactionList;
	}
	
	private static ArrayList<Interaction> formatInteractionListWithoutOverelyConnectedProteins(HashMap<String, Double> confidentInteractions, 
			HashMap<String, ArrayList<String>> mapProtToRefSeqIds, HashSet<String> proteinsToRemove){

		ArrayList<Interaction> interactionList = new ArrayList<>(); // list to contain formated interactions

		for(String interaction: confidentInteractions.keySet()) {

			/* bait and prey identifiers for significant saint PPI */
			String prot1 = interaction.split("\\_")[0];
			String prot2 = interaction.split("\\_")[1];

			/* add interaction if both protein interactors aren't involved in more than the max number of acceptable interactions */
			if(!proteinsToRemove.contains(prot1) && !proteinsToRemove.contains(prot2)) {
				/* Interaction constructor: 
				 * Interaction(String _Protein1, String _Protein2, List<String> _ID1, List<String> _ID2, double _w)
				 * Interaction(Bait-gene-symbol, Prey-gene-symbol, Bait-RefSeq-mRNA-Acc, Prey-RefSeq-mRNA-Acc, weight)*/
				Interaction ppi = new Interaction(prot1, prot2, mapProtToRefSeqIds.get(prot1), mapProtToRefSeqIds.get(prot2), (1-confidentInteractions.get(interaction)));
				interactionList.add(ppi);
			}
			
			
			
		}

		return interactionList;
	}

	
	/**
	 * Generate the list of proteins in the network from the list of confident interactions >> for testing purposes
	 * @param confidentInteractions
	 * @return
	 */
	@SuppressWarnings("unused")
	private static HashSet<String> getConfidentProteins(HashMap<String, Double> confidentInteractions){
		HashSet<String> confidentProteinSet = new HashSet<>();
		
		for(String interaction: confidentInteractions.keySet()) {
			
			String prot1 = interaction.split("\\_")[0];
			String prot2 = interaction.split("\\_")[1];
			
			confidentProteinSet.add(prot2);
			confidentProteinSet.add(prot1);
		}
		
		return confidentProteinSet;
	}
	
	@SuppressWarnings("unused")
	private static void printProteinsInNetwork(HashSet<String> confidentProteinSet, String outputFile) {
		
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			
			for(String protein: confidentProteinSet) {
				out.write(protein + "\n");
				out.flush();
			}
			
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
}
