package utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import graph.Interaction;
import graph.Protein;

public class NetworkProteins {

	public static ArrayList<Protein> getProteinsInNetwork(ArrayList<Interaction> InteractionList) {
		/* Find all proteins in the network form the Interaction list and store them as protein
		 * objects in an ArrayList */
		ArrayList<Protein> proteinList = new ArrayList<Protein>();

		/* Find all proteins involved in protein network and store in a HashSet.
		 * HashSets can only contain unique values, we obtain a list of unique protein names */
		HashMap<String, ArrayList<String>> proteinNetworkSet = new HashMap<String, ArrayList<String>>();

		for (int i = 0; i < InteractionList.size(); i++) {
			Interaction inter = InteractionList.get(i);

			proteinNetworkSet.put(inter.getProtein1(), inter.getID1());
			proteinNetworkSet.put(inter.getProtein2(), inter.getID2());
		}

		/*  Store protein names and IDs from HashSet in the ArrayList */
		for(String protein: proteinNetworkSet.keySet()) {

			Protein protein1 = new Protein(protein, proteinNetworkSet.get(protein)); // call Protein Class constructor, store protein name (iterator.next())

			proteinList.add(protein1); // add protein object to protein list
		}


		/*  try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File ("C://Users//Rachel//Documents//PIGNON//proteinsInStringNetwork2.tsv")));

			for(String id: proteinNetworkSet ){
				out.write(id + "\n");
				out.flush();
			}

			out.close();		
		}catch (Exception e) {
			e.printStackTrace();
		} */

		return proteinList;
	} // close getProteinsInNetwork


	/***
	 * Modify the initial list of proteins in the network based on the number the MaxValue counts in the distance matrix. 
	 * Proteins that are disconnected from the network will be removed; =proteins are disconnected when half their distance values are MAX_Values
	 *
	 * @param networkProteinsList
	 * @param proteinsToKeep
	 * @return
	 */
	public static ArrayList<Protein> modifyNetworkProteinsList(ArrayList<Protein> networkProteinsList, boolean[] proteinsToKeep, String infoFile) {


		ArrayList<Protein> networkProteinsListUpdate = new ArrayList<Protein>();

		// Run through the rows of the matrix
		for (int i = 0; i < proteinsToKeep.length; i++) {

			// if the row has less than 1/2 it's values at max value it will remain in the
			// distance matrix
			if (proteinsToKeep[i]) {
				// create new array list that will be stored in the main array
				networkProteinsListUpdate.add(networkProteinsList.get(i));
			}
		}

		File f = new File(infoFile);
		if(!f.exists() && !f.isDirectory()) {
			System.out.println("Printing protein info file");
			printProteinInfo(networkProteinsListUpdate, infoFile);
		}
		return networkProteinsListUpdate;
	}

	private static void printProteinInfo(ArrayList<Protein> proteinList, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			for(int i=0; i<proteinList.size(); i++) {
				Protein prot = proteinList.get(i);
				System.out.println(i);
				out.write(prot.getProteinName() + "\t");

				if(prot.getProteinId() != null) {
					for(String id: prot.getProteinId()) {
						out.write(id + "|");
					}
				}
				out.write("\n");
				out.flush();
			}

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}


	}

	/**
	 * Generate easy look up for the index of protein in the existing network. 
	 * 
	 * @param proteinsInNetworkList 
	 * @return proteinsIdxInNetworkMap
	 */
	public static HashMap<String, Integer> getProteinIdexInNetwork(ArrayList<Protein> proteinsInNetworkList){

		HashMap<String, Integer> proteinsIdxInNetworkMap = new HashMap<String, Integer>();

		for(int i=0; i<proteinsInNetworkList.size(); i++) {
			/* Map protein idx (value) to protein name (key) */
			proteinsIdxInNetworkMap.put(proteinsInNetworkList.get(i).getProteinName(), i); 
		}

		return proteinsIdxInNetworkMap;
	}

}
