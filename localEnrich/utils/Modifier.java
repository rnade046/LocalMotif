package utils;

import java.util.ArrayList;
import java.util.HashMap;

import graph.Annotation;
import graph.Interaction;
import graph.Protein;

public class Modifier {


	

	
	public static void modifyInteractionWeight(ArrayList<Protein> networkProteins, ArrayList<Interaction> networkInteractions) {

		int countNotWeightedInteractions = 0;
		int countWeightedInteractions = 0;

		for(int i=0; i<networkInteractions.size(); i++){

			Interaction inter = networkInteractions.get(i);
			String iprot1 = inter.getProtein1();
			String iprot2 = inter.getProtein2();


			/* Find the fold change associated to each protein in the interaction */
			boolean prot1found = false;
			boolean prot2found = false;
			int protCount = 0;

			double foldchange1 = 1.0;
			double foldchange2 = 1.0; 

			/* Search through the list of network proteins, for the two proteins in this interaction 
			 * proteinCount < size of protein list exist since disconnect proteins have been removed from the
			 * network  */
			while ((prot1found == false || prot2found == false) && protCount < networkProteins.size()) {

				Protein prot1 = networkProteins.get(protCount);

				if ((prot1.getProteinName().equals(iprot1) || prot1.getProteinName().equals(iprot2)) && (prot1found == false && prot2found == false)) {
					prot1found = true;
					foldchange1 = prot1.getFoldChange();
				} else if ((prot1.getProteinName().equals(iprot1) || prot1.getProteinName().equals(iprot2)) && (prot1found == true && prot2found == false) ) {
					prot2found = true;
					foldchange2 = prot1.getFoldChange();
				} 
				protCount++;
			}

			/* Compute strength of interaction based on fold change. 
			 * If the average fold change between proteins is > 1, the strength becomes 1/average_fc
			 * If the average fold change between proteins is =< 1, the strength remains the average_fc */

			double average_fc = (foldchange1 + foldchange2)/2; // compute average fold change

			if (average_fc > 1) {
				inter.setWeight(1/average_fc);
				countWeightedInteractions++;
			} else if (average_fc < 1) {
				inter.setWeight(average_fc);
				countWeightedInteractions++;
			}

			if(average_fc == 1) {
				countNotWeightedInteractions++;
			}

		}

	System.out.println("Weighted interactions " + countWeightedInteractions + "; total interactions " + networkInteractions.size() +
				"; percentage of network weighted "+ (networkInteractions.size() - countNotWeightedInteractions)/(double)networkInteractions.size() + "\n");

	} // end modifyInteractionWeight
	



	public static void setClusterTPD(ArrayList<Annotation> clusterList, double[][] distance_matrix) {
		/**********************************************************************************
		 * Computes the total pairwise distance from a list of protein indexes of
		 * interest corresponding to proteins in the network, using the distances found
		 * in the distance matrix
		 ***********************************************************************************/
		for (int h = 0; h < clusterList.size(); h++) {
			Annotation cluster = clusterList.get(h);
			double tpd = TopPercentPairwiseDistance.computeTPD(cluster.getIdxProteinsList(), distance_matrix);
			cluster.setTPD(tpd);
			if(h%100 == 0) {
				System.out.print(h + "|");
			}
			if(h%1000 == 0 && h!=0) {
				System.out.println();
			}
		}
		System.out.println();
	}

	/**
	 * Sets TTPD for the cluster array.
	 *
	 * @param clusters clusters array
	 * @param distance_matrix distance matrix to calculate TTPD
	 * @param l parameter to calculate a criterion for top and core proteins selection
	 */
	public static void setClustersTTPD(ArrayList<Annotation> clusters, double[][] distance_matrix, double l){
		for (Annotation cluster: clusters){
			setClusterTTPD(cluster, distance_matrix, l);
		}
	}

	/**
	 * Sets a TTPD of the cluster.
	 *
	 * @param cluster target cluster
	 * @param distance_matrix distance matrix to calculate ttpd
	 * @param l parameter to calculate a criterion for top and core proteins selection
	 */
	public static void setClusterTTPD(Annotation cluster, double[][] distance_matrix, double l){
		// maps proteins to distance sums between the protein and its neighbourhood (also called D^l_m(u) in LESMON paper
		HashMap<Integer, Double> proteinDistancesSumMapper = new HashMap<>();
		for (Integer protein: cluster.getIdxProteinsList()){
			// sum of all distances between current protein and all other proteins in the GoTerm
			double distancesSum = 0;
			for (Integer neighbourProtein: cluster.getIdxProteinsList()){
				if (neighbourProtein != protein){
					distancesSum += distance_matrix[neighbourProtein][protein];
				}
			}
			// distance criterion
			double distance = distancesSum * l;
			// sum of all distances between the protein and its neighbourhood (also called N^l_m(u) in LESMON paper)
			double closestDistancesSum = 0;
			for (Integer neighbourProtein: cluster.getIdxProteinsList()){
				if (!neighbourProtein.equals(protein) && distance_matrix[neighbourProtein][protein] <= distance){
					closestDistancesSum += distance_matrix[neighbourProtein][protein];
				}
			}
			proteinDistancesSumMapper.put(protein, closestDistancesSum);
		}
		double distancesSum = 0;
		for (Integer key: proteinDistancesSumMapper.keySet()){
			distancesSum += proteinDistancesSumMapper.get(key);
		}
		double coreDistancesSum = distancesSum * l;
		double ttpd = 0;
		for (Integer key: proteinDistancesSumMapper.keySet()){
			if (proteinDistancesSumMapper.get(key) <= coreDistancesSum){
				ttpd += proteinDistancesSumMapper.get(key);
			}
		}
		cluster.setTTPD(ttpd);
	}

}
