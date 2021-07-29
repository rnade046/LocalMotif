package utils;

import java.security.acl.LastOwnerException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;

public class TopPercentPairwiseDistance {

	
    /**
     * Computes the total pairwise distance from a list of protein indexes of
     * interest corresponding to proteins in the network, using the distances found
     * in the distance matrix
     * 
     * @param distance_matrix
     * @param proteinIdxList : list of protein indexes in distance matrix
     * @return total pairwise distance
     **/
    public static double computeTPD(ArrayList<Integer> proteinIdxList, double[][] distance_matrix) {
        

        double distance = 0; // initialize distance

        // read distance matrix to get sum of distances between significant proteins
        for (int i = 0; i < proteinIdxList.size(); i++) {
            for (int j = i + 1; j < proteinIdxList.size(); j++) {
            	
                if (i > distance_matrix.length || j > distance_matrix.length) {
                    System.out.println("error");
                }
                
                if (distance_matrix[proteinIdxList.get(i)][proteinIdxList.get(j)] == Double.MAX_VALUE || distance == Double.MAX_VALUE) {
                    distance = Double.MAX_VALUE;
                } else {
                    // indexes of proteins of interests are found in the array list indexProt
                    distance += distance_matrix[proteinIdxList.get(i)][proteinIdxList.get(j)];
                }
                
            }
        }
        if(distance == 0.0) {
        	System.out.println("tpd = 0");
        }
        
        return (double) Math.round(distance * 100d) / 100d;

    }
	
	/**
	 * For a given list of protein indexes within the network, identify the core proteins, and compute their top % pairwise distance. 
	 * 
	 * @param annotatedProteinsIndexes	List<Integer> index of proteins annotated in the network
	 * @param distanceMatrix			double[][] distance matrix of shortest weighted path of nodes in the network
	 * @param percentThreshold			double percent threshold cut off set by user
	 * 
	 * @return tppd						double top percent pairwise distance
	 */
	public static double getTPPD(ArrayList<Integer> annotatedProteinsIndexes, double[][] distanceMatrix, double percentThreshold) {
		
		System.out.println("**GET TPPD**");
		/* Identify core nodes */
		int nodeThreshold = (int) Math.ceil(annotatedProteinsIndexes.size()*percentThreshold); // determine number of core proteins
		long startTime = System.currentTimeMillis();
		ArrayList<Integer> coreProteins = getCoreProteins(annotatedProteinsIndexes, distanceMatrix, nodeThreshold); 	// Identify core proteins
		long endTime = System.currentTimeMillis();
		System.out.println("*CoreProteins total: " + (endTime - startTime) + "\n");
		
		/* List sorted paths involving  */
		startTime = System.currentTimeMillis();
		ArrayList<Double> sortedPaths = getPaths(coreProteins, annotatedProteinsIndexes, distanceMatrix);
		endTime = System.currentTimeMillis();
		System.out.println("Getting + sorting shortest paths total : " + (endTime - startTime) + "\n");
		
		/* Compute top percent pairwise distance */
		startTime = System.currentTimeMillis();
		double tppd = computeTPPD(sortedPaths, percentThreshold);
		endTime = System.currentTimeMillis();
		System.out.println("TPPD total: " + (endTime - startTime) + "\n");
		
		return (double) Math.round(tppd * 1000d) / 1000d;
	}

	/**
	 * For a given list of protein indexes within the network, identify the core proteins, and compute their total pairwise distances. 
	 * 
	 * @param annotatedProteinsIndexes	List<Integer> index of proteins annotated in the network
	 * @param distanceMatrix			double[][] distance matrix of shortest weighted path of nodes in the network
	 * @param percentThreshold			double percent threshold cut off set by user
	 * 
	 * @return coreTPD					double core total pairwise distance
	 */
	public static double getCoreTPD(ArrayList<Integer> annotatedProteinsIndexes, double[][] distanceMatrix, double percentThreshold) {
		
		/* Identify core nodes */
		int nodeThreshold = (int) Math.ceil(annotatedProteinsIndexes.size()*percentThreshold); // determine number of core proteins
		ArrayList<Integer> coreProteins = getCoreProteins(annotatedProteinsIndexes, distanceMatrix, nodeThreshold); 	// Identify core proteins
		
		/* compute total pairwise distance for core proteins */
		double coreTPD = computeCoreTPD(coreProteins, distanceMatrix);
		
		return (double) Math.round(coreTPD * 1000d) / 1000d;
	}
	
	/**
	 * Identify core proteins from a list of proteins in the network. Core proteins are the top % nodes with the lowest sum of pairwise distance 
	 * with all other annotated proteins. The list of core proteins indexes in the network is returned.
	 * 
	 * @param proteinIndexes	List<Integer> index of annotated proteins in the network
	 * @param distanceMatrix	double[][]	distance matrix of shortest paths
	 * @param nodeThreshold		int threshold for protein core
	 * 
	 * @return coreProteinsList	List<Integer> indexes of proteins in the core (most clustered annotated proteins)
	 */
	private static ArrayList<Integer> getCoreProteins(ArrayList<Integer> proteinIndexes, double[][] distanceMatrix, int nodeThreshold){
		ArrayList<Integer> coreProteinsList = new ArrayList<>(); 

		/* Calculate the sum of interaction weights for the protein to all proteins annotated */
		long startTime = System.currentTimeMillis();
		HashMap<Integer, Double> mapOfSumInteractionWeights = new HashMap<>();
		for(int proteinIdx : proteinIndexes) {
			double sumOfWeights = 0;

			for(int j: proteinIndexes) {
				sumOfWeights += distanceMatrix[proteinIdx][j];			
			}

			mapOfSumInteractionWeights.put(proteinIdx, sumOfWeights);
		}
		long endTime = System.currentTimeMillis();
		System.out.println("coreProts sum: " + (endTime - startTime));
		
		/* Sort map by value */
		
		startTime = System.currentTimeMillis();
		List<Entry<Integer, Double>> rankedProteinslist = new ArrayList<>(mapOfSumInteractionWeights.entrySet());
		rankedProteinslist.sort(Entry.comparingByValue());
		endTime = System.currentTimeMillis();
		
		System.out.println("rank sums: " + (endTime - startTime));
		
		startTime = System.currentTimeMillis();
		/* Identify top % nodes */ 
		for(int i=0; i<nodeThreshold; i++) {
			coreProteinsList.add(rankedProteinslist.get(i).getKey());
		}
		endTime = System.currentTimeMillis();
		System.out.println("get core nodes: " + (endTime -startTime));
		
		return coreProteinsList;
	}

	/**
	 * Given a list of core proteins compute their total pairwise distance (sum of shortest paths between all proteins).
	 * 
	 * @param coreNodesIdxs		List<Integer> indexes of core annotated proteins
	 * @param distanceMatrix	double[][] distance matrix
	 * 
	 * @return coreTPD	double the total pairwise distance computed for the set of core proteins
	 */
	private static double computeCoreTPD(ArrayList<Integer> coreNodesIdxs, double[][] distanceMatrix) {

		/* Get all paths for core nodes */
		double coreTPD=0;

		for(int i=0; i<coreNodesIdxs.size(); i++) {
			for(int j=i+1; j<coreNodesIdxs.size(); j++) {
				coreTPD += distanceMatrix[coreNodesIdxs.get(i)][coreNodesIdxs.get(j)];
			}
		}
		
		return coreTPD;
	}
	
	
	/**
	 * Identify the list of paths the core proteins are part of with all annotated proteins. 
	 * The list is sorted such that smaller paths are at the top.
	 * 
	 * @param coreIdxs			List<Integer> indexes of core annotated proteins
	 * @param allIdxs			List<Integer> index of annotated proteins in the network
	 * @param distanceMatrix	double[][] distance matrix
	 * 
	 * @return allPaths		List<Double> - sorted list of paths
	 */
	private static ArrayList<Double> getPaths(ArrayList<Integer> coreIdxs, ArrayList<Integer> allIdxs, double[][] distanceMatrix){

		ArrayList<Double> allPaths = new ArrayList<>();
		
		/* sort allIdx */
		long startTime = System.currentTimeMillis();
		ArrayList<Integer> orderedIdxs = new ArrayList<>();
		orderedIdxs.addAll(coreIdxs);
		
		for(int idx: allIdxs) {
			if(!coreIdxs.contains(idx)) {
				orderedIdxs.add(idx);
			}
		}
		
		long endTime = System.currentTimeMillis();
		System.out.println("order idxs: " + (endTime - startTime));
		
		startTime = System.currentTimeMillis();
		/* get all paths */
		for(int i=0; i<coreIdxs.size(); i++) {
			for(int j=i+1; j<orderedIdxs.size(); j++) {
				allPaths.add(distanceMatrix[coreIdxs.get(i)][orderedIdxs.get(j)]);
			}
		}
		endTime = System.currentTimeMillis();
		System.out.println("get paths: " + (endTime - startTime));
		
		startTime = System.currentTimeMillis();
		Collections.sort(allPaths);   
		endTime = System.currentTimeMillis();
				
		System.out.println("sort paths: " + (endTime-startTime));
		
		return allPaths;
	}
	
	private static ArrayList<Double> getTopPaths(ArrayList<Integer> coreIdxs, ArrayList<Integer> allIdxs, double[][] distanceMatrix, double threshold){
		ArrayList<Double> sortedPaths = new ArrayList<>();
		
		/* sort allIdx */
		ArrayList<Integer> orderedIdxs = new ArrayList<>();
		HashSet<Integer> coreIdxSet = new HashSet<>(coreIdxs);
		orderedIdxs.addAll(coreIdxs);
		
		for(int idx: allIdxs) {
			if(!coreIdxSet.contains(idx)) {
				orderedIdxs.add(idx);
			}
		}
		
		/* Calculate total paths */
		int totalPaths = 0;
		for(int i=0; i<coreIdxs.size(); i++) {
			totalPaths += (orderedIdxs.size()-(i+1));
		}
		
		/* determine max number of paths that will be considered */
		int maxPaths = (int) Math.ceil(totalPaths * threshold);
		
		for(int i=0; i<coreIdxs.size(); i++) {
			for(int j=i+1; j<orderedIdxs.size(); j++) {
				
				double currentPath = distanceMatrix[coreIdxs.get(i)][orderedIdxs.get(j)];
				
				if(sortedPaths.size() < 1 ) {
					sortedPaths.add(currentPath);
				} else if(sortedPaths.size() < maxPaths) {
					/* add elements to list, ensuring it stays sorted with binary search */
					int idx = binarySearch_leftmostIdx(sortedPaths, currentPath);
					sortedPaths.add(idx, currentPath);
				} else if(currentPath <= sortedPaths.get(maxPaths-1)) {
					/* add element in a sorted fashion if the current paths value is smaller than the max path + remove max element */
					int idx = binarySearch_leftmostIdx(sortedPaths, currentPath);
					sortedPaths.add(idx, currentPath);
					sortedPaths.remove(maxPaths-1);
				}
			}
		}
		return sortedPaths;
	}
	
	private static int binarySearch_leftmostIdx(ArrayList<Double> sortedPaths, double currentPath) {
		
		int l = 0;
		int r = sortedPaths.size();

		/* protein corresponding to selected random weight will be the first protein to be greater or equal to the weight */
		while(l < r) {
			int m = (int) Math.floor((l+r)/ (double)2);
			if(sortedPaths.get(m) < currentPath) {
				l = m + 1;
			} else { 
				r = m;
			}
		}
		
		return l;
	}
	
	/**
	 * Compute the top percent pairwise distance. Measure the sum of the top % paths.
	 * 
	 * @param sortedPaths		List<Double> list of paths between core and all other proteins in annotation
	 * @param percentThreshold	double percentage threshold to cut off interaction
	 * 
	 * @return tppd		double top percent pairwise distance
	 */
	private static double computeTPPD(ArrayList<Double> sortedPaths, double percentThreshold) {
		
		double tppd = 0;
		int threshold = (int) Math.ceil(sortedPaths.size()*percentThreshold);
		
		for(int i=0; i<threshold; i++) {
			tppd += sortedPaths.get(i);
		}
		
		return tppd;
	}
	
}
