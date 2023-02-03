package metrics;


import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class Motif {

	private String motif;
	private int motifNumber;
	private int numProts;
	private String[] protList;
	private List<Integer> proteinIdxs;
	private double diameter;

	public Motif(String m, int number, String[] prots, HashMap<String, Integer> protIdxMap, double[][] dm) {
		motif = m;
		motifNumber = number;
		protList = prots;
		numProts= prots.length;
		
		proteinIdxs = setProteinIndexes(protIdxMap);
		diameter = determineDiameter(dm);
	}

	public String getMotif() {
		return this.motif;
	}

	public int getNumberOfProteins() {
		return this.numProts;
	}

	public String getProteinList() {
		String proteins = "";

		for(String prot: protList) {
			proteins += prot + "|";
		}
		return proteins;
	}
	
	public int getFamilyNumber() {
		return this.motifNumber;
	}
	
	public double getDiameter() {
		return this.diameter;
	}

	private List<Integer> setProteinIndexes(HashMap<String, Integer> protIdxMap){
		
		List<Integer> proteinIndexes = new ArrayList<>();
		
		for(String prot: this.protList) {
			if(protIdxMap.containsKey(prot)) {
				proteinIndexes.add(protIdxMap.get(prot));
			}
		}
		return proteinIndexes;
	}
	
	private double determineDiameter(double[][] dm) {
		
		double d = 0;
	
		for(int i=0; i<this.proteinIdxs.size(); i++) {
			for (int j=i+1; j<this.proteinIdxs.size(); j++) {
				
				double currentDiameter = dm[proteinIdxs.get(i)][proteinIdxs.get(j)];
				if(currentDiameter > d) {
					d = currentDiameter;
				}
			}
		}
		
		return d;
	}
}
