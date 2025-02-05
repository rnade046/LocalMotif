package summarizeMotifs;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import ClusteredMotifs.Family;

public class ClusteredMotif {

	private String motif;
	private double clustering;
	private double pval;
	private String fdr;
	private int numProts;
	private String[] protList;
	private boolean family;
	private int familyNumber;

	public ClusteredMotif(String[] info, ArrayList<Family> motifFamilies, List<Double[]> fdrInfo, HashMap<String, String> proteinMap) {
		motif = info[0];
		numProts = (int) Math.ceil(0.4*Integer.parseInt(info[1]));
		clustering = Double.parseDouble(info[2]);
		pval = Double.parseDouble(info[3]);

		setFamily(motifFamilies);
		setFDR(fdrInfo);

		setProteinList(proteinMap);
	}

	public String getMotif() {
		return this.motif;
	}

	public double getClusteringMeasure() {
		return this.clustering;
	}

	public double getPval() {
		return this.pval;
	}

	public String getFDR() {
		return this.fdr;
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

	public Integer getFamily() {
		return this.familyNumber;
	}

	private void setFamily(ArrayList<Family> motifFamilies) {

		for(Family f: motifFamilies) {
			
			if(f.getAllMotifs().contains(this.motif)) {
			
				this.familyNumber = f.getNumber();
				if(f.getRepresentativeMotif().equals(this.motif)) {
					this.family = true;
				} else {
					this.family = false;			
				}
				break;
			}
		}
	}

	private void setFDR(List<Double[]> fdrInfo) {

		for(int i=0; i<fdrInfo.size(); i++) {

			if(this.pval < fdrInfo.get(i)[0]) {

				double currentFDR = fdrInfo.get(i)[1];

				if(currentFDR == 0) {
					/* search for none zero FDR */ 
					double nonZeroFDR = getNonZeroFDR(fdrInfo);
					this.fdr = "<" + Double.toString(nonZeroFDR);

				} else {
					this.fdr = Double.toString(currentFDR);
				}
				break;
			}
		}
	}
	
	public boolean isRepresentative() {
		return this.family;
	}

	private double getNonZeroFDR(List<Double[]> fdrInfo) {

		double fdr = 0;
		for (int i=0; i<fdrInfo.size(); i++) {

			if(fdrInfo.get(i)[1] > 0) {
				fdr = fdrInfo.get(i)[1];
				break;
			}
		}
		return fdr;
	}
	
	private void setProteinList(HashMap<String, String> proteinMap) {
		
		if(proteinMap.containsKey(this.motif)) {
			this.protList = proteinMap.get(this.motif).split("\\|");
		}
	}
}
