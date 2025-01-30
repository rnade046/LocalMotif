package ClusteredMotifs;

import java.util.HashSet;

public class Family {
	
	private String representativeMotif;
	private int number;
	private HashSet<String> allMotifs;
	
	public Family(String m, int n, HashSet<String> motifSet) {
		this.representativeMotif = m;
		this.number = n;
		this.allMotifs = motifSet;
	}
	
    // Getter for all attributes
    public String getFamilyDetails() {
    	
    	String motifs = "";
    	for(String m: allMotifs) {
    		motifs+=m + "|";
    	}
        return String.format("%s\t%d\t%s", representativeMotif, number, motifs);
    }

    // Individual getters (optional)
    public String getRepresentativeMotif() {
        return representativeMotif;
    }

    public int getNumber() {
        return number;
    }

    public HashSet<String> getAllMotifs() {
        return allMotifs;
    }

}
