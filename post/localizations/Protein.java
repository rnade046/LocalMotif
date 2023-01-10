package localizations;

import java.util.HashSet;

public class Protein {
	
	private String proteinName;
	private int safeDomain;
	private int nmfRegions;
	private HashSet<Integer> listOfMotifs;
	
	
	public Protein(String name, int domain, int region) {
		this.proteinName = name;
		this.safeDomain = domain;
		this.nmfRegions = region;
		this.listOfMotifs = new HashSet<>();
	}
	
	public String getName() {
		return this.proteinName;
	}
	
	public int getSafeDomain() {
		return this.safeDomain;
	}
	
	public int getNMFregion() {
		return this.nmfRegions;
	}
	
	public void addMotif(int motif) {
		this.listOfMotifs.add(motif);
	}
	
	public boolean containsMotif(int motif) {
		return this.listOfMotifs.contains(motif);
	}
}
