package functionalEnrichment;

public class GeneOntology {

	private String name;
	private String description;
	private double pval;
	private double enrichmentScore;
	
	public GeneOntology(String goName, String goDescription, double goPval, int studyTerm, int studyTotal) {
		
		this.name = goName;
		this.description = goDescription; 
		this.pval = goPval;
		this.enrichmentScore = studyTerm / (double) studyTotal;
	}
	
	public String getName() {
		return this.name;
	}
	
	public String getDescription() {
		return this.description;
	}
	
	public double getPval() {
		return this.pval;
	}

	public double getScore() {
		return this.enrichmentScore;
	}
}
