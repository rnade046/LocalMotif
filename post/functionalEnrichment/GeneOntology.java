package functionalEnrichment;

public class GeneOntology {

	private String name;
	private String description;
	private double pval;
	private double enrichmentScore;
	
//	public GeneOntology(String goName, String goDescription, double goPval, int studyTerm, int studyTotal) {
//		
//		this.name = goName;
//		this.description = goDescription; 
//		this.pval = goPval;
//		this.enrichmentScore = studyTerm / (double) studyTotal;
//	}
	
	public GeneOntology(String[] col) {
		
		this.name = col[0];
		this.description = col[12]; 
		this.pval = Double.parseDouble(col[10]);
		
		
		double studyTerm = Double.parseDouble(col[4]);
		double studyTotal = Double.parseDouble(col[3]);
		double popTerm = Double.parseDouble(col[2]);
		double popTotal = Double.parseDouble(col[1]);
		
		this.enrichmentScore = (studyTerm/studyTotal) / (popTerm/popTotal);
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
