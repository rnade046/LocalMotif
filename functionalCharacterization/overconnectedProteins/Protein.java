package overconnectedProteins;

public class Protein {

	private String name;
	private String nmf;
	private String safe; 
	
	
	public Protein(String id, int s, int n) {
		this.name = id;
		this.nmf = String.valueOf(n);
		this.safe = String.valueOf(s);
	}
	
	public void setNMF(String value) {
		this.nmf = value;
	}
	
	public void setSafe(String value) {
		this.safe = value;
	}
	
	public String[] getProteinInfo() {
		return new String[] {this.name, this.nmf, this.safe};
	}
	
}
