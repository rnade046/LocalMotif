package localizations;

public class MotifLocalizations {

	private String motif;
	private int motifCount;
	private int[] safeDomains; 
	
	public MotifLocalizations(String m, int count, int[] domains) {
		motif = m;
		motifCount = count;
		safeDomains = domains;
	}
	
	public String getMotif() {
		return this.motif;
	}
	
	public int getCount() {
		return this.motifCount;
	}
	
	public int[] getDomainCount() {
		return this.safeDomains;
	}
}
