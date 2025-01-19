package annotateMotifs;

import java.util.HashSet;

public class Motif {

	private String motifIUPAC;
	private String motifRegex;
	private HashSet<String> proteins;
	
	public Motif(String m, String regex) {
		this.motifIUPAC = m;
		this.motifRegex = regex;
		this.proteins = new HashSet<>();
	}
	
	public String getMotifIUPAC() {
		return this.motifIUPAC;
	}
	
	public String getRegexMotif() {
		return this.motifRegex;
	}
	
	public void addProtein(String p) {
		this.proteins.add(p);
	}
	
	public HashSet<String> getProteins(){
		return this.proteins;
	}
}

