package formatHCM;
public class Interaction {

	// fields
	private String Protein1;
	private String Protein2;
	private double correlation;

	// constructors
	public Interaction(String _Protein1, String _Protein2, double _w) {
		Protein1 = _Protein1;
		Protein2 = _Protein2;
		correlation = _w;
	}

	// get
	public String getProtein1() {
		return Protein1;
	}

	public String getProtein2() {
		return Protein2;
	}

	public double getWeight() {
		return this.correlation;
	}

	public void setWeight(double w) {
		this.correlation = w;
	}
}
