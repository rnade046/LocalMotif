package graph;

import java.util.ArrayList;

public class Interaction {
		
	private final static double DEFAULT_WEIGHT = 1.0;
	
		// fields
		private String Protein1;
		private String Protein2;
		private ArrayList<String> ID1;
		private ArrayList<String> ID2;
		private double weight;

		// constructors
		public Interaction(String _Protein1, String _Protein2, ArrayList<String> _ID1, ArrayList<String> _ID2) {
			this(_Protein1, _Protein2, _ID1, _ID2, DEFAULT_WEIGHT);
		}
		
		public Interaction(String _Protein1, String _Protein2, ArrayList<String> _ID1, ArrayList<String> _ID2, double _w) {
			Protein1 = _Protein1;
			Protein2 = _Protein2;
			ID1 = _ID1;
			ID2 = _ID2;
			weight = _w;
		}

		// get
		public String getProtein1() {
			return Protein1;
		}

		public String getProtein2() {
			return Protein2;
		}

		public ArrayList<String> getID1() {
			return ID1;
		}

		public ArrayList<String> getID2() {
			return ID2;
		}

		public double getWeight() {
			return this.weight;
		}
		
		public void setWeight(double w) {
			this.weight = w;
		}
}

