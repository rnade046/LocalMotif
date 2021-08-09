package opt;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;

import graph.Interaction;
import graph.Protein;

public class CheckDegreeDistributions {

	public static void assessDegreeDistribution(ArrayList<Protein> proteinsInNetwork, ArrayList<Interaction> interactionList, String outputFile) {
		
		/* Make Map of protein names in network */
		HashMap<String, Integer> mapOfProteinDegrees = new HashMap<>();
		for(Protein prot: proteinsInNetwork) {
			mapOfProteinDegrees.put(prot.getProteinName(), 0);
		}
		
		/* Count # of interactions each protein is involved in */ 
		for(Interaction inter: interactionList) {
			
			/* insure both interactors are in the final network */
			if(mapOfProteinDegrees.containsKey(inter.getProtein1()) && mapOfProteinDegrees.containsKey(inter.getProtein2())) {
				mapOfProteinDegrees.put(inter.getProtein1(), mapOfProteinDegrees.get(inter.getProtein1()) + 1);
				mapOfProteinDegrees.put(inter.getProtein2(), mapOfProteinDegrees.get(inter.getProtein2()) + 1);
			}
		}
			
		/* print map */
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			
			out.write("Protein\t#Degrees\n");
			for(Entry<String, Integer> entry: mapOfProteinDegrees.entrySet()) {
				out.write(entry.getKey() + "\t" + entry.getValue() + "\n");
				out.flush();
			}
			
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	
}
