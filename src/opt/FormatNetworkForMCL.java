package opt;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

import graph.Interaction;
import graph.Protein;

public class FormatNetworkForMCL {
	
	public static void formatMCLnetwork(ArrayList<Protein> proteinList, ArrayList<Interaction> interactionList, String networkOutputFile) {
		
		/* transform protein list to set<ProteinName> */ 
		HashSet<String> proteinSet = getProteinSet(proteinList);
		int count=0;
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(networkOutputFile)));
			
			for(Interaction inter: interactionList) {
			
				String prot1 = inter.getProtein1();
				String prot2 = inter.getProtein2();
				
				/* print interaction if both interactors (proteins) are found in the final network */
				if(proteinSet.contains(prot1) && proteinSet.contains(prot2)) {
					out.write(prot1 + " " + prot2 + " " + inter.getWeight() + "\n");
					out.flush();
					
					count++;
				}
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		System.out.println("Number of interactions in final network: " + count);
		System.out.println("Number of proteins in final network: " + proteinSet.size());
	}
	
	private static HashSet<String> getProteinSet(ArrayList<Protein> proteinList){
		HashSet<String> proteinSet = new HashSet<>();
		
		for(Protein prot: proteinList) {
			proteinSet.add(prot.getProteinName());
		}
		return proteinSet;
	}
	
}
