package localizations;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

public class FormatHCMlocalizationForEnrichment {

	public static void main(String[] args) {
		
		String wd = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/";
		
		String proteinsInNetworkFile = wd + "corrNetTop2-400_proteinsInNetwork_info.tsv";
		String preyFile = wd + "localization/preys-latest.txt"; 
		String motifsInterestFile = wd + "corrNetTop2-400_coreTPD_p0.4_coreProteins_h0.7_motifFamiliesInfo.tsv"; 
		String annotationFile = wd + "corrNetTop2-400_coreTPD_p0.4_coreProteinsByMotif.tsv";
		String outputFile = wd + "localization/SAFE/safeEnrichmentInput_motif";
		String outputNMFFile = wd + "localization/NMF/nmfEnrichmentInput_motif";
		
		
		formatLocalizations(proteinsInNetworkFile, preyFile, motifsInterestFile, annotationFile, outputFile, outputNMFFile);
		
	}


	public static void formatLocalizations(String proteinsInNetworkFile, String preyFile, String motifsInterestFile, String annotationFile, String outputFile, String outputNMFFile) {

		/* get proteins in network SAFE domain */
		System.out.println("** load prey info **");
		List<Protein> proteinList = loadPreyInfo(preyFile, loadProteinsInNetwork(proteinsInNetworkFile));
		
		/* determine output of each motif */
		
		/* load motifs */
		System.out.println("** load motifs to test **");
		HashMap<String, Integer> motifs = getMotifsToTest(motifsInterestFile);
		
		/* for each motif get list of proteins */ 
		System.out.println("** get motif info **");
		proteinList = determineMotifsAssociatedToProteins(annotationFile, motifs, proteinList);
		
		/* print protein summary per motif */
		System.out.println("** print motifs **");
		printSafeDomains(outputFile, motifs, proteinList);
		printNMFregions(outputNMFFile, motifs, proteinList);
	}

	private static HashSet<String> loadProteinsInNetwork(String proteinsInNetworkFile){

		HashSet<String> proteinSet = new HashSet<>();
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(proteinsInNetworkFile))));
			
			String line = in.readLine();
			while(line !=null) {
				
				proteinSet.add(line.split("\t")[0]);
				line = in.readLine();
			}
			
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return proteinSet;
	}
	
	private static List<Protein> loadPreyInfo(String preyFile, HashSet<String> proteinSet){ 
		
		List<Protein> proteinList = new ArrayList<>();
		
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(preyFile))));
			
			String line = in.readLine(); // header 
			line = in.readLine();
			
			while(line!=null) {
				
				String[] col = line.split("\t");
				
				if(proteinSet.contains(col[0])) {
					
					int safe = 25;
					int nmf = 19;
					
					if(!col[3].equals("-")) {
						safe = Integer.parseInt(col[3]);
					}
					
					if(!col[1].equals("19")) {
						nmf = Integer.parseInt(col[1]);
					}
					
					proteinList.add(new Protein(col[0], safe, nmf));
				}
				
				line = in.readLine();
			}
			
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	
		return proteinList;
	}

	private static HashMap<String, Integer> getMotifsToTest(String motifsInterestFile){
		
		HashMap<String, Integer> motifSet = new HashMap<>();

		try {
			InputStream in = new FileInputStream(new File(motifsInterestFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // header
			line = input.readLine();
			int count = 1;

			while(line != null) {

				motifSet.put(line.split("\t")[0], count);
				line = input.readLine();
				count++;
			}

			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return motifSet;
	}
	
	private static List<Protein> determineMotifsAssociatedToProteins(String annotationFile, HashMap<String, Integer> motifFamilies, List<Protein> proteinList){
		
		/* index protein list */
		HashMap<String, Integer> proteinIdxMap = new HashMap<>();
		for(int i=0; i<proteinList.size(); i++) {
			proteinIdxMap.put(proteinList.get(i).getName(), i);
		}
		
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(annotationFile))));
			
			String line = in.readLine();
			while(line!=null) {
				
				String motif = line.split("\t")[0];
				if(motifFamilies.containsKey(motif)) {
					
					int motifIdx = motifFamilies.get(motif);
;					String[] proteins = line.split("\t")[2].split("\\|");
					
					for(String p : proteins) {
						if(proteinIdxMap.containsKey(p)) {
							proteinList.get(proteinIdxMap.get(p)).addMotif(motifIdx);
						}
					}
					
				}
				
				line = in.readLine();
			}				
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return proteinList;
	}
	
	private static void printSafeDomains(String outputFile, HashMap<String, Integer> motifFamilies, List<Protein> proteinList) {
		
		for(int motif: motifFamilies.values()) {
			
			try {
				BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile + motif + ".tsv")));

				out.write("protein\tgroup\tSAFE\n");
				
				for(Protein prot : proteinList) {
					out.write(prot.getName() + "\t");
					
					if(prot.containsMotif(motif)) {
						out.write("Motif\t");
					} else {
						out.write("NotMotif\t");
					}
					
					out.write(prot.getSafeDomain() + "\n");
					out.flush();
				}
				
				out.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
			
		}
		
	}
	
private static void printNMFregions(String outputFile, HashMap<String, Integer> motifFamilies, List<Protein> proteinList) {
		
		for(int motif: motifFamilies.values()) {
			
			try {
				BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile + motif + ".tsv")));

				out.write("protein\tgroup\tSAFE\n");
				
				for(Protein prot : proteinList) {
					out.write(prot.getName() + "\t");
					
					if(prot.containsMotif(motif)) {
						out.write("Motif\t");
					} else {
						out.write("NotMotif\t");
					}
					
					out.write(prot.getNMFregion() + "\n");
					out.flush();
				}
				
				out.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
			
		}
		
	}
}
