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

public class SummarizeHCMLocalizations {

	public static void main(String[] args) {
		
		String wd = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/";
		
		String preyFile = wd + "localization/preys-latest.txt";
		
		String motifFamilies = wd + "corrNetTop2-400_coreTPD_p0.4_coreProteins_h0.7_motifFamiliesInfo.tsv";
		String annotationFile = wd + "corrNetTop2-400_coreTPD_p0.4_coreProteinsByMotif.tsv";
		
		String outputFile = wd + "localization/hcm_SAFEdomains_summary_normalized_coreProts.tsv";
		
		assessSAFEdomains(preyFile, motifFamilies, annotationFile, outputFile);
	}
	
	
	public static void assessSAFEdomains(String preyFile, String motifsInterestFile, String annotationFile, String outputFile) {
		
		/* load prey - domains, protein = SAFE domain */
		HashMap<String, Integer> mapPreyLocalizations = loadPreyInfo(preyFile);
		
		//int[] totalProtCount = countProteinsPerDomain(mapPreyLocalizations);
		
		/* load motifs to test */ 
		HashMap<String, Integer> motifFamilies = getMotifsToTest(motifsInterestFile);
		
		/* for each motif determine SAFE domains */ 
		ArrayList<MotifLocalizations> domainSummary = new ArrayList<>();
		
		int[] coreProtCount = new int[24];
		HashSet<String> coreProtsToIgnore = new HashSet<>();
		
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(annotationFile))));
			
			String line = in.readLine();
			while(line!=null || domainSummary.size() < motifFamilies.size()) {
				
				String motif = line.split("\t")[0];
				if(motifFamilies.containsKey(motif)) {
					
					int[] domainCount = new int[24];
					String[] proteins = line.split("\t")[2].split("\\|");
					
					for(String p : proteins) {
						if(mapPreyLocalizations.containsKey(p)) {
							int domain = mapPreyLocalizations.get(p);
							domainCount[domain] += 1; 	
							
							if(!coreProtsToIgnore.contains(p)) {
								coreProtCount[domain] += 1;
								coreProtsToIgnore.add(p);
							}
						}
						
					}
					
					domainSummary.add(new MotifLocalizations(motif, motifFamilies.get(motif), domainCount));
				}
				
				line = in.readLine();
			}				
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		/* print to file */
		printCountSummary(outputFile, domainSummary, coreProtCount);
	}
	
	private static HashMap<String, Integer> loadPreyInfo(String preyFile){
		
		HashMap<String, Integer> mapPreyLocalizations = new HashMap<>();
		
		try {
		BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(preyFile))));
		
			String line = in.readLine(); // header 
			line = in.readLine();
			
			while(line!=null) {
				
				String[] col = line.split("\t");
				
				if(col[3].equals("-")) {
					mapPreyLocalizations.put(col[0], 0);
				} else {
					mapPreyLocalizations.put(col[0], Integer.parseInt(col[3])-1);
				}
				
				line = in.readLine();
			}
		
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return mapPreyLocalizations;
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
	
	private static void printCountSummary(String outputFile, ArrayList<MotifLocalizations> domainSummary, int[] protCount) {
		
		/* determine order of motifs */
		int[] order = new int[domainSummary.size()];
		for(int i=0; i<domainSummary.size(); i++) {
			int motifCount = domainSummary.get(i).getCount();
			order[motifCount-1] = i; 
		}
		
		/* print summary */
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			
			//header
			out.write("Family\tMotif\t");
			for(int i=1;i<=24; i++) {
				out.write(i+"\t");
			}
			out.write("\n");
			
			// summary
			for(int i=0; i<order.length; i++) {
				
				MotifLocalizations local = domainSummary.get(order[i]);
				out.write(local.getCount() + "\t" + local.getMotif() + "\t");
				
				int[] summary = local.getDomainCount();
				for(int j=0; j<summary.length; j++) {
					if(protCount[j]!=0) {
						out.write((summary[j] / (double) protCount[j]) + "\t");	
					} else {
						out.write(summary[j] + "\t");
					}
					
				}
				
				out.write("\n");
				out.flush();
			}
			
			out.write("\t\t");
			for(int i=0; i<protCount.length; i++) {
				out.write(protCount[i] + "\t");
			}
			
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	@SuppressWarnings("unused")
	private static int[] countProteinsPerDomain(HashMap<String, Integer> mapPreyLocalizations) {
		
		int[] proteinCount = new int[24];
		for(Integer domain : mapPreyLocalizations.values()) {
			
			proteinCount[domain] += 1;
		}
		
		return proteinCount;
	}
}
