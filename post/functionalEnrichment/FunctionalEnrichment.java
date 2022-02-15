package functionalEnrichment;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.Properties;

public class FunctionalEnrichment {

	public static void main(String[] args) throws FileNotFoundException, IOException {

		System.out.println("**Loading parameters file** \n");
		Properties params = new Properties();
		params.load(new FileInputStream(args[0]));		

		String wd = params.getProperty("working_directory");
		String networkName = params.getProperty("network_name");

		int clusteringMeasure = Integer.parseInt(params.getProperty("clusteringMeasure", "0"));
		double percentThreshold = Double.parseDouble(params.getProperty("percentThreshold", "0.2"));

		String clusteringName = "";

		switch(clusteringMeasure) {
		case 0: clusteringName = "_TPD";
		break;
		case 1: clusteringName = "_TPPD_p" + percentThreshold;
		break;
		case 2: clusteringName = "_coreTPD_p" + percentThreshold;
		break;
		}

		File directory = new File(wd + "/Ontologizer/"); 
		if (! directory.exists()){
			System.out.println("creating directory: Ontologizer/");
			directory.mkdir();
		}

		String directory2 = "/" + networkName + clusteringName + "/";
		File dir2 = new File(wd + "/Ontologizer/" + directory2);
		if(! dir2.exists()) {
			System.out.println("creating directory: Ontologizer/" + directory2);
			dir2.mkdir();
		}

		String allProtDir = "";
		String coreProtDir = "/coreProteins/";
		if(clusteringMeasure == 1 || clusteringMeasure == 2) {
			File dir3 = new File(wd+ "/Ontologizer/" + directory2 + "/allProteins/");
			if(! dir3.exists()) {
				System.out.println("creating directory: Ontologizer/" + directory2 + "/allProteins/");
				dir3.mkdir();
			}
			allProtDir = "/allProteins/";
			
			File dir4 = new File(wd + "/Ontologizer/" + directory2 + coreProtDir);
			if(! dir4.exists()) {
				System.out.println("creating directory: Ontologizer/" + directory2 + coreProtDir);
				dir4.mkdir();
			}
		}

		if(Boolean.parseBoolean(params.getProperty("prepFunctionalEnrichment"))) {
			/* Generate file of all proteins in network for analysis - if it doesn't exist */ 
			String protAnnotationFreqFile = wd + networkName + "_protFreqAnnotation.tsv";
			String proteinsInNetworkFile = wd + "Ontologizer/" +  networkName + "_proteinsInNetwork.txt";

			File f = new File(proteinsInNetworkFile);
			if(!f.exists() && !f.isDirectory()) {
				/* Create list of proteins in network from protein frequency file */
				System.out.println("**Format proteins in network**");
				formatProteinsInNetworkFile(protAnnotationFreqFile, proteinsInNetworkFile);
			}

			/* Load motif family info file to guide further analysis */
			String height = params.getProperty("height");
			String motifFamilyInfoFile = wd + networkName + clusteringName + "_h" + height+ "_motifFamiliesInfo.tsv" ;

			System.out.println("**Loading representative motifs**");
			HashSet<String> motifSet = loadRepresentativeMotifs(motifFamilyInfoFile);
			System.out.println("Loaded motifs: " + motifSet.size());

			String extractedAnnotationsFile = wd + networkName + clusteringName + "_annotationSubset.tsv";
			String corePorteinsFile = wd + networkName + clusteringName +"_coreProteinsByMotif.tsv";

			String annotatedProteinsPrefix = wd + "Ontologizer/" +  directory2 + allProtDir + networkName + clusteringName + "_annotatedProteinsByMotif_";

			/* Create 1 file per motif with the list of its annotated proteins for Ontologizer analysis */
			System.out.println("Format annotated proteins by motif:");
			formatAnnotatedProteinsByMotif(motifSet, extractedAnnotationsFile, annotatedProteinsPrefix);

			/* Go enrichment of core proteins */
			if(clusteringMeasure == 1 || clusteringMeasure == 2) {

				String coreHeight = params.getProperty("coreHeight");
				String coreMotifFamilyInfoFile = wd + networkName + clusteringName + "_coreProteins_h" + coreHeight+ "_motifFamiliesInfo.tsv" ;

				System.out.println("**Loading representative motifs**");
				motifSet = loadRepresentativeMotifs(coreMotifFamilyInfoFile);
				System.out.println("Loaded motifs: " + motifSet.size());

				String coreProteinsForOntologizerPrefix = wd + "Ontologizer/" + directory2 + coreProtDir +networkName + clusteringName + "_coreProteinsByMotif_";
				/* Create 1 file per motif with the list of annotated core proteins */
				System.out.println("Format core proteins by motif:");
				formatAnnotatedProteinsByMotif(motifSet, corePorteinsFile, coreProteinsForOntologizerPrefix);

			}
		}

		if(Boolean.parseBoolean(params.getProperty("postFunctionalEnrichment"))) {
			/* Load table files - get GO-terms associated to an adjusted p-value < 0.05 */
			String tablePrefix = wd + "Ontologizer/"+ directory2 + allProtDir + "table-"+ networkName + clusteringName + "_annotatedProteinsByMotif_";
			String outputPrefix = wd + "Ontologizer/" + directory2 + allProtDir + networkName + clusteringName + "_significantGOterms_allProteins_";
			
			System.out.println("**Extracting significant GO-terms from all proteins**");
			getSignificantGOterms(tablePrefix, outputPrefix, Integer.parseInt(params.getProperty("motifFamilyGroups")));
			
			if(clusteringMeasure == 1 || clusteringMeasure == 2) {
				tablePrefix = wd + "Ontologizer/"+ directory2 + coreProtDir + "table-"+ networkName + clusteringName + "_coreProteinsByMotif_";
				outputPrefix = wd + "Ontologizer/" + directory2 + coreProtDir + networkName + clusteringName + "_significantGOterms_coreProteins_";
				
				System.out.println("**Extracting significant GO-terms from core proteins**");
				getSignificantGOterms(tablePrefix, outputPrefix, Integer.parseInt(params.getProperty("coreFamilyGroups")));
			}
			 
		}
	}

	private static void formatProteinsInNetworkFile(String proteinAnnotatedFreqFile, String proteinNetworkOutputFile) {

		try {

			InputStream in = new FileInputStream(new File(proteinAnnotatedFreqFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			BufferedWriter out = new BufferedWriter(new FileWriter(new File(proteinNetworkOutputFile)));

			String line = input.readLine();

			while(line!=null) {

				out.write(line.split("\t")[0] + "\n"); // protein name
				out.flush();
				line = input.readLine();
			}
			input.close();
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	/**
	 * Load representative motifs identified in the motif family process
	 *  
	 * @param motifFamilyFile	String - file containing all the representative motifs
	 * @return motifSet			Set<String> - set of representative motifs
	 */
	private static HashSet<String> loadRepresentativeMotifs(String motifFamilyFile){

		HashSet<String> motifSet = new HashSet<>();
		try {

			InputStream in = new FileInputStream(new File(motifFamilyFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // header
			line = input.readLine();

			while(line!=null) {
				motifSet.add(line.split("\t")[1]); // [1] representative motif
				line = input.readLine();
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return motifSet;
	}	

	private static void formatAnnotatedProteinsByMotif(HashSet<String> repMotifSet, String extractedAnnotationFile, String annotatedProteinByMotifPrefix) {

		try {

			InputStream in = new FileInputStream(new File(extractedAnnotationFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // motif \t Prot1|Prot2|Prot3
			int motifCount = 1; 

			while(line!=null && motifCount<= repMotifSet.size()) {

				/* output necessary info when line in file corresponds to representative motif */
				if(repMotifSet.contains(line.split("\t")[0])) {

					System.out.println(motifCount);
					BufferedWriter out = new BufferedWriter(new FileWriter(new File(annotatedProteinByMotifPrefix + motifCount)));
					motifCount++;

					String[] protList = line.split("\t")[2].split("\\|");

					for(String prot: protList) {
						out.write(prot + "\n"); // protein name 
						out.flush();
					}
					out.close();
				}
				line = input.readLine(); // next line 
			}
			input.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private static void getSignificantGOterms(String tableInputFilePrefix, String outputPrefix, int numberOfFiles) {


		for(int i=1; i<= numberOfFiles; i++) {
			HashMap<String, String> goMap = new HashMap<>();
			try {

				InputStream in = new FileInputStream(new File(tableInputFilePrefix + i + "-Parent-Child-Union-Benjamini-Hochberg.txt"));
				BufferedReader input = new BufferedReader(new InputStreamReader(in));

				String line = input.readLine(); // header
				line = input.readLine(); 
				
				while(line!=null) {
					String[] col = line.split("\t"); // [0] = GO-term; [10] = adjusted p-val
					if(Double.parseDouble(col[10]) < 0.05) {
						goMap.put(col[0], col[10]);
					}
					
					line = input.readLine(); // next line 
				}
				input.close();
				
				if(!goMap.isEmpty()) {
					BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputPrefix + i)));
					
					for(Entry<String, String> e: goMap.entrySet()) {
						out.write(e.getKey() + "\t" + e.getValue() + "\n");
						out.flush();
					}
				
					out.close();
				}
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
}
