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

		if(Boolean.parseBoolean(params.getProperty("summarizeCellularComponents"))) {
			
			int numFamilies = Integer.parseInt(params.getProperty("motifFamilyGroups"));
			String revigoPrefix = wd + "Ontologizer/" + directory2 + allProtDir + networkName + clusteringName + "_revigo_allProteins_";
			String significantPrefix = wd + "Ontologizer/" + directory2 + allProtDir + networkName + clusteringName + "_significantGOterms_allProteins_";
			String outputFile = wd + "Ontologizer/" + directory2 + allProtDir + networkName + clusteringName + "_allProteins_cellularComponentSummary_matrix.tsv";
			
			/* Get significant cellular components */
			System.out.println("**Obtain significant cellular components - all proteins**");
			HashMap<String, String> ccInfoMap = getAllCellularComponents(revigoPrefix, numFamilies);
			
			/* Summarize cellular components & print */
			System.out.println("**Summarize cellular compartment significance scores - all proteins**");
			HashMap<String, Double[]> significantCCMap = summarizeCellularComponents(ccInfoMap, significantPrefix, numFamilies);
			printCellularComponentSummary(ccInfoMap, significantCCMap, numFamilies, outputFile);
			
			if(clusteringMeasure == 1 || clusteringMeasure == 2) {
		
				numFamilies = Integer.parseInt(params.getProperty("coreFamilyGroups"));
				revigoPrefix = wd + "Ontologizer/" + directory2 + coreProtDir + networkName + clusteringName + "_revigo_coreProteins_";
				significantPrefix = wd + "Ontologizer/" + directory2 + coreProtDir + networkName + clusteringName + "_significantGOterms_coreProteins_";
				outputFile = wd + "Ontologizer/" + directory2 + coreProtDir + networkName + clusteringName + "_coreProteins_cellularComponentSummary_matrix.tsv";
				
				/* Get significant cellular components */
				System.out.println("**Obtain significant cellular components - core proteins**");
				ccInfoMap = getAllCellularComponents(revigoPrefix, numFamilies);
				System.out.println("ccMap : " + ccInfoMap.size());
				/* Summarize cellular components & print */
				System.out.println("**Summarize cellular compartment significance scores - core proteins**");
				significantCCMap = summarizeCellularComponents(ccInfoMap, significantPrefix, numFamilies);
				System.out.println("summaryMap : " + significantCCMap.size());
				printCellularComponentSummary(ccInfoMap, significantCCMap, numFamilies, outputFile);
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

	private static HashMap<String, String> getAllCellularComponents(String revigoPrefixFile, int numFamilies) {

		HashMap<String, String> ccInfo = new HashMap<>(); //GO:XX, Name

		for(int i=1; i<= numFamilies; i++) {
			System.out.print(i + ".");
			try {

				InputStream in = new FileInputStream(new File(revigoPrefixFile+ i));
				BufferedReader input = new BufferedReader(new InputStreamReader(in));

				String line = input.readLine(); // header
				for(int lineCount=1; lineCount<=5; lineCount ++) {
					line = input.readLine(); 
				}

				while(line!=null) {

					String[] col = line.split(","); // [0] = GO-term; [1] = name
					String term = col[0].split("\"")[1];
					String name = col[1].split("\"")[1];
					if(!ccInfo.containsKey(term)) {
						ccInfo.put(term, name);
					}
					line = input.readLine(); // next line 
				}
				input.close();

			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		System.out.println("Done");
		return ccInfo;
	}		

	private static HashMap<String, Double[]> summarizeCellularComponents(HashMap<String, String> ccInfoMap, String significantGOsPrefixFile, int numberFamily) {

		HashMap<String, Double[]> significantCCMap = new HashMap<>(); //GO:XXD, Double[family1, 2, .., x] (adjusted p-value)

		for(int i=1; i<=numberFamily; i++) {
			System.out.print(i + ".");
			try {

				InputStream in = new FileInputStream(new File(significantGOsPrefixFile+ i));
				BufferedReader input = new BufferedReader(new InputStreamReader(in));

				String line = input.readLine(); // header

				while(line!=null) {
					String[] col = line.split("\t"); // [0] = GO-term; [1] = adjusted p-val
					
					/* store info if significant GO:XX is a CC */
					if(ccInfoMap.containsKey(col[0])) {

						/* create entry if it's the first time this GO:XX is seen */ 
						if(!significantCCMap.containsKey(col[0])) {
							significantCCMap.put(col[0], new Double[numberFamily]);
						}
						/* update adjusted p-value at appropriate position */
						Double[] array = significantCCMap.get(col[0]);
						array[i-1] = Double.parseDouble(col[1]);
						
						significantCCMap.put(col[0], array);
					}
					line = input.readLine(); // next line 
				}
				input.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		System.out.println("Done");
		return significantCCMap;
	}

	private static void printCellularComponentSummary(HashMap<String, String> ccInfoMap, HashMap<String, Double[]> significantCCMap, int numFamilies, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));
			
			// header 
			out.write("GO-term\tGO-name\t#Groups\t");
			for(int i=1; i<= numFamilies; i++) {
				out.write(i + "\t");
			}
			out.write("\n");
			
			for(Entry<String, Double[]> goEntry: significantCCMap.entrySet()) {
				String go = goEntry.getKey();
				Double[] scores = goEntry.getValue();
				
				int count = 0;
				for(int i=0; i<scores.length; i++) {
					if(scores[i] != null) {
						count++;
					}
				}
				
				out.write(go + "\t" + ccInfoMap.get(go) + "\t" + count + "\t");
				for(int i=0; i<scores.length; i++) {
					out.write(scores[i] + "\t");
				}
				out.write("\n");
				out.flush();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

}
