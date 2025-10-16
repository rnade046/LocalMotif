package opt;

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
import java.util.Map.Entry;
import java.util.Set;

public class CompareMotifAnnotations {

	/* args[0] = annotation directory 1
	 * args[1] = annotation directory 2
	 * args[3] = file number
	 */
	public static void main(String[] args) {

		int File = Integer.parseInt(args[4]);

		String annotationFile1 = args[0] + "annotation_" + File + ".tsv";
		String annotationFile2 = args[1] + "annotation_" + File + ".tsv";
		String annotatedProteinsFile = args[2];
		String wd = args[3];
		
		String degreeFile = wd + "corrNetTop2-400_degreesInNetwork.tsv";
		String proteinChangesFile = wd + "annotationChanges/ProteinAnnotationChanges_" + File;

		HashSet<String> proteinSet = loadProteinsInNetwork(annotatedProteinsFile);
		List<HashSet<String>> motifSets = initializeMotifsInBatches(annotationFile1);
		HashMap<String, int[]> proteinModifications = initializeProteinModsMap(proteinSet);
		
		List<String> output = new ArrayList<>();

		for(HashSet<String> motifs : motifSets) {

			/* map {key = motif ; value = Set<annotated proteins> } */
			HashMap<String, HashSet<String>> annotations1 = getAnnotationsFromSet1(annotationFile1, proteinSet, motifs);
			HashMap<String, HashSet<String>> annotations2 = getAnnotationsFromSet1(annotationFile2, proteinSet, motifs);

			for(Entry<String, HashSet<String>> currentMotif : annotations1.entrySet()) {

				String motif = currentMotif.getKey();
				if(annotations2.containsKey(motif)){
					
					HashSet<String> annotatedProteins1 = currentMotif.getValue();
					HashSet<String> annotatedProteins2 = annotations2.get(motif);
					
					double overlap = measureProteinOverlap(annotatedProteins1, annotatedProteins2);
					double diff = annotatedProteins1.size() - annotatedProteins2.size();
					proteinModifications = evaluateProteinChanges(annotatedProteins1, annotatedProteins2, proteinModifications);
					
					output.add(currentMotif.getKey() + "\t" + currentMotif.getValue().size() + "\t" + annotations2.get(currentMotif.getKey()).size() + "\t" + overlap + "\t" + diff + "\n");
				}
			}
		}
		printProteinChanges(proteinModifications, degreeFile, proteinChangesFile);
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(wd + "/compareAnnotations/comparedAnnotations_" + File)));

			for(String motifInfo : output) {
				out.write(motifInfo);
				out.flush();	
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private static List<HashSet<String>> initializeMotifsInBatches(String annotationFile) {

		/* obtain original list of motifs */
		List<String> motifList = new ArrayList<>();
		try {
			FileInputStream in = new FileInputStream(new File(annotationFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String motif;
			while ((motif = input.readLine()) != null) {

				motifList.add(motif.split("\t")[0]);
				motif = input.readLine();
			}
			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		/* batch motifs into chunks of 1000 for sequential annotations */
		List<HashSet<String>> motifListBatches = new ArrayList<>();
		int batchSize=1000;

		for(int i=0; i<motifList.size(); i+=batchSize) {
			int endIdx = Math.min(i + batchSize, motifList.size());
			motifListBatches.add(new HashSet<>(motifList.subList(i, endIdx)));
		}

		return motifListBatches;
	}


	private static HashMap<String, HashSet<String>> getAnnotationsFromSet1(String annotationFile, HashSet<String> proteinSet, HashSet<String> motifs){

		HashMap<String, HashSet<String>> motifAnnotations = new HashMap<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(annotationFile))));
			String line = in.readLine(); 

			while((line = in.readLine()) != null){

				String[] col=line.split("\t");

				/* load annotations for current batch of motifs */
				if(motifs.contains(col[0])) {
					motifAnnotations.put(col[0], keepProteinsInNetwork(col[2].split("\\|"), proteinSet));
				}
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return motifAnnotations;
	}

	private static HashSet<String> loadProteinsInNetwork(String file){

		HashSet<String> proteinSet = new HashSet<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(file))));
			String line; 

			while((line = in.readLine()) != null){

				String[] col=line.split("\t");
				if(col.length > 1) {
					proteinSet.add(col[0]);
				}
			}
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return proteinSet;
	}

	private static HashSet<String> keepProteinsInNetwork(String[] annotatedProteins, HashSet<String> proteinSetInNetwork){

		HashSet<String> finalAnnotatedProteins = new HashSet<>();

		for(String prot: annotatedProteins) {
			if(proteinSetInNetwork.contains(prot)) {
				finalAnnotatedProteins.add(prot);
			}
		}
		return finalAnnotatedProteins;
	}

	private static double measureProteinOverlap(HashSet<String> annotations1, HashSet<String> annotations2) {

		int overlap = 0;
		for(String protein: annotations1) {
			if(annotations2.contains(protein)) {
				overlap++;
			}
		}

		double overlap_norm = overlap / (double) Math.max(annotations1.size(), annotations2.size());
		return overlap_norm;
	}

	private static HashMap<String, Integer> loadMap(String inputFile) {

		HashMap<String, Integer> map = new HashMap<>();

		InputStream in;
		try {
			in = new FileInputStream(new File(inputFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));

			String line = input.readLine(); // header
			line = input.readLine();
			while(line!=null) {

				String[] col = line.split("\t"); // [0] = protein Name, [1] = integer (#degrees or #motif)
				map.put(col[0], Integer.parseInt(col[1]));

				line = input.readLine();
			}

			input.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return map;
	}

	private static HashMap<String, int[]> initializeProteinModsMap(Set<String> proteinNames){
		
		HashMap<String, int[]> proteinModifications = new HashMap<>();

		for(String p: proteinNames) {
			proteinModifications.put(p, new int[2]);
		}
		return proteinModifications;
	}
	
	private static HashMap<String, int[]> evaluateProteinChanges(Set<String> annotatedProteins1, Set<String> annotatedProteins2, HashMap<String, int[]> proteinModifications){
		
		// removed proteins in Set2
		for(String p: annotatedProteins1) {
			if(!annotatedProteins2.contains(p)) {
				int[] mods = proteinModifications.get(p);
				mods[0]++;
				proteinModifications.put(p, mods); 
			}
		}
		
		// added proteins in Set2
		for(String p: annotatedProteins2) { 
			if(!annotatedProteins1.contains(p)) {
				int[] mods = proteinModifications.get(p);
				mods[1]++;
				proteinModifications.put(p, mods); 
			}
		}
		return proteinModifications;
	}
	
	private static void printProteinChanges(HashMap<String, int[]> proteinModifications, String degreeFile, String output) {
		
		HashMap<String, Integer> proteinDegreeMap = loadMap(degreeFile);
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(output)));
			out.write("ProteinName\t#Degrees\t#ProteinsRemoved\t#ProteinsAdded");
			
			for(Entry<String, int[]> proteinEntry : proteinModifications.entrySet()) {
				
				out.write(proteinEntry.getKey() + "\t");
				
				if(proteinDegreeMap.containsKey(proteinEntry.getKey())){
					out.write(proteinDegreeMap.get(proteinEntry.getKey()) + "\t");
				} else {
					out.write("NA\t");
				}
				out.write(proteinEntry.getValue()[0] + "\t" + proteinEntry.getValue()[1] + "\n");
				out.flush();
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
