package formatHCM;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;

public class FormatHumanCellMap {

	public static void main(String[] args) {

		String inputNetworkFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/input_files/correlation_v2.tsv";
		double percentThreshold = 0.02;
		int maxInteractors = 400;

		String outputNetwork = "/Users/rnadeau2/Documents/Structures/hcm/corrNet2-400_formattedNetwork_test.tsv";
		String removedInteractors = "/Users/rnadeau2/Documents/Structures/hcm/corrNet2-400_removedInteractors_test.tsv";
		/* load full HCM network */
		List<Interaction> hcm = loadHCM(inputNetworkFile);
		System.out.println("Original network - interactions: " + hcm.size() + "\n");

		/* filter interactions for top percent interactions */
		hcm = assessPercentageScore(hcm, percentThreshold);
		System.out.println("Top X% network - interactions: " + hcm.size());
		System.out.println("Top X% - predicted interactions: " + Math.round(hcm.size() * percentThreshold) + "\n");

		/* determine overly connected interactors */
		HashSet<String> interactorsToRemove = determineOverConnectedInteractors(hcm, maxInteractors);
		System.out.println("Removed interactors: " + interactorsToRemove.size() + "\n");
		/* remove over connected interactors */
		hcm = removeOverConnectedInteractors(hcm, interactorsToRemove);
		System.out.println("TopX - max interactions network - interactions: " + hcm.size());
		/* print updated network */
		printRemovedInteractions(interactorsToRemove, removedInteractors);
		printUpdatedNetwork(hcm, outputNetwork);
	}

	public static ArrayList<Interaction> loadHCM(String inputNetworkFile) {

		ArrayList<Interaction> hcmInteractions = new ArrayList<>();

		try {
			BufferedReader in = new BufferedReader(
					new InputStreamReader(new FileInputStream(new File(inputNetworkFile))));

			String line = in.readLine(); // header
			line = in.readLine();

			while (line != null) {

				String[] elements = line.split("\t");
				hcmInteractions.add(new Interaction(elements[0], elements[1], Double.parseDouble(elements[2])));
				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return hcmInteractions;
	}

	public static List<Interaction> assessPercentageScore(List<Interaction> hcm, double percentThreshold) {

		ArrayList<Interaction> updatedHCM = new ArrayList<>();

		/* determine number of interactions */
		int numInteractions = (int) Math.round(hcm.size() * percentThreshold);

		/* make list of correlation score and sort */
		List<Double> scores = new ArrayList<>();
		for (Interaction inter : hcm) {
			scores.add(inter.getWeight());
		}

		Collections.sort(scores);
		Collections.reverse(scores);

		double minimumScore = scores.get(numInteractions);

		for (Interaction inter : hcm) {
			if (inter.getWeight() >= minimumScore) {
				updatedHCM.add(inter);
			}
		}
		return updatedHCM;
	}

	public static HashSet<String> determineOverConnectedInteractors(List<Interaction> hcm, int threshold) {

		/* count number of interactions */
		HashMap<String, Integer> proteinInteractionCount = new HashMap<>();

		for (Interaction inter : hcm) {
			if (proteinInteractionCount.containsKey(inter.getProtein1())) {
				proteinInteractionCount.put(inter.getProtein1(), proteinInteractionCount.get(inter.getProtein1()) + 1);
			} else {
				proteinInteractionCount.put(inter.getProtein1(), 1);
			}

			if (proteinInteractionCount.containsKey(inter.getProtein2())) {
				proteinInteractionCount.put(inter.getProtein2(), proteinInteractionCount.get(inter.getProtein2()) + 1);
			} else {
				proteinInteractionCount.put(inter.getProtein2(), 1);
			}
		}

		/* determine interactors with more than required count */
		HashSet<String> interactorsToRemove = new HashSet<>();
		for (Entry<String, Integer> entry : proteinInteractionCount.entrySet()) {

			if (entry.getValue() >= threshold) { // number interactions
				interactorsToRemove.add(entry.getKey());
			}
		}
		return interactorsToRemove;
	}

	public static List<Interaction> removeOverConnectedInteractors(List<Interaction> hcm,
			HashSet<String> interactorsToRemove) {

		List<Interaction> updatedHCM = new ArrayList<>();

		for (Interaction inter : hcm) {
			if (!interactorsToRemove.contains(inter.getProtein1())
					|| !interactorsToRemove.contains(inter.getProtein2())) {
				updatedHCM.add(inter);
			}
		}

		return updatedHCM;
	}

	public static void printUpdatedNetwork(List<Interaction> hcm, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			out.write("Protein1\tProtein2\tScore\n");
			for (Interaction inter : hcm) {
				out.write(inter.getProtein1() + "\t" + inter.getProtein2() + "\t" + inter.getWeight() + "\n");
				out.flush();
			}

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static void printRemovedInteractions(HashSet<String> interactorsToRemove, String outputFile) {
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			for (String prot : interactorsToRemove) {
				out.write(prot + "\n");
				out.flush();
			}

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
