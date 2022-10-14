import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

public class GenerateCorrelationMatrixes {

	public static void main(String[] args) {

		String corrNetFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/input_files/correlation_v2.tsv";
		String proteinInfo = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/corrNetTop2-400_proteinsInNetwork_info.tsv";

		String formatOutputFile = "/Users/rnadeau2/Documents/LESMoNlocal/figures/network/corrNet2-400_HeatMap.tsv";
		String protOutFile = "/Users/rnadeau2/Documents/LESMoNlocal/figures/network/corrNet2-400_listOfProts.tsv";

		formatCorrelationNetwork(corrNetFile, proteinInfo, formatOutputFile, protOutFile);
	}

	public static void formatCorrelationNetwork(String corrNetFile, String proteinFile, String formatOutput, String protOutput) {

		/* get list of proteins */
		HashSet<String> proteinSet = getListOfProteins(proteinFile);

		List<String> proteinList = new ArrayList<>(proteinSet);

		/* get correlations */ 
		double[][] corrNet = loadCorrelationNetwork(corrNetFile, proteinList);

		printNetwork(corrNet, formatOutput);
		printProteinOrder(proteinList, protOutput);

	}

	/* load list of proteins in the final network */
	private static HashSet<String> getListOfProteins(String proteinInfoFile){

		HashSet<String> proteinSet = new HashSet<>();

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(proteinInfoFile))));

			String line = in.readLine();
			while(line!=null) {

				String[] col = line.split("\t");
				proteinSet.add(col[0]);

				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return proteinSet;
	}

	private static double[][] loadCorrelationNetwork(String corrNetFile, List<String> proteinList){

		double[][] corrNet = new double[proteinList.size()][proteinList.size()];

		/* index protein list */
		HashMap<String, Integer> proteinIdxMap = new HashMap<>();
		for(int i=0; i<proteinList.size(); i++) {
			proteinIdxMap.put(proteinList.get(i), i);
		}

		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(corrNetFile))));

			String line = in.readLine(); // header
			line = in.readLine();
			while(line!=null) {

				String[] col = line.split("\t");

				/* store correlations if they were kept in network */
				if(proteinIdxMap.containsKey(col[0]) && proteinIdxMap.containsKey(col[1])) {

					double corr = Double.parseDouble(col[2]);
					if(corr >=  0.574171) {
						corrNet[proteinIdxMap.get(col[0])][proteinIdxMap.get(col[1])] = corr;
						corrNet[proteinIdxMap.get(col[1])][proteinIdxMap.get(col[0])] = corr;
					}
				}
				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return corrNet;
	}

	private static void printNetwork(double[][] corrNet, String outputFile) {

		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			for(int i=0; i<corrNet.length; i++) {
				for(int j=0; j<corrNet.length; j++) {
					out.write(corrNet[i][j] + "\t");
				}
				out.write("\n");
				out.flush();
			}

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	private static void printProteinOrder(List<String> proteinList, String outputFile) {
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFile)));

			for(int i=0; i<proteinList.size(); i++) {
				out.write(proteinList.get(i) + "\n");
				out.flush();
			}

			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
