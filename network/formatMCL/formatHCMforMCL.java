package formatMCL;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashSet;

public class formatHCMforMCL {

	public static void main(String[] args) {

		String hcmNetworkFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/hcm/corrNet2-400_formattedNetwork.tsv";
		String proteinInfoFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/corrNetTop2-400_proteinsInNetwork_info.tsv";

		String mclOutputFile = "/Users/rnadeau2/Documents/LESMoNlocal/analysis/hcm/corrNet2-400_mclNetwork_corrCoef.tsv";

		/* load list of proteins in final network */
		HashSet<String> proteinSet = getProteinsInNetwork(proteinInfoFile);

		/* iterate through HCM network and output to ABC format with correlation scores */
		try {
			BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(new File(hcmNetworkFile))));
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(mclOutputFile)));

			String line = in.readLine();
			while(line != null) {

				String[] values = line.split("\t"); // [0] = protein1, [1] = protein2, [3] = correlation coefficient
				if(proteinSet.contains(values[0]) && proteinSet.contains(values[1])) {
					if(!values[0].equals(values[1])) {
						out.write(values[0] + " " + values[1] + " " + values[2] + "\n");
						out.flush();
					}
				}
				line = in.readLine();
			}
			out.close();
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	public static HashSet<String> getProteinsInNetwork(String proteinInfoFile){

		HashSet<String> proteinSet = new HashSet<>();

		BufferedReader in;
		try {
			in = new BufferedReader(new InputStreamReader(new FileInputStream(proteinInfoFile)));

			String line = in.readLine();

			while(line!=null) {

				proteinSet.add(line.split("\t")[0]);
				line = in.readLine();
			}

			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return proteinSet;
	}
}
