package positionConservation;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Properties;

public class AssessMotifPositionConservation {

	public static void main(String[] args) throws FileNotFoundException, IOException {

		System.out.println("**Loading parameters file** \n");
		Properties params = new Properties();
		params.load(new FileInputStream(args[0]));

		String wd = params.getProperty("working_directory");

		File directory = new File(wd + "/MotifPosition/");
		if (! directory.exists()){
			System.out.println("creating directory: MotifPosition/");
			directory.mkdir();
		}

		String fastaFile = wd + params.getProperty("fastaFile");
		
		String motifFamilies = wd + args[3];
		String motifOutputPrefixFile = wd + args[4];

		PositionConservation p = new PositionConservation(fastaFile, 8, Integer.parseInt(args[1]));
		p.getMotifPositionsFromLongestSequences(motifFamilies, fastaFile, motifOutputPrefixFile, Integer.parseInt(args[5]));


	}

}
