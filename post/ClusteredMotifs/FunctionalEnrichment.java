package ClusteredMotifs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

public class FunctionalEnrichment {

	public static void formatFilesForOntologizer(String proteinAnnotatedFreqFile, String extractedAnnotationFile, String proteinNetworkOutFile,	String annotatedProteinByMotifPrefix) {
		
		/* Create list of proteins in network from protein frequency file */
		System.out.println("Format proteins in network");
		formatProteinsInNetworkFile(proteinAnnotatedFreqFile, proteinNetworkOutFile);
		
		/* Create 1 file per motif with the list of its annotated proteins */
		System.out.println("Format annotated proteins by motif:");
		formatAnnotatedProteinsByMotif(extractedAnnotationFile, annotatedProteinByMotifPrefix);
	}
	
	public static void formatCoreFilesForOntologizer(String coreProteinAnnotationFile, String outputFilePrefix) {
		
		/* Create 1 file per motif with the list of annotated core proteins */
		System.out.println("Format core proteins by motif:");
		formatAnnotatedProteinsByMotif(coreProteinAnnotationFile, outputFilePrefix);
		
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
	
	private static void formatAnnotatedProteinsByMotif(String extractedAnnotationFile, String annotatedProteinByMotifPrefix) {
		
		try {

			InputStream in = new FileInputStream(new File(extractedAnnotationFile));
			BufferedReader input = new BufferedReader(new InputStreamReader(in));
			
			String line = input.readLine(); // motif \t Prot1|Prot2|Prot3
			int motifCount = 1; 
			
			while(line!=null) {
				
				System.out.println(motifCount);
				BufferedWriter out = new BufferedWriter(new FileWriter(new File(annotatedProteinByMotifPrefix + motifCount)));

				String[] protList = line.split("\t")[2].split("\\|");
 				
				for(String prot: protList) {
					out.write(prot + "\n"); // protein name
					out.flush();
				}
				
				line = input.readLine();
				motifCount++;
				
				out.close();
			}
			input.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
}
