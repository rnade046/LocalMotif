package randomSeq;

public class RandomMain {
	
	public static void main(String[] args) {
		
		String wd = "C:\\Users\\rnade046\\Documents\\LESMoNlocal\\analysis\\";
		
		String fastaFile = wd + "input_files\\human_3UTRsequences.txt";
		String rnaIdListFile = wd + "motif_enumeration\\BiomaRt_MappingRefSeqIdsToGeneSymbol_corrNet.tsv";
		String randomFastaFile = wd + "input_files\\random_humanCellMap_3UTRsequences.txt";
		
		System.out.println("**Generating randomized fasta sequences**");
		RandomizeSequences.generateRandomizedFasta(rnaIdListFile, fastaFile, randomFastaFile);
	}


}
