public class RandomMain {
	
	public static void main(String[] args) {
		
		/* .fasta containing sequences to randomize */
		String fastaFile = args[0];
		
		/* mapping file {proteinName = refSeqId1|refSeqId2|...|} */
		String rnaIdListFile = args[1];
		
		/* output file */
		String randomFastaFile = args[0] + "_winShuffled.fasta";
		
		System.out.println("**Generating randomized fasta sequences**");
		RandomizeSequences.generateRandomizedFasta(rnaIdListFile, fastaFile, randomFastaFile);
	}


}
