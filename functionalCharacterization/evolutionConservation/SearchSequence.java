package evolutionConservation;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;

public class SearchSequence {

	private String sequence;
	private String refSeqId;
	private boolean positive;
	private String chromosome;
	private Integer startPosition;
	private Integer endPosition;
	private List<Integer> searchPositions;
	private HashSet<Integer> foundPositions;

	// constructor
	public SearchSequence(String seq, String header) {

		this.sequence = seq;
		this.refSeqId = header.split("[\\_\\s++\\.]")[2] + "_"+ header.split("[\\_\\s++\\.]")[3];
		this.positive = header.contains("+");
		this.chromosome = header.split(" ")[1].split("=")[1].split(":")[0];
		this.startPosition =  Integer.parseInt(header.split(" ")[1].split("=")[1].split(":")[1].split("-")[0]);
		this.endPosition =  Integer.parseInt(header.split(" ")[1].split("=")[1].split(":")[1].split("-")[1]);
		this.searchPositions = setList(seq);
		this.foundPositions = new HashSet<>();

	}	

	public String getSequence() {
		return sequence;
	}

	public String getId() {
		return refSeqId;
	}

	public boolean isPositive() {
		return positive;
	}

	public String getChromosome() {
		return chromosome;
	}

	public int getStartPosition() {
		return startPosition;
	}

	public int getEndPosition() {
		return endPosition;
	}

	public List<Integer> setList(String seq){
		List<Integer> positionsToSearch = new ArrayList<>();
		for(int i=0; i<seq.length(); i++) {
			positionsToSearch.add(i);
		}
		return positionsToSearch;
	}

	public List<Integer> getPositionsToSearch(){
		return searchPositions;
	}

	public void setPositionsToSearch(List<Integer> positions) {
		this.searchPositions = positions;
	}

	public HashSet<Integer> getFoundPositions(){
		return this.foundPositions;
	}

	public void setFoundPosition(int pos) {		
		for(int i=pos; i<pos + 8; i++) {
			this.foundPositions.add(i);
		}
	}

	public void setPosition(List<Integer> positions) {
		this.foundPositions.addAll(positions);
	}

	public void removePositionToSearch(List<Integer> idx) {
		this.searchPositions.removeAll(idx);
	}

	public void updatePositionsToSearch() {

		HashSet<Integer> positionsToKeep = new HashSet<>();	
		
		List<Integer> foundPositionsCopy = new ArrayList<>();
		foundPositionsCopy.addAll(this.foundPositions);
		Collections.sort(foundPositionsCopy);
		
		if(!this.searchPositions.isEmpty()) {
			/* for each positions remaining to search; check if it is contained in the coverage of motifs */
			for(int pos : this.searchPositions) {
				/* if position was not found to be covered */
				if(!foundPositions.contains(pos)) {
					/* for positions remained uncovered; add itself +/- 7 indexes but only if they are not in the current position to search */
					for(int i=pos-7; i<= pos + 7; i++) {
						if(this.searchPositions.contains(i)) {
							positionsToKeep.add(i);
						}
					}

				}
			}

			/* update list of remaining positions */
			List<Integer> newPositions = new ArrayList<>();
			newPositions.addAll(positionsToKeep);
			Collections.sort(newPositions);
			this.searchPositions = newPositions;
		}
	}
}
