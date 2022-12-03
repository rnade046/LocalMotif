package sequenceLength;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

public class Motif {

	
	String motif;
	int order;
	HashSet<String> proteins;
	HashSet<String> refSeqIds;
	List<Integer> sequenceLengths;
	
	
	public Motif(String name, int position) {
		motif = name;
		order = position;
	}
	
	public void setProteins(HashSet<String> proteinSet) {
		this.proteins = proteinSet;
	}
	
	public void setRefSeqIds(HashMap<String, String[]> proteinToIdMap) {
		HashSet<String> idSet = new HashSet<>();
		
		for(String prot : this.proteins) {
			
			if(proteinToIdMap.containsKey(prot)) {
				for(String id: proteinToIdMap.get(prot)) {
					idSet.add(id);
				}
				
			}
		}
		
		this.refSeqIds = idSet;
	}
	
	public void  setSequenceLenghts(List<Integer> lengths) {
		this.sequenceLengths = lengths;
	}
	
	public String getMotifName() {
		return this.motif;
	}
	
	public int getMotifOrder() {
		return this.order;
	}
	
	public HashSet<String> getProteins() {
		return this.proteins;
	}
	
	public HashSet<String> getIds(){
		return this.refSeqIds;
	}
	
	public List<Integer> getLengths(){
		return this.sequenceLengths;
	}
	
}
