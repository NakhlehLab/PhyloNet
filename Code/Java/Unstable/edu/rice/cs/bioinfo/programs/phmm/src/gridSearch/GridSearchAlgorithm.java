package gridSearch;

import java.util.ArrayList;
import java.util.Map;
import java.util.Set;

import phylogeny.EvoTree;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.Observation;
import be.ac.ulg.montefiore.run.jahmm.phmm.HiddenState;
import be.ac.ulg.montefiore.run.jahmm.phmm.TransitionProbabilityParameters;


public class GridSearchAlgorithm<O extends Observation> {

	private int gBranch;
	private int gRecombination;
	private int gHybridization;
	private int gBaseSub;
	
	private double branchMin;
	private double branchMax;
	private double recombinationMin;
	private double recombinationMax;
	private double hybridizationMin;
	private double hybridizationMax;
	private double baseSubMin;
	private double baseSubMax;
	
	private ArrayList<Nob> nobs;
	
	public GridSearchAlgorithm(int gBranchIn, int gRecombinationIn, int gHybridizationIn, int gBaseSubIn,
			double branchMinIn, double branchMaxIn, double recombinationMinIn,
			double recombinationMaxIn, double hybridizationMinIn, double hybridizationMaxIn, 
			double baseSubMinIn, double baseSubMaxIn) {
		
		this.gBranch = gBranchIn;
		this.gRecombination = gRecombinationIn;
		this.gHybridization = gHybridizationIn;
		this.gBaseSub = gBaseSubIn;
		
		this.branchMin = branchMinIn;
		this.branchMax = branchMaxIn;
		this.recombinationMin = recombinationMinIn;
		this.recombinationMax = recombinationMaxIn;
		this.hybridizationMin = hybridizationMinIn;
		this.hybridizationMax = hybridizationMaxIn;
		this.baseSubMin = baseSubMinIn;
		this.baseSubMax = baseSubMaxIn;
		
		this.nobs = new ArrayList<Nob>();
		
	}
	
	public void runGridSearch(Hmm<O> hmm, TransitionProbabilityParameters tpp, O observation, ArrayList<HiddenState> trees_states, Map<EvoTree,Set<HiddenState>> parentalTreeClasses) {
		
		//add base sub nob
		nobs.add(new BaseSubNob(gBaseSub, branchMin, branchMax));
		
		//add transition probability recombination
		nobs.add(new RecombinationFreqNob(gRecombination, recombinationMin, recombinationMax, hmm, tpp, trees_states, parentalTreeClasses));
		
		//add transition probability hybridization
		nobs.add(new RecombinationFreqNob(gHybridization, hybridizationMin, hybridizationMax, hmm, tpp, trees_states, parentalTreeClasses));
		
		//add tree branches
		int currentparent = -1;
		for (int i = 0; i < hmm.nbStates(); i++) {
			// get the parent tree
			// get id
			//if (parenttreeid != currentparent) {
				// do something
				//currentparent = parenttreeid;
				
			//}
			// get geneology tree
			
		}
		
		
	}
	
	


}
