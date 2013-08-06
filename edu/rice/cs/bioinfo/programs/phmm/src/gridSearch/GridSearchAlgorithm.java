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
		nobs.add(new RecombinationFreqNob(gRecombination, recombinationMin, recombinationMax, tpp));
		
		//add transition probability hybridization
		nobs.add(new RecombinationFreqNob(gHybridization, hybridizationMin, hybridizationMax, tpp));
		
		//add tree branches
		
		
	}
	
	

	
	
	
	 /**
     * Calculate initial transition probability matrix a_{ij}.
     * See revised writeup for details.
     *
     * Call this after changing parental tree branch lengths
     */
    private double[][] calculateAij (ArrayList<HiddenState> trees_states, double recombinationFreq, double hybridizationFreq, Map<EvoTree,Set<HiddenState>> parentalTreeClasses) {
	double[][] a = new double[trees_states.size()][trees_states.size()];
	for (int i = 0; i < a.length; i++) {
	    HiddenState si = trees_states.get(i);
	    double totalNonSelfTransitionProbabilities = 0.0;
	    for (int j = 0 ; j < a[i].length; j++) {
		// set self-transition probability at the end
		if (i == j) {
		    continue;
		}
		
		HiddenState sj = trees_states.get(j);
		a[i][j] = sj.calculateProbabilityOfGeneGenealogyInParentalTree();
		if (checkSameParentalClass(si, sj, parentalTreeClasses)) {
		    a[i][j] *= recombinationFreq;
		}
		else {
		    a[i][j] *= hybridizationFreq;
		}

		totalNonSelfTransitionProbabilities += a[i][j];
	    }
	    
	    // now set self-transition probability
	    a[i][i] = 1.0 - totalNonSelfTransitionProbabilities;
	}
	
	// strict!
	if (!verifyAij(a)) {
	    System.err.println ("ERROR: verifyAij() failed. Returning null to signal error.");
	    return (null);
	}

	return (a);
    }
	
    
    protected boolean checkSameParentalClass (HiddenState si, HiddenState sj, Map<EvoTree,Set<HiddenState>> parentalTreeClasses) {
    	Set<HiddenState> sic = parentalTreeClasses.get(si.getParentalTree());
    	return (sic.contains(sj));
        }
	
    /**
     * By construction, rows of a_ij matrix sum to one.
     */
    protected boolean verifyAij (double[][] a) {
	for (int i = 0; i < a.length; i++) {
	    for (int j = 0; j < a[i].length; j++) {
		if ((a[i][j] < 0.0) || (a[i][j] > 1.0)) {
		    System.err.println ("ERROR: entry in a_ij transition matrix is invalid. " + a[i][j]);
		    return (false);
		}
	    }
	}

	return (true);
    }
    
}
