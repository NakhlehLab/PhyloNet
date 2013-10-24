package gridSearch;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import runHmm.runHmm;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.phmm.HiddenState;
import be.ac.ulg.montefiore.run.jahmm.phmm.ObservationMap;
import be.ac.ulg.montefiore.run.jahmm.phmm.TransitionProbabilityParameters;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbabilityYF.CoalescePattern;



public class GridSearchAlgorithm {

    private int gBranch;
    //private int gRecombination;
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

    protected runHmm runHmmObject;

    private ArrayList<Nob> nobs;

    // int gRecombinationIn,
    //             double recombinationMinIn, double recombinationMaxIn,
    public GridSearchAlgorithm(int gBranchIn, 
            int gHybridizationIn, int gBaseSubIn,
            double branchMinIn, double branchMaxIn,
            double hybridizationMinIn, double hybridizationMaxIn,
			       double baseSubMinIn, double baseSubMaxIn,
			       runHmm inRunHmmObject
			       ) {

        this.gBranch = gBranchIn;
        //this.gRecombination = gRecombinationIn;
        this.gHybridization = gHybridizationIn;
        this.gBaseSub = gBaseSubIn;

        this.branchMin = branchMinIn;
        this.branchMax = branchMaxIn;
        // this.recombinationMin = recombinationMinIn;
        // this.recombinationMax = recombinationMaxIn;
        this.hybridizationMin = hybridizationMinIn;
        this.hybridizationMax = hybridizationMaxIn;
        this.baseSubMin = baseSubMinIn;
        this.baseSubMax = baseSubMaxIn;

	runHmmObject = inRunHmmObject;
    }

    // throws CloneNotSupportedException

    /**
     * Builds the Nobs array that contain all parameters
     * that need to be learned
     * GridSearch Algorithm will later traverse through this array in
     * order to find the optimal value for each parameter
     * @param hmm
     * @param tpp
     * @param observation
     * @param trees_states
     * @param parentalTreeClasses
     * @throws CloneNotSupportedException
     */
    private void initializeGridSearch(Hmm<ObservationMap> hmm,
            TransitionProbabilityParameters tpp,
            ArrayList<HiddenState> trees_states,
            Map<Network<CoalescePattern[]>,Set<HiddenState>> parentalTreeClasses) {
        // Initialize Nobs array
        nobs = new ArrayList<Nob>();

        //add base sub nob
        nobs.add(new BaseSubNob(gBaseSub, baseSubMin, baseSubMax));

        //add transition probability recombination
        //nobs.add(new RecombinationFreqNob(gRecombination, recombinationMin,
	//recombinationMax, hmm, tpp, trees_states, parentalTreeClasses, runHmmObject));

        //add transition probability hybridization
        nobs.add(new HybridizationFreqNob(gHybridization, hybridizationMin,
					  hybridizationMax, hmm, tpp, trees_states, parentalTreeClasses, runHmmObject));

    	// kliu - need to use PhyloNet tree/network structures
    	// walk through unique parental tree objects
    	for (Network<CoalescePattern[]> parentalTree : parentalTreeClasses.keySet()) {
    	    for (NetNode<CoalescePattern[]> node : parentalTree.dfs()) {
    		// root has no incoming edge to optimize
    		if (node.isRoot()) {
    		    continue;
    		}

    		// create a knob
    		nobs.add(new ParentalTreeNob(gBranch, branchMin, branchMax, node, runHmmObject));
    	    }
    	}

    	// ditto for all gene genealogies
    	for (HiddenState hiddenState : trees_states) {
    	    for (TNode node : hiddenState.getRootedGeneGenealogy().postTraverse()) {
    		// root has no incoming edge to optimize
    		if (node.isRoot()) {
    		    continue;
    		}

    		nobs.add(new GeneGenealogyNob(gBranch, branchMin, branchMax, node));
    	    }
    	}

    }




    public void runGridSearch(List<ObservationMap> observation, Hmm<ObservationMap> hmm,
            TransitionProbabilityParameters tpp,
            ArrayList<HiddenState> trees_states,
            Map<Network<CoalescePattern[]>,Set<HiddenState>> parentalTreeClasses) {


        initializeGridSearch(hmm, tpp, trees_states, parentalTreeClasses);

        double[] sampleInterval;
        double curMaxProb;
        double tempProb;

        for (Nob curNob: nobs) {
            curMaxProb = hmm.probability(observation);
            System.out.println("curMaxProb : " + curMaxProb);
            sampleInterval = curNob.getSamples();
            for (int i = 0; i < sampleInterval.length; i++) {
                curNob.set_param(sampleInterval[i]);
                tempProb = hmm.probability(observation);
                if (tempProb <= curMaxProb)
                    curNob.restoreParameterValue();
                else {
                    curMaxProb = tempProb;
                    System.out.println("tempMaxProb : " + tempProb);
                }
            }
        }

    }


}
