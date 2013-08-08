package gridSearch;

import java.util.ArrayList;
import java.util.Map;
import java.util.Set;

//import phylogeny.EvoTree;
//import phylogeny.Node;
import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.Observation;
import be.ac.ulg.montefiore.run.jahmm.phmm.HiddenState;
import be.ac.ulg.montefiore.run.jahmm.phmm.OpdfMap;
import be.ac.ulg.montefiore.run.jahmm.phmm.TransitionProbabilityParameters;
import edu.rice.cs.bioinfo.library.programming.BijectiveHashtable;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;


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

    public GridSearchAlgorithm(int gBranchIn, int gRecombinationIn,
            int gHybridizationIn, int gBaseSubIn,
            double branchMinIn, double branchMaxIn,
            double recombinationMinIn, double recombinationMaxIn,
            double hybridizationMinIn, double hybridizationMaxIn,
			       double baseSubMinIn, double baseSubMaxIn
			       ) {

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


    }

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
     */
    private void initializeGridSearch(Hmm<O> hmm,
            TransitionProbabilityParameters tpp,
            ArrayList<HiddenState> trees_states,
            BijectiveHashtable<Network<Double>,Set<HiddenState>> parentalTreeClasses) {
        // Initialize Nobs array
        nobs = new ArrayList<Nob>();

        //add base sub nob
        nobs.add(new BaseSubNob(gBaseSub, branchMin, branchMax));

        //add transition probability recombination
        nobs.add(new RecombinationFreqNob(gRecombination, recombinationMin,
                recombinationMax, hmm, tpp, trees_states, parentalTreeClasses));

        //add transition probability hybridization
        nobs.add(new RecombinationFreqNob(gHybridization, hybridizationMin,
                hybridizationMax, hmm, tpp, trees_states, parentalTreeClasses));

	// kliu - need to use PhyloNet tree/network structures
	// walk through unique parental tree objects
	for (Network<Double> parentalTree : parentalTreeClasses.keys()) {
	    for (NetNode<Double> node : parentalTree.dfs()) {
		// root has no incoming edge to optimize
		if (node.isRoot()) {
		    continue;
		}

		// create a knob
		nobs.add(new ParentalTreeNob(gBranch, branchMin, branchMax, node));
	    }
	}

	// ditto for all gene genealogies
	for (HiddenState hiddenState : trees_states) {
	    for (TNode node : hiddenState.getGeneGenealogy().postTraverse()) {
		// root has no incoming edge to optimize
		if (node.isRoot()) {
		    continue;
		}

		nobs.add(new GeneGenealogyNob(gBranch, branchMin, branchMax, node));
	    }
	}

        //add tree branches
        // for (int i = 0; i < hmm.nbStates(); i++) {
        //     // get the current hidden state
        //     HiddenState currentState = ((OpdfMap)hmm.getOpdf(i)).getHiddenState();

        //     // Only if it's a parent that has not been encountered
        //     if (currentState.getParentalTreeID() != currentparentID) {
        //         //add all parent branches into the NOBs array
        //         getBranches(currentState.getParentalTree().getRoot());
        //     }

        //     // set the current Id
        //     currentparentID = currentState.getParentalTreeID();

        //     // add all gene tree branches into the NOBs array
        //     getBranches(currentState.getGeneGenealogy().getRoot());

        // }


    }


    /**
     * Helper function that aids in getting all the branches
     * in a tree.
     * Wraps each branch in a Nob Class and inserts this Nob object
     * into a the Nobs Array
     *
     * @param aNode
     */
    // private void getBranches(Node aNode) {
    //     if (!aNode.isRoot()) {
    //         //if Node is not a root then do this

    //         // add self to NOBs array
    //         nobs.add(new EvoTreeNob(gBranch, branchMin, branchMax, aNode));

    //         if (!aNode.isLeaf()) {
    //             // if node is not a leaf then run this method on its children
    //             ArrayList<Node> children = aNode.getChildren();
    //             for (int i = 0; i < children.size(); i++) {
    //                 getBranches(children.get(i));
    //             }
    //         }
    //     }
    // }


    public void runGridSearch(O observation, Hmm<O> hmm,
            TransitionProbabilityParameters tpp,
            ArrayList<HiddenState> trees_states,
            BijectiveHashtable<Network<Double>,Set<HiddenState>> parentalTreeClasses) {

        initializeGridSearch(hmm, tpp, trees_states, parentalTreeClasses);



    }


}