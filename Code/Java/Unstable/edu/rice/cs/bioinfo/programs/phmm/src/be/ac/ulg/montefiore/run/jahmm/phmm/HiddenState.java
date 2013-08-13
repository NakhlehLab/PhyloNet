/**
 * Container for data associated with a hidden state.
 * Neater to do it this way.
 * Modularity good - may want to extend it later.
 */

package be.ac.ulg.montefiore.run.jahmm.phmm;

import java.util.Map;
import java.util.Vector;
import java.util.List;
import java.util.Set;
import java.io.StringReader;
import java.io.StringWriter;
import java.io.IOException;
//import phylogeny.EvoTree;

// kliu - Phylonet support libraries
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.ExNewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.RnNewickPrinter;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbability;


public class HiddenState
{
    public static final String HIDDEN_STATE_NAME_DELIMITER = ",";

    // why not give it a unique name
    protected String name;
    protected Network<Double> parentalTree;
    protected Tree geneGenealogy;

    // hmm, why not keep caches in hidden state?

    // maintain mapping between taxa in parentalTree and geneGenealogy
    // only store reference to a shared object
    protected Map<String,String> alleleToSpeciesMapping;    

    // kliu - provide functionality to check if two HiddenState objects share a parental tree (Network<Double> object)
    protected Set<HiddenState> parentalTreeEquivalenceClass;

    /**
     * For coalescent model calculations.
     */
    protected GeneTreeProbability gtp;

    // cache it
    protected RnNewickPrinter<Double> rnNewickPrinter;

    public HiddenState (String inName, Network<Double> inParentalTree, Tree inGeneGenealogy, Map<String,String> inputAlleleToSpeciesMapping, Set<HiddenState> inParentalTreeEquivalenceClass) {
	setName(inName);
	setParentalTree(inParentalTree);
	setGeneGenealogy(inGeneGenealogy);
	setAlleleToSpeciesMapping(inputAlleleToSpeciesMapping);
	this.parentalTreeEquivalenceClass = inParentalTreeEquivalenceClass;
	gtp = new GeneTreeProbability();
	rnNewickPrinter = new RnNewickPrinter<Double>();
    }
    
    // For test
    // public HiddenState (EvoTree inGeneGenealogy) {
    // 	setGeneGenealogy(inGeneGenealogy);
    // }

    public void setName (String inName) {
	this.name = inName;
    }

    public String getName () {
	return (name);
    }

    public Network<Double> getParentalTree () {
    	return (parentalTree);
    }

    public Tree getGeneGenealogy () {
    	return (geneGenealogy);
    }

    public Map<String,String> getAlleleToSpeciesMapping () {
    	return (alleleToSpeciesMapping);
    }

    public void setParentalTree (Network<Double> inParentalTree) {
    	this.parentalTree = inParentalTree;
    }

    public void setGeneGenealogy (Tree inGeneGenealogy) {
    	this.geneGenealogy = inGeneGenealogy;
    }

    public void setAlleleToSpeciesMapping (Map<String,String> map) {
    	this.alleleToSpeciesMapping = map;
    }

    /**
     * Checks to see if another HiddenState object hs shares a parental tree object
     * (hs.getParentalTree()) with this HiddenState object.
     */
    public boolean checkSharedParentalTree (HiddenState hs) {
	return (this.parentalTreeEquivalenceClass.contains(hs));
    }

    /**
     * Convenience function.
     */
    protected String getParentalTreeString () {
	// man, how do I empty the buffer/reset the writer?
	// don't think there's support for this
	StringWriter sw = new StringWriter();
	rnNewickPrinter.print(parentalTree, sw);
	String result = sw.toString();
	return (result);
    }

    /**
     * To String Method --> optional arguments
     * @param displayBranchLengths
     * @param displayInternalNodeNames
     * @return
     */
    public String toString () {
		return ("Parental tree:\n" +
			getParentalTreeString() + "\n" +
			"Gene genealogy:\n" +
			geneGenealogy.toNewickWD() + "\n");
    }
    

    /**
     * Expose this method to other classes, since transition probability calculation will need to access this.
     * Default to no debug messages.
     */
    public double calculateProbabilityOfGeneGenealogyInParentalTree () {
	return (calculateProbabilityOfGeneGenealogyInParentalTree(false));
    }

    /**
     * Perform standard coalescent model calculation to obtain
     * probability P[g(s_i) | T(s_i), c_{T(s_i)}] of observing a gene genealogy given a parental tree.
     *
     * See writeup for details.
     *
     * Might move this later to HiddenState.
     *
     * Use code from ComputeGTProb.
     * WARNING - returns likelihood, *NOT* log likelihood!
     *
     */
    protected double calculateProbabilityOfGeneGenealogyInParentalTree (boolean debugFlag) {
	Map<String,String> alleleToSpeciesMapping = getAlleleToSpeciesMapping();
	// A list with one element. Inefficient - consider doing a one-shot approach later.
	Vector<Tree> geneGenealogies = new Vector<Tree>();
	geneGenealogies.add(geneGenealogy);

	gtp.emptyState();
	
	// look like the calculation is proceeding OK
	//
	// calculation under model from Yu et al. 2012
	// this method requires Network<Double>
	// Yun uses Double to store hybridization probabilities during calculation
        List<Double> probList = gtp.calculateGTDistribution(parentalTree, geneGenealogies, alleleToSpeciesMapping, debugFlag);

	// should only be a single entry
	if (probList.size() != 1) {
	    System.err.println ("ERROR: GeneTreeProbability.calculateGTDistribution(...) didn't return exactly one probability. Returning -1 to signal error.");
	    return (-1.0);
	}

        return (probList.get(0));
    }


    
}
