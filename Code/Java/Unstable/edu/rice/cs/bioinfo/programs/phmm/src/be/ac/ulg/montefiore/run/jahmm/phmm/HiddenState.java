/**
 * Container for data associated with a hidden state.
 * Neater to do it this way.
 * Modularity good - may want to extend it later.
 * 
 * Heaviest calculations go in here: 
 * - P[g|T] probability of gene genealogy given parental tree
 *   - parental tree -> gene genealogy -> double
 *     - Change one parental tree's branch length -> update its entry
 * - Conversion from GTR (rate matrix, base frequency vector, length of time) to 
 *   substitution probability matrix
 *   - substitution model object -> gene genealogy node -> double[][]
 *     - Change one gene genealogy's branch length -> update its entry.
 *     - Change any other substitution model parameter -> update entire cache.
 * - Emission probability calculation. Remember to use above substitution probability
 *   matrix cache when re-computing this one.
 *   - substitution model object -> gene genealogy -> ObservationMap -> double
 *     - Change one gene genealogy's branch length -> update its entry.
 *     - Change any other substitution model parameter -> update entire cache.
 *     - ObservationMap objects never change.
 *
 * All of the heaviest calculations must go through the cache.
 * Expensive to optimize substitution model parameters?
 * 
 * Perform caching for each parameter globally. Each HiddenState will have access to a global cache object.
 */

package be.ac.ulg.montefiore.run.jahmm.phmm;

import java.util.Map;
import java.util.Vector;
import java.util.List;
import java.util.ArrayList;
import java.util.Set;
import java.util.Iterator;
import java.io.StringReader;
import java.io.StringWriter;
import java.io.IOException;
//import phylogeny.EvoTree;
import phylogeny.Felsenstein;
import optimize.CalculationCache;
import substitutionModel.SubstitutionModel;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.ExNewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.RnNewickPrinter;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbabilityYF;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbabilityYF.CoalescePattern;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;
//import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbability;
import util.TreeUtils;

public class HiddenState
{
    public static final String HIDDEN_STATE_NAME_DELIMITER = ",";
    public static final String EQUIVALENCE_CLASS_NAME_DELIMITER = "!";

    // why not give it a unique name
    protected String name;
    protected Network<CoalescePattern[]> parentalTree;
    // For coalescent model calculations. 
    // Current model disregards gene genealogy branch lengths.
    protected Tree rootedGeneGenealogy; 
    // For substitution model calculations.
    // Shared among hidden states whose rootedGeneGenealogy
    // objects are equivalent after unrooting.
    //
    // Rooting of this object has no meaning.
    protected Tree unrootedGeneGenealogy; 

    // add a check to ensure that unrooted and rooted phylogeny
    // have Robinson-Foulds distance zero.
   
    // hmm, why not keep caches in HiddenState?
    // keep substitution probability cache in OpdfMap
    // keep transition probability cache in MultivariateOptimizer??

    // maintain mapping between taxa in parentalTree and geneGenealogy
    // only store reference to a shared object
    //
    // maps from species-name to list of allele-names
    // e.g. M -> AAAAAAAAQI,AAAAAAAAOG,AAAAAAAANQ
    // etc.
    protected Map<String,List<String>> alleleToSpeciesMapping;    

    // kliu - provide functionality to check if two HiddenState objects share a parental tree (Network<CoalescePattern[]> object)
    //protected Set<HiddenState> parentalTreeEquivalenceClass;

    // for emission probability calculation
    protected SubstitutionModel substitutionModel;

    // shared reference to a single global CalculationCache object
    protected CalculationCache calculationCache;

    /**
     * For emission probability calculation.
     * Use a standard phylogenetic substitution model.
     */
    protected Felsenstein felsensteinCalculator;

    /**
     * For coalescent model calculations.
     */
    protected GeneTreeProbabilityYF gtpyf;
    //protected GeneTreeProbability gtp;

    // cache it
    protected RnNewickPrinter<CoalescePattern[]> rnNewickPrinter;

    public HiddenState (String inName, 
			Network<CoalescePattern[]> inParentalTree, 
			Tree inRootedGeneGenealogy, 
			Tree inUnrootedGeneGenealogy, 
			Map<String,List<String>> inputSpeciesToAllelesMapping, 
			//Set<HiddenState> inParentalTreeEquivalenceClass, 
			SubstitutionModel inSubstitutionModel,
			CalculationCache inCalculationCache) {
	setName(inName);
	setParentalTree(inParentalTree);
	setRootedAndUnrootedGeneGenealogy(inRootedGeneGenealogy, inUnrootedGeneGenealogy);
	setSpeciesToAllelesMapping(inputSpeciesToAllelesMapping);
	//this.parentalTreeEquivalenceClass = inParentalTreeEquivalenceClass;
	setSubstitutionModel(inSubstitutionModel);
	this.calculationCache = inCalculationCache;
	// keep our own Felsenstein calculator
	// everything shares the same inCalculationCache anyways
	felsensteinCalculator = new Felsenstein(getSubstitutionModel(), calculationCache);
	gtpyf = new GeneTreeProbabilityYF();
	rnNewickPrinter = new RnNewickPrinter<CoalescePattern[]>();
    }
    
    public void setName (String inName) {
	this.name = inName;
    }

    public String getName () {
	return (name);
    }

    public Network<CoalescePattern[]> getParentalTree () {
    	return (parentalTree);
    }

    public Tree getRootedGeneGenealogy () {
    	return (rootedGeneGenealogy);
    }

    public Tree getUnrootedGeneGenealogy () {
    	return (unrootedGeneGenealogy);
    }

    public Map<String,List<String>> getSpeciesToAllelesMapping () {
    	return (alleleToSpeciesMapping);
    }

    public SubstitutionModel getSubstitutionModel () {
	return (substitutionModel);
    }

    public void setParentalTree (Network<CoalescePattern[]> inParentalTree) {
    	this.parentalTree = inParentalTree;
    }

    public void setRootedAndUnrootedGeneGenealogy (Tree inRootedGeneGenealogy, Tree inUnrootedGeneGenealogy) {
	this.rootedGeneGenealogy = inRootedGeneGenealogy;
    	this.unrootedGeneGenealogy = inUnrootedGeneGenealogy;
	// strict!
	if (!verifyRootedAndUnrootedGeneGenealogy()) {
	    throw (new RuntimeException ("ERROR: called HiddenState.setRootedAndUnrootedGeneGenealogy(...) with rooted phylogeny and unrooted phylogeny that had non-zero Robinson-Foulds distance."));
	}
    }

    /**
     * Make sure that Robinson-Foulds distance between rooted phylogeny and its unrooted
     * version is zero.
     */
    protected boolean verifyRootedAndUnrootedGeneGenealogy () {
	return (TreeUtils.calculateRobinsonFouldsDistance(rootedGeneGenealogy, unrootedGeneGenealogy) == 0);
    }

    public void setSpeciesToAllelesMapping (Map<String,List<String>> map) {
    	this.alleleToSpeciesMapping = map;
    }

    public void setSubstitutionModel (SubstitutionModel inSubstitutionModel) {
	this.substitutionModel = inSubstitutionModel;
    }

    /**
     * Checks to see if another HiddenState object hs shares a parental tree object
     * (hs.getParentalTree()) with this HiddenState object.
     */
    // public boolean checkSharedParentalTree (HiddenState hs) {
    // 	return (this.parentalTreeEquivalenceClass.contains(hs));
    // }

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
			"Rooted gene genealogy:\n" +
			// kliu - was toNewickWD - I don't think that was correct
			rootedGeneGenealogy.toNewick() + "\n" +
			"Unrooted gene genealogy:\n" +
			// kliu - was toNewickWD - I don't think that was correct
			unrootedGeneGenealogy.toNewick() + "\n");
    }
    

    /**
     * Expose this method to other classes, since transition probability calculation will need to access this.
     * Default to no debug messages.
     */
    public double calculateProbabilityOfRootedGeneGenealogyInParentalTree () {
	return (calculateProbabilityOfRootedGeneGenealogyInParentalTree(false));
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
    protected double calculateProbabilityOfRootedGeneGenealogyInParentalTree (boolean debugFlag) {
	// use cache if it exists
	if (calculationCache.cacheProbabilityOfGeneGenealogyInParentalTree.contains(getParentalTree(), getRootedGeneGenealogy())) {
	    return (calculationCache.cacheProbabilityOfGeneGenealogyInParentalTree.get(getParentalTree(), getRootedGeneGenealogy()).doubleValue());
	}

	// otherwise compute and cache
	double result = computeGTProb();

	// not quite right - while parental tree objects are unique by topology, gene genealogies aren't
	// under the current model
	// meh
	calculationCache.cacheProbabilityOfGeneGenealogyInParentalTree.put(getParentalTree(), getRootedGeneGenealogy(), new Double(result));
        return (result);
    }

    protected double computeGTProb () {
	return (computeGTProb(false));
    }

    /**
     * No real need to support non-binary trees or more than one tree.
     * Doesn't make sense in the PhyloNet-HMM model.
     * Eh, just leave it be.
     * Try to make as few changes to PhyloNet code as possible.
     */
    protected double computeGTProb (boolean debugFlag) {
	gtpyf.emptyState();

	// A list with one element.
	Vector<Tree> rootedGeneGenealogies = new Vector<Tree>();
	rootedGeneGenealogies.add(getRootedGeneGenealogy());

	// need to use these!
	List<Tree> bGeneTrees = new ArrayList<Tree>();
        List<List<Integer>> nbTree2bTrees = new ArrayList<List<Integer>>();
        for(Tree nbgt: rootedGeneGenealogies){
            List<Integer> bTrees = new ArrayList<Integer>();
            for(Tree bgt: Trees.getAllBinaryResolution(nbgt)){
                int index = 0;
                for(Tree exBgt: bGeneTrees){
                    if(Trees.haveSameRootedTopology(bgt,exBgt)){
                        break;
                    }
                    index++;
                }
		// bTrees is a list of all trees with identical rooted topology
		// *AND* multiplexed among all possible binary resolutions of nbgt
                if(index==bGeneTrees.size()){
                    bGeneTrees.add(bgt);
                }
                bTrees.add(index);
            }
	    // kliu
	    // nbTree2bTrees is a list of all possible resolutions of nbgt into binary trees
	    // bTrees is a list of indices
	    // add a whole list
            nbTree2bTrees.add(bTrees);
        }

        List<Double> probList;
	probList = gtpyf.calculateGTDistribution(parentalTree, bGeneTrees, getSpeciesToAllelesMapping(), 0);

	// look like the calculation is proceeding OK
	//
	// calculation under model from Yu et al. 2012
	// this method requires Network<CoalescePattern[]>
	// Yun uses Double to store hybridization probabilities during calculation
	// should only be a single entry
	if (probList.size() != bGeneTrees.size()) {
	    System.err.println ("ERROR: GeneTreeProbabilityYF.calculateGTDistribution(...) didn't return the same number of probabilities as the cardinality of the input set of gene trees. Returning -1 to signal error. " + probList.size() + " " + bGeneTrees.size());
	    return (-1.0);
	}

	    //}

	// only a single tree
	Vector<Double> geneTreeCountsVec = new Vector<Double>();
	geneTreeCountsVec.add(1.0);
	Iterator<Double> nbCounterIt = geneTreeCountsVec.iterator();
        //Iterator<Double> nbCounterIt = _geneTreeCounts.iterator();
        Iterator<List<Integer>> bGTIDs = nbTree2bTrees.iterator();
        double total = 0;
        for(Tree nbgt: rootedGeneGenealogies){
            // for(TNode node: nbgt.getNodes()){
	    // 	// kliu - hrm
	    // 	// this is no good
            //     node.setParentDistance(TNode.NO_DISTANCE);
            // }
            double maxProb = 0;
            for(int id: bGTIDs.next()){
                maxProb = Math.max(maxProb, probList.get(id));
            }
            double weight = nbCounterIt.next();
            total += Math.log(maxProb)*weight;
            if (debugFlag) {
		System.out.println("\n[x" + weight + "] " + nbgt.toString() + " : " + maxProb);
	    }
        }

	// bleh - make as few changes to computeGTProb(...) as possible
	return (Math.exp(total));
        //result.append("\n" + "Total log probability: " + total);
    }


    // /**
    //  * Perform standard coalescent model calculation to obtain
    //  * probability P[g(s_i) | T(s_i), c_{T(s_i)}] of observing a gene genealogy given a parental tree.
    //  *
    //  * See writeup for details.
    //  *
    //  * Might move this later to HiddenState.
    //  *
    //  * Use code from ComputeGTProb.
    //  * WARNING - returns likelihood, *NOT* log likelihood!
    //  *
    //  */
    // protected double calculateProbabilityOfGeneGenealogyInParentalTree (boolean debugFlag) {
    // 	// use cache if it exists
    // 	if (calculationCache.cacheProbabilityOfGeneGenealogyInParentalTree.contains(getParentalTree(), getGeneGenealogy())) {
    // 	    return (calculationCache.cacheProbabilityOfGeneGenealogyInParentalTree.get(getParentalTree(), getGeneGenealogy()).doubleValue());
    // 	}

    // 	// otherwise compute and cache

    // 	Map<String,String> alleleToSpeciesMapping = getAlleleToSpeciesMapping();
    // 	// A list with one element. Inefficient - consider doing a one-shot approach later.
    // 	Vector<Tree> geneGenealogies = new Vector<Tree>();
    // 	geneGenealogies.add(geneGenealogy);

    // 	gtp.emptyState();
	
    // 	// look like the calculation is proceeding OK
    // 	//
    // 	// calculation under model from Yu et al. 2012
    // 	// this method requires Network<Double>
    // 	// Yun uses Double to store hybridization probabilities during calculation
    //     List<Double> probList = gtp.calculateGTDistribution(parentalTree, geneGenealogies, alleleToSpeciesMapping, debugFlag);

    // 	// should only be a single entry
    // 	if (probList.size() != 1) {
    // 	    System.err.println ("ERROR: GeneTreeProbability.calculateGTDistribution(...) didn't return exactly one probability. Returning -1 to signal error.");
    // 	    return (-1.0);
    // 	}

    // 	double result = probList.get(0).doubleValue();
    // 	// not quite right - while parental tree objects are unique by topology, gene genealogies aren't
    // 	// under the current model
    // 	// meh
    // 	calculationCache.cacheProbabilityOfGeneGenealogyInParentalTree.put(getParentalTree(), getGeneGenealogy(), new Double(result));
    //     return (result);
    // }




    public double calculateEmissionProbability (ObservationMap o) {
	// if cache entry exists, use it
	if (calculationCache.cacheSubstitutionProbability.contains(getUnrootedGeneGenealogy(), o)) {
	    return (calculationCache.cacheSubstitutionProbability.get(getUnrootedGeneGenealogy(), o));
	}
	
	double result = felsensteinCalculator.getLikelihoodtree(getUnrootedGeneGenealogy(), o);
	calculationCache.cacheSubstitutionProbability.put(getUnrootedGeneGenealogy(), o, new Double(result));
	return (result);
    }
    
}
