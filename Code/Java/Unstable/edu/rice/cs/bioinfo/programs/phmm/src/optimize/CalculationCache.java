/**
 * Heaviest calculations go in here: 
 * - P[g|T] probability of gene genealogy given parental tree
 *   - parental tree -> gene genealogy -> double
 *     - get/set: HiddenState.calculateProbabilityOfGeneGenealogyInParentalTree(...) DONE
 *     - clear:
 *       - Change one parental tree's branch length -> clear its entry DONE
 * - Conversion from GTR (rate matrix, base frequency vector, length of time) to 
 *   substitution probability matrix
 *   - substitution model object -> TNode -> double[][]
 *     - get/set: Felsenstein.calculatePij(...) DONE
 *     - clear:
 *       - Change one gene genealogy's branch length -> clear its entry. DONE
 *       - Change any other substitution model parameter -> clear entire cache. DONE
 * - Emission probability calculation. Remember to use above substitution probability
 *   matrix cache when re-computing this one.
 *   - substitution model object -> gene genealogy -> ObservationMap -> double
 *     - get/set: HiddenState.calculateEmissionProbability(...) DONE
 *     - clear:
 *       - Change one gene genealogy's branch length -> clear its entry. DONE
 *       - Change any other substitution model parameter -> clear entire cache. DONE
 *       - ObservationMap objects never change.
 *
 * All of the heaviest calculations must go through the cache.
 * Expensive to optimize substitution model parameters?
 * 
 * One global CalculationCache object will need to be shared amongst
 * the HiddenState objects.
 *
 * Perform caching for each parameter globally. Each HiddenState will have access to a global cache object.
 *
 * No value guards in here.
 */

package optimize;

import java.util.Hashtable;
import util.MapOfMap;
import be.ac.ulg.montefiore.run.jahmm.phmm.ObservationMap;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

public class CalculationCache {
    // for convenience, just provide full external access
    public MapOfMap<Network<Double>,Tree,Double> cacheProbabilityOfGeneGenealogyInParentalTree;
    public Hashtable<TNode,double[][]> cacheSubstitutionProbabilityMatrix;
    public MapOfMap<Tree,ObservationMap,Double> cacheSubstitutionProbability;

    public CalculationCache () {
    	cacheProbabilityOfGeneGenealogyInParentalTree = new MapOfMap<Network<Double>,Tree,Double>();
    	cacheSubstitutionProbabilityMatrix = new Hashtable<TNode,double[][]>();
    	cacheSubstitutionProbability = new MapOfMap<Tree,ObservationMap,Double>();
    }

}