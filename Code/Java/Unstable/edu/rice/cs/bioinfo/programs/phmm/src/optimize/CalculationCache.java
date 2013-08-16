/**
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
 * One global CalculationCache object will need to be shared amongst
 * the HiddenState objects.
 *
 * Perform caching for each parameter globally. Each HiddenState will have access to a global cache object.
 */

package optimize;

public class CalculationCache {
    public CalculationCache () {
	// TODO
    }
}
