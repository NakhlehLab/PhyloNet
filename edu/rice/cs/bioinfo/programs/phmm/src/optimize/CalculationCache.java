/**
 * This file is part of PhyloNet-HMM.
 *
 * Copyright © 2013-2014 Kevin Liu, Jingxuan Dai, Kathy Truong, 
 * Ying Song, Michael H. Kohn, and Luay Nakhleh. <http://bioinfo.cs.rice.edu/>
 * 
 * PhyloNet-HMM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * PhyloNet-HMM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * Heaviest calculations go in here: 
 * - P[g|T] probability of gene genealogy given parental tree
 *   - parental tree -> gene genealogy -> double
 *     - get/set: HiddenState.calculateProbabilityOfRootedGeneGenealogyInParentalTree(...) DONE
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
import be.ac.ulg.montefiore.run.jahmm.phmm.HiddenState;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbabilityYF.CoalescePattern;

public class CalculationCache {
    // for convenience, just provide full external access
    public MapOfMap<Network<CoalescePattern[]>,Tree,Double> cacheProbabilityOfGeneGenealogyInParentalTree;
    public Hashtable<TNode,double[][]> cacheSubstitutionProbabilityMatrix;
    public MapOfMap<Tree,ObservationMap,Double> cacheSubstitutionProbability;
    // check non-null and non-empty to check cache ready state
    //public MapOfMap<HiddenState,HiddenState,Double> cacheSwitchingFrequencyMap;

    public CalculationCache () {
    	cacheProbabilityOfGeneGenealogyInParentalTree = new MapOfMap<Network<CoalescePattern[]>,Tree,Double>();
    	cacheSubstitutionProbabilityMatrix = new Hashtable<TNode,double[][]>();
    	cacheSubstitutionProbability = new MapOfMap<Tree,ObservationMap,Double>();
	// special case, since we use null-ness (as well as 
	// empty-ness) to signal cached quantities not ready for
	// switching frequencies
	//cacheSwitchingFrequencyMap = null;
    }

}
