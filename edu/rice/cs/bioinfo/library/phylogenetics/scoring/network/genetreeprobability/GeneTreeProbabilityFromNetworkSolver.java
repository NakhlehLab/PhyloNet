package edu.rice.cs.bioinfo.library.phylogenetics.scoring.network.genetreeprobability;

import edu.rice.cs.bioinfo.library.phylogenetics.Graph;

import java.util.Map;
import java.util.Set;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/26/12
 * Time: 3:16 PM
 * To change this template use File | Settings | File Templates.
 */
public interface GeneTreeProbabilityFromNetworkSolver<NN,NE,TN,TE,P>
{
    public Map<Graph<TN,TE>,P> computeGeneTreeProbabilities(Graph<NN,NE> network, Set<Graph<TN,TE>> geneTrees);
}
