package edu.rice.cs.bioinfo.programs.phylonet.algos.search.HillClimbing;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkNeighbourhoodGenerator;

import java.util.Comparator;
/*
 *@ClassName: SimulatedAnnealingTA
 *@Description: This is the search for network inference based on a tree
 *@Author: Zhen Cao
 *@Date:  2019-06-11 17:04
 *@Version: 1.0
 */
public class SimpleHillClimbingTA extends SimpleHillClimbing {
    public SimpleHillClimbingTA(Comparator<Double> scoreComparator, NetworkNeighbourhoodGenerator generator) {
        super(scoreComparator, generator);

    }

    protected void search(Network currentNetwork, Func1<Network,Double> getScore, double initialScore)
    {
        _networkGenerator.resetList();
        super.search(currentNetwork, getScore, initialScore);
    }

}
