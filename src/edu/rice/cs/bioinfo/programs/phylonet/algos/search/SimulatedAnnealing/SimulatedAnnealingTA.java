package edu.rice.cs.bioinfo.programs.phylonet.algos.search.SimulatedAnnealing;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Ref;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement.NetworkNeighbourhoodGenerator;
/*
 *@ClassName: SimulatedAnnealingTA
 *@Description: This is the search for network inference based on a tree
 *@Author: Zhen Cao
 *@Date:  2019-06-11 17:04
 *@Version: 1.0
 */
import java.util.Comparator;

public class SimulatedAnnealingTA extends SimulatedAnnealingSalterPearL {
    private int [] _acceptCount;


    public SimulatedAnnealingTA(Comparator<Double> scoreComparator, NetworkNeighbourhoodGenerator generator, Long seed) {
        super(scoreComparator, generator, seed);
        _acceptCount = new int[]{0, 0, 0, 0, 0, 0, 0};
    }

    protected void search(Network currentNetwork, Func1<Network,Double> getScore, double initialScore)
    {

        _networkGenerator.resetList();
        super.search(currentNetwork, getScore, initialScore);
    }

    public void handleAcceptCase(Network proposedNetwork, Ref<Double> currentScore, double newScore){
        super.handleAcceptCase(proposedNetwork, currentScore, newScore);
//        System.out.println("id:"+_networkGenerator.getOperationID());
        int id = _networkGenerator.getOperationID();
        if(id != -1){
            _acceptCount[_networkGenerator.getOperationID()]++;

        }
    }

    public int[] GetAcceptCount() {
        return _acceptCount;
    }
}
