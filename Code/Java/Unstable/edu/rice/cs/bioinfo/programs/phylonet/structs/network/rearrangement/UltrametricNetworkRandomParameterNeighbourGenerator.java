package edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.UnorderedPair;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: yy9
 * Date: 6/7/12
 * Time: 9:26 PM
 * To change this template use File | Settings | File Templates.
 */
public class UltrametricNetworkRandomParameterNeighbourGenerator extends NetworkRandomParameterNeighbourGenerator {
    public UltrametricNetworkRandomParameterNeighbourGenerator(Map<UnorderedPair, Double> pairwiseTimeLimit, Func1<Integer, Integer> numParametersToChange, Long seed) {
        super(numParametersToChange, seed);
        _lengthChanger = new NodeHeightChange(pairwiseTimeLimit);
    }

    public UltrametricNetworkRandomParameterNeighbourGenerator(Map<UnorderedPair, Double> pairwiseTimeLimit,Func1<Integer, Integer> numParametersToChange) {
        this(pairwiseTimeLimit, numParametersToChange, null);
    }


    public UltrametricNetworkRandomParameterNeighbourGenerator(Map<UnorderedPair, Double> pairwiseTimeLimit) {
        this(pairwiseTimeLimit, new Func1<Integer, Integer>() {
            @Override
            public Integer execute(Integer input) {
                return 2;
            }
        });
    }



    protected void getAllParameters(Network network, ArrayList<Tuple<Tuple<NetNode, NetNode>, Boolean>> allEdges) {
        for (Object nodeO : Networks.postTraversal(network)) {
            NetNode node = (NetNode) nodeO;
            if(!node.isLeaf()){
                Tuple<NetNode, NetNode> edge = new Tuple<>(null, node);
                allEdges.add(new Tuple<Tuple<NetNode, NetNode>, Boolean>(edge, true));
            }
            if (node.isNetworkNode()) {
                Tuple<NetNode, NetNode> edge = new Tuple<>(null, node);
                allEdges.add(new Tuple<Tuple<NetNode, NetNode>, Boolean>(edge, false));
            }
        }
    }




}