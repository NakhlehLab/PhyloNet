package edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.UnorderedPair;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;

/**
 * Created by Yun Yu
 *
 * This class is a subclass of NetworkRandomParameterNeighbourGenerator.
 * It generates a random neighbor of a given network by changing its parameters
 * Note that in this class both the original network and the resulting network need to be ultrametric
 */
public class UltrametricNetworkRandomParameterNeighbourGenerator extends NetworkRandomParameterNeighbourGenerator {

    /**
     * Constructor of this class
     *
     * @param numParametersToChange     the function that calculates the number of parameters to change in one move
     * @param pairwiseTimeLimit         the temporal constraints
     * @param seed                      the seed for controlling randomness
     */
    public UltrametricNetworkRandomParameterNeighbourGenerator(Map<UnorderedPair, Double> pairwiseTimeLimit, Func1<Integer, Integer> numParametersToChange, Long seed) {
        super(numParametersToChange, seed);
        _lengthChanger = new NodeHeightChange(pairwiseTimeLimit);
    }


    /**
     * Constructor of this class
     *
     * @param numParametersToChange     the function that calculates the number of parameters to change in one move
     * @param pairwiseTimeLimit         the temporal constraints
     */
    public UltrametricNetworkRandomParameterNeighbourGenerator(Map<UnorderedPair, Double> pairwiseTimeLimit,Func1<Integer, Integer> numParametersToChange) {
        this(pairwiseTimeLimit, numParametersToChange, null);
    }


    /**
     * Constructor of this class
     * Note that here the default is set to change 2 parameters at a time
     *
     * @param pairwiseTimeLimit         the temporal constraints
     */
    public UltrametricNetworkRandomParameterNeighbourGenerator(Map<UnorderedPair, Double> pairwiseTimeLimit) {
        this(pairwiseTimeLimit, new Func1<Integer, Integer>() {
            @Override
            public Integer execute(Integer input) {
                return 2;
            }
        });
    }


    /**
     * This function is to get all parameters of the network
     *
     * @param network   the species network
     * @param allEdges  all edges in the species network
     */
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