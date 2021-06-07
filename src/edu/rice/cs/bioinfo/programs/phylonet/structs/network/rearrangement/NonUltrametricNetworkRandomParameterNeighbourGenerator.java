package edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;

/**
 * Created by Yun Yu
 *
 * This class is a subclass of NetworkRandomParameterNeighbourGenerator.
 * It generates a random neighbor of a given network by changing its parameters
 * Note that in this class neither the original network nor the resulting network need to be ultrametric
 */
public class NonUltrametricNetworkRandomParameterNeighbourGenerator extends NetworkRandomParameterNeighbourGenerator {
    private Set<String> _singleAlleleSpecies;


    /**
     * Constructor of this class
     *
     * @param numParametersToChange     the function that calculates the number of parameters to change in one move
     * @param singleAlleleSpecies       all edges that has only one leaf under it who is in this set do not need to change their lengths
     * @param seed                      the seed for controlling randomness
     */
    public NonUltrametricNetworkRandomParameterNeighbourGenerator(Func1<Integer, Integer> numParametersToChange, Set<String> singleAlleleSpecies, Long seed) {
        super(numParametersToChange, seed);
        _singleAlleleSpecies = singleAlleleSpecies;
        _lengthChanger = new EdgeLengthChange();
    }


    /**
     * Constructor of this class
     *
     * @param numParametersToChange     the function that calculates the number of parameters to change in one move
     * @param singleAlleleSpecies       all edges that has only one leaf under it who is in this set do not need to change their lengths
     */
    public NonUltrametricNetworkRandomParameterNeighbourGenerator(Func1<Integer, Integer> numParametersToChange, Set<String> singleAlleleSpecies) {
        this(numParametersToChange, singleAlleleSpecies, null);
    }



    /**
     * Constructor of this class
     * Note that here the default is set to change 2 parameters at a time
     *
     * @param singleAlleleSpecies       all edges that has only one leaf under it who is in this set do not need to change their lengths
     */
    public NonUltrametricNetworkRandomParameterNeighbourGenerator(Set<String> singleAlleleSpecies) {
        this(new Func1<Integer, Integer>() {
            @Override
            public Integer execute(Integer input) {
                return 2;
            }
        }, singleAlleleSpecies);
    }


    /**
     * Constructor of this class
     * Note that here the default is set to change 2 parameters at a time
     */
    public NonUltrametricNetworkRandomParameterNeighbourGenerator() {
        this(new HashSet<String>());

    }


    /**
     * This function is to get all parameters of the network
     *
     * @param network   the species network
     * @param allEdges  all edges in the species network
     */
    protected void getAllParameters(Network network, ArrayList<Tuple<Tuple<NetNode, NetNode>, Boolean>> allEdges) {
        Map<NetNode, Set<String>> node2leaves = new HashMap<>();
        for (Object nodeO : Networks.postTraversal(network)) {
            NetNode node = (NetNode) nodeO;
            Set<String> leaves = new HashSet<>();
            if(node.isLeaf()){
                leaves.add(node.getName());
            }
            for (Object childO : node.getChildren()) {
                NetNode childNode = (NetNode) childO;
                Set<String> childLeaves = node2leaves.get(childNode);
                leaves.addAll(childLeaves);
                Tuple<NetNode, NetNode> edge = new Tuple<>(node, childNode);
                if (childLeaves.size() != 1 || !_singleAlleleSpecies.containsAll(childLeaves)) {
                    allEdges.add(new Tuple<Tuple<NetNode, NetNode>, Boolean>(edge, true));
                }
            }
            if (node.isNetworkNode()) {
                Tuple<NetNode, NetNode> edge = new Tuple<>((NetNode) node.getParents().iterator().next(), node);
                allEdges.add(new Tuple<Tuple<NetNode, NetNode>, Boolean>(edge, false));
            }
            node2leaves.put(node, leaves);
        }

    }


    public void resetList(){}
    public int getOperationID(){
        return -1;
    }


}