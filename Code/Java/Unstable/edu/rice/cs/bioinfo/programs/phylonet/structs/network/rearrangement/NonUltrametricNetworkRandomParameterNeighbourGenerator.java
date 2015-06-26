package edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Tuple;
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
public class NonUltrametricNetworkRandomParameterNeighbourGenerator extends NetworkRandomParameterNeighbourGenerator {
    private Set<String> _singleAlleleSpecies;

    public NonUltrametricNetworkRandomParameterNeighbourGenerator(Func1<Integer, Integer> numParametersToChange, Set<String> singleAlleleSpecies, Long seed) {
        super(numParametersToChange, seed);
        _singleAlleleSpecies = singleAlleleSpecies;
        _lengthChanger = new EdgeLengthChange();
    }

    public NonUltrametricNetworkRandomParameterNeighbourGenerator(Func1<Integer, Integer> numParametersToChange, Set<String> singleAlleleSpecies) {
        this(numParametersToChange, singleAlleleSpecies, null);
    }


    public NonUltrametricNetworkRandomParameterNeighbourGenerator(Set<String> singleAlleleSpecies) {
        this(new Func1<Integer, Integer>() {
            @Override
            public Integer execute(Integer input) {
                return 2;
            }
        }, singleAlleleSpecies);
    }


    public NonUltrametricNetworkRandomParameterNeighbourGenerator() {
        this(new HashSet<String>());

    }

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




}