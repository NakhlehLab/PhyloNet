package edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement;

import edu.rice.cs.bioinfo.library.programming.Func;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Ref;
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
public class NetworkRandomParameterNeighbourGenerator extends NetworkNeighbourhoodGenerator {
    private EdgeLengthChange _lengthChanger = new EdgeLengthChange();
    private EdgeInheritanceProbabilityChange _inheriProbChanger = new EdgeInheritanceProbabilityChange();
    private Func1<Integer, Integer> _numParametersToChange;
    private Set<String> _singleAlleleSpecies;
    private List<Tuple<Tuple<Tuple<NetNode, NetNode>, Boolean>, Double>> _edgeChanged = new ArrayList<>();
    private Random _random;


    public NetworkRandomParameterNeighbourGenerator(Func1<Integer, Integer> numParametersToChange, Set<String> singleAlleleSpecies, Long seed) {
        _numParametersToChange = numParametersToChange;
        _singleAlleleSpecies = singleAlleleSpecies;
        if (seed != null) {
            _random = new Random(seed);
        } else {
            _random = new Random();
        }
    }

    public NetworkRandomParameterNeighbourGenerator(Func1<Integer, Integer> numParametersToChange, Set<String> singleAlleleSpecies) {
        this(numParametersToChange, singleAlleleSpecies, null);
    }


    public NetworkRandomParameterNeighbourGenerator(Set<String> singleAlleleSpecies) {
        this(new Func1<Integer, Integer>() {
            @Override
            public Integer execute(Integer input) {
                return 2;
            }
        }, singleAlleleSpecies, null);
    }

    public void setParametersToChange(Func1<Integer, Integer> numParametersToChange){
        _numParametersToChange = numParametersToChange;
    }

    public NetworkRandomParameterNeighbourGenerator() {
        this(new HashSet<String>());

    }

    public void setWindowSize(double window){
        _lengthChanger.setWindowSize(window);
        _inheriProbChanger.setWindowSize(window);
    }


    public void mutateNetwork(Network network) {
        if(_printDetails){
            System.out.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
            System.out.println("Before rearrangement: "+ network.toString());
            System.out.println("Select operation: "+ this.getClass().getSimpleName());
        }
        ArrayList<Tuple<Tuple<NetNode, NetNode>, Boolean>> allEdges = new ArrayList<>();
        getAllParameters(network, allEdges);
        _edgeChanged.clear();
        int numEdgeToChange = _numParametersToChange.execute(allEdges.size());

        for (int i = 0; i < numEdgeToChange; i++) {
            Tuple<Tuple<NetNode, NetNode>, Boolean> selectEdgeInfo = allEdges.remove(_random.nextInt(allEdges.size()));
            Tuple<NetNode, NetNode> selectEdge = selectEdgeInfo.Item1;
            if (selectEdgeInfo.Item2) {
                _edgeChanged.add(new Tuple<Tuple<Tuple<NetNode, NetNode>, Boolean>, Double>(selectEdgeInfo, selectEdge.Item2.getParentDistance(selectEdge.Item1)));
                _lengthChanger.setParameters(null, selectEdge, -1, -1);
                _lengthChanger.performOperation();
                if (_printDetails) {
                    System.out.println("Change length of: " + printEdge(selectEdgeInfo.Item1));
                }
            } else {
                _edgeChanged.add(new Tuple<Tuple<Tuple<NetNode, NetNode>, Boolean>, Double>(selectEdgeInfo, selectEdge.Item2.getParentProbability(selectEdge.Item1)));
                _inheriProbChanger.setParameters(null, selectEdge, -1, -1);
                _inheriProbChanger.performOperation();
                if (_printDetails) {
                    System.out.println("Change probability of: " + printEdge(selectEdgeInfo.Item1));
                }
            }
        }
        if (_printDetails) {
            System.out.println("After rearrangement: " + network.toString());
            System.out.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
        }
        //rearrangementComputed.execute(network,_operationID, new Tuple3(_targetEdge,_sourceEdge,_destinationEdge));
        //_networkOperators[_operationID].undoOperation();

    }

    public void undo() {
        for(Tuple<Tuple<Tuple<NetNode, NetNode>, Boolean>, Double> selectEdgeInfo: _edgeChanged){
            if(selectEdgeInfo.Item1.Item2){
                _lengthChanger.setParameters(null, selectEdgeInfo.Item1.Item1, selectEdgeInfo.Item2, -1);
                _lengthChanger.performOperation();
            }
            else{
                _inheriProbChanger.setParameters(null, selectEdgeInfo.Item1.Item1, -1, selectEdgeInfo.Item2);
                _inheriProbChanger.performOperation();
            }
        }
    }


    private void getAllParameters(Network network, ArrayList<Tuple<Tuple<NetNode, NetNode>, Boolean>> allEdges) {
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