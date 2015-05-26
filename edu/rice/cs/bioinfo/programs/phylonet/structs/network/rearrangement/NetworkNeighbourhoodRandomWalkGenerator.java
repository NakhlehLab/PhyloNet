package edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement;

import edu.rice.cs.bioinfo.library.programming.*;
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
public class NetworkNeighbourhoodRandomWalkGenerator extends NetworkNeighbourhoodGenerator{
    private double[] _operationProbabilities;
    private int _operationID;
    private Tuple<NetNode,NetNode> _targetEdge;
    private Tuple<NetNode,NetNode> _sourceEdge;
    private Tuple<NetNode,NetNode> _destinationEdge;
    private int[][] _nodeDistanceMatrix;
    private int _moveDiameterLimit;
    private int _reticulationDiameterLimit;
    private Map<NetNode, Integer> _node2ID;
    private Random _random;
    private NetworkRearrangementOperation[] _networkOperators;
    private boolean _printDetails = false;
    private int _maxReticulations;
    private Set<String> _singleAlleleSpecies;


    public NetworkNeighbourhoodRandomWalkGenerator(double[] probabilities, int maxReticulations, int moveDiameterLimit, int reticulationDiameterLimit, double brlenChangeWindow, double inheriProbChangeWindow, Set<String> singleAlleleSpecies, Long seed)
    {
        _networkOperators = new NetworkRearrangementOperation[7];
        _networkOperators[0] = new ReticulationEdgeAddition();
        _networkOperators[1] = new ReticulationEdgeDeletion();
        _networkOperators[2] = new ReticulationEdgeDestinationChange();
        _networkOperators[3] = new EdgeSourceChange();
        _networkOperators[4] = new ReticulationFlip();
        EdgeLengthChange edgeLengthChange = new EdgeLengthChange();
        edgeLengthChange.setWindowSize(brlenChangeWindow);
        _networkOperators[5] = edgeLengthChange;
        EdgeInheritanceProbabilityChange edgeInheritanceProbabilityChange = new EdgeInheritanceProbabilityChange();
        edgeInheritanceProbabilityChange.setWindowSize(inheriProbChangeWindow);
        _networkOperators[6] = edgeInheritanceProbabilityChange;
        _reticulationDiameterLimit = reticulationDiameterLimit;
        _moveDiameterLimit = moveDiameterLimit;
        _operationProbabilities = new double[_networkOperators.length];
        for(int i=0; i<_operationProbabilities.length; i++){
            if(i==0){
                _operationProbabilities[i] = probabilities[i];
            }
            else{
                _operationProbabilities[i] = _operationProbabilities[i-1] + probabilities[i];
            }
        }
        if(_operationProbabilities[_operationProbabilities.length-1]!=1){
            for(int i=0; i<_operationProbabilities.length; i++){
                _operationProbabilities[i] = _operationProbabilities[i]/_operationProbabilities[_operationProbabilities.length-1];
            }
        }
        if(seed!=null) {
            _random = new Random(seed);
        }else{
            _random = new Random();
        }
        _maxReticulations = maxReticulations;
        _singleAlleleSpecies = singleAlleleSpecies;
    }

    public NetworkNeighbourhoodRandomWalkGenerator(double[] probabilities, int maxReticulations, int moveDiameterLimit, int reticulationDiameterLimit, Set<String> singleAlleleSpecies, Long seed){
        this(probabilities, maxReticulations, moveDiameterLimit, reticulationDiameterLimit, 0.1, 0.1, singleAlleleSpecies, seed);
    }

    public NetworkNeighbourhoodRandomWalkGenerator(double[] probabilities, int maxReticulations, int moveDiameterLimit, int reticulationDiameterLimit, Long seed){
        this(probabilities, maxReticulations, moveDiameterLimit, reticulationDiameterLimit, new HashSet<String>(), seed);
    }

    public NetworkNeighbourhoodRandomWalkGenerator(double[] probabilities, int maxReticulations, Long seed){
        this(probabilities, maxReticulations, -1, -1, seed);
    }

    public NetworkNeighbourhoodRandomWalkGenerator(double[] probabilities, int maxReticulations){
        this(probabilities, maxReticulations, null);
    }


    public void mutateNetwork(Network network){
        ArrayList<Tuple<NetNode, NetNode>> allEdges = new ArrayList<>();
        ArrayList<Tuple<NetNode, NetNode>> allEdgesNeedBrlens = new ArrayList<>();
        ArrayList<Tuple<NetNode, NetNode>> allReticulationEdges = new ArrayList<>();
        ArrayList<Tuple<NetNode, NetNode>> allRemovableReticulationEdges = new ArrayList<>();
        Set<String> taxa = new HashSet<>();
        Ref<Boolean> incrementHybrid = new Ref<>(false);
        getNetworkInfo(network, allEdges, allEdgesNeedBrlens, allReticulationEdges, allRemovableReticulationEdges, taxa, incrementHybrid);

        if(_reticulationDiameterLimit!=-1 || _moveDiameterLimit!=-1){
            _node2ID = new HashMap<NetNode, Integer>();
            int id = 0;
            for(Object node: Networks.postTraversal(network)){
                _node2ID.put((NetNode)node, id++);
            }
            computeNodeDistances(network);
        }
        //System.out.println("#reticulationEdges:"+allReticulationEdges.size());
        if(_printDetails){
            System.out.println("Before rearrangement: "+ network.toString());
        }
        boolean successRearrangment;
        boolean[] triedOperations = new boolean[_operationProbabilities.length];
        do{
            do {
                _operationID = getNextOperationID(incrementHybrid.get(), allRemovableReticulationEdges.size()!=0);
            }while(triedOperations[_operationID]);
            Set<Integer> previousTriedEdges = new HashSet<>();
            successRearrangment = false;
            if(_printDetails){
                System.out.println("Select operation: "+ _networkOperators[_operationID].getClass().getName());
            }
            while(!successRearrangment && setNextMove(network, allEdges, allEdgesNeedBrlens, allReticulationEdges, allRemovableReticulationEdges, previousTriedEdges)){
                if (_networkOperators[_operationID].performOperation()) {

                    if (Networks.hasCycle(network) || !isNetworkValid(network, taxa)) {
                        if(_printDetails){
                            System.out.println(": Invalid");
                        }
                        successRearrangment = false;
                        _networkOperators[_operationID].undoOperation();
                    }
                    else{
                        successRearrangment = true;
                        if(_printDetails){
                            System.out.println();
                        }
                    }
                }
            }
            if(!successRearrangment){
                triedOperations[_operationID] = true;
            }
        }while(!successRearrangment);

        if(_printDetails){
            System.out.println("After rearrangement: " + network.toString());
        }

        //rearrangementComputed.execute(network,_operationID, new Tuple3(_targetEdge,_sourceEdge,_destinationEdge));
        //_networkOperators[_operationID].undoOperation();

    }

    public void undo(){
        _networkOperators[_operationID].undoOperation();
    }


    private void computeNodeDistances(Network network){
        _nodeDistanceMatrix = new int[_node2ID.size()][_node2ID.size()];
        HashMap<NetNode, List<Tuple<NetNode,Integer>>> node2children = new HashMap<NetNode, List<Tuple<NetNode,Integer>>>();
        for(Object nodeO: Networks.postTraversal(network)){
            NetNode node = (NetNode)nodeO;
            List<Tuple<NetNode, Integer>> child2depth = new ArrayList<Tuple<NetNode, Integer>>();
            child2depth.add(new Tuple<NetNode,Integer>(node, 0));
            int nodeID = _node2ID.get(node);
            List<NetNode> children = new ArrayList<NetNode>();
            for(Object child: node.getChildren()){
                children.add((NetNode)child);
                for(Tuple<NetNode, Integer> tuple: node2children.get(child)){
                    int depth = tuple.Item2+1;
                    child2depth.add(new Tuple<NetNode, Integer>(tuple.Item1, depth));
                    int childID = _node2ID.get(tuple.Item1);
                    int dist = _nodeDistanceMatrix[nodeID][childID];
                    if(dist==0 || dist>depth){
                        _nodeDistanceMatrix[nodeID][childID] = depth - 1;
                        _nodeDistanceMatrix[childID][nodeID] = depth - 1;
                    }
                }
            }
            if(children.size()==2){
                for(Tuple<NetNode, Integer> tuple1: node2children.get(children.get(0))){
                    int child1ID = _node2ID.get(tuple1.Item1);
                    int depth1 = tuple1.Item2;
                    for(Tuple<NetNode, Integer> tuple2: node2children.get(children.get(1))){
                        int child2ID = _node2ID.get(tuple2.Item1);
                        int depth2 = tuple2.Item2;

                        if(child1ID == child2ID){
                            continue;
                        }
                        int dist = _nodeDistanceMatrix[child1ID][child2ID];
                        int newDist = depth1 + depth2;
                        if(dist==0 || dist>newDist){
                            _nodeDistanceMatrix[child1ID][child2ID] = newDist;
                            _nodeDistanceMatrix[child2ID][child1ID] = newDist;
                        }
                    }
                }
            }
            node2children.put(node, child2depth);
        }
    }


    private int getNextOperationID(boolean incrementHybrid, boolean decrementHybrid){

        boolean stop;
        int operationID = -1;
        do{
            double random = _random.nextDouble();
            stop = true;
            for(int i=0; i<_operationProbabilities.length; i++){
                if(random<_operationProbabilities[i]){
                    operationID = i;
                    break;
                }
            }

            if((operationID == 0 && !incrementHybrid) || ((operationID == 1||operationID == 2) && !decrementHybrid)){
                stop = false;
            }

        }while(!stop);
        return operationID;
    }




    private void getNetworkInfo(Network network, ArrayList<Tuple<NetNode, NetNode>> allEdges, ArrayList<Tuple<NetNode, NetNode>> edgesNeedBrlens, ArrayList<Tuple<NetNode, NetNode>> allReticulationEdges, ArrayList<Tuple<NetNode, NetNode>> removableReticulationEdges, Set<String> leafSet, Ref<Boolean> increaseReticulations){
        int numReticulations = 0;
        Map<NetNode,Set<String>> node2leaves = new HashMap<>();
        for(Object nodeO: Networks.postTraversal(network)){
            NetNode node = (NetNode)nodeO;
            Set<String> leaves = new HashSet<>();
            if(node.isLeaf()){
                leafSet.add(node.getName());
                leaves.add(node.getName());
            }
            else if(node.isNetworkNode()){
                numReticulations++;
            }
            for(Object childO: node.getChildren()){
                NetNode childNode = (NetNode)childO;
                Set<String> childLeaves = node2leaves.get(childNode);
                leaves.addAll(childLeaves);
                Tuple<NetNode,NetNode> edge = new Tuple<>(node, childNode);
                allEdges.add(edge);
                if(childLeaves.size()!=1 || !_singleAlleleSpecies.containsAll(childLeaves)){
                    edgesNeedBrlens.add(edge);
                }
                if(childNode.isNetworkNode()) {
                    if (node.isTreeNode()) {
                        removableReticulationEdges.add(edge);
                    }
                    allReticulationEdges.add(edge);
                }

            }
            node2leaves.put(node,leaves);
        }
        if(numReticulations>=_maxReticulations){
            increaseReticulations.set(false);
        }
        else{
            increaseReticulations.set(true);
        }
    }

    private boolean setNextMove(Network network, ArrayList<Tuple<NetNode,NetNode>> allEdges, ArrayList<Tuple<NetNode,NetNode>> allEdgesNeedBrlens, ArrayList<Tuple<NetNode,NetNode>> allReticulationEdges,ArrayList<Tuple<NetNode,NetNode>> allRemovableReticulationEdges, Set<Integer> previousTriedEdges){
        boolean successMove = true;
        switch(_operationID){
            case 0:
                setParametersForReticulationEdgeAddition(allEdges, previousTriedEdges);
                break;
            case 1:
                successMove = setParametersForReticulationEdgeDeletion(allRemovableReticulationEdges, previousTriedEdges);
                break;
            case 2:
                successMove = setParametersForReticulationEdgeDestinationChange(allRemovableReticulationEdges,allEdges, previousTriedEdges);
                break;
            case 3:
                successMove = setParametersForEdgeSourceChange(allEdges, previousTriedEdges);
                break;
            case 4:
                successMove = setParametersForReticulationFlip(allRemovableReticulationEdges, previousTriedEdges);
                break;
            case 5:
                successMove = setParametersForEdgeLengthChange(allEdgesNeedBrlens);
                break;
            case 6:
                successMove = setParametersForEdgeInheritanceProbabilityChange(allReticulationEdges);
                break;
        }
        if(successMove){
            _networkOperators[_operationID].setParameters(network, _targetEdge, _sourceEdge, _destinationEdge);
        }
        return successMove;
    }


    private void setParametersForReticulationEdgeAddition(ArrayList<Tuple<NetNode,NetNode>> allEdges, Set<Integer> edgesTried){
        _targetEdge = null;
        int size = allEdges.size();
        boolean endSampling;
        do{
            endSampling = true;
            int sourceID = _random.nextInt(size);
            int destinationID = sourceID;
            while(sourceID==destinationID){
                destinationID = _random.nextInt(size);
            }
            int tupleID = (int)(Math.pow(10,new String(allEdges.size()+"").length())*sourceID) + destinationID;
            if(edgesTried.contains(tupleID)){
                endSampling = false;
            }
            else{
                _sourceEdge = allEdges.get(sourceID);
                _destinationEdge = allEdges.get(destinationID);
                edgesTried.add(tupleID);
                if(_reticulationDiameterLimit!=-1){
                    if(_nodeDistanceMatrix[_node2ID.get(_sourceEdge.Item2)][_node2ID.get(_destinationEdge.Item2)]>_reticulationDiameterLimit){
                        endSampling = false;
                    }
                }
            }
        }while(!endSampling);
        if(_printDetails){
            System.out.print("Add " + printEdge(_sourceEdge) + " to " + printEdge(_destinationEdge));
        }

    }

    private boolean setParametersForReticulationEdgeDeletion(ArrayList<Tuple<NetNode,NetNode>> allRemovableReticulationEdges, Set<Integer> edgesTried){
        _sourceEdge = null;
        _destinationEdge = null;
        int size = allRemovableReticulationEdges.size();
        int randomID;
        do{
            if(edgesTried.size()==size){
                return false;
            }
            randomID = _random.nextInt(size);
        }while(edgesTried.contains(randomID));
        _targetEdge = allRemovableReticulationEdges.get(randomID);
        edgesTried.add(randomID);
        if(_printDetails){
            System.out.print("Remove " + printEdge(_targetEdge));
        }
        return true;
    }


    private boolean setParametersForReticulationEdgeDestinationChange(ArrayList<Tuple<NetNode,NetNode>> allRemovableReticulationEdges, ArrayList<Tuple<NetNode,NetNode>> allEdges, Set<Integer> edgesTried){
        _sourceEdge = null;
        int reticulationEdgeSize = allRemovableReticulationEdges.size();
        int edgeSize = allEdges.size();

        boolean endSampling;
        do{
            if(edgesTried.size()==reticulationEdgeSize*edgeSize){
                return false;
            }

            endSampling = true;
            int targetID = _random.nextInt(reticulationEdgeSize);
            _targetEdge = allRemovableReticulationEdges.get(targetID);

            int destinationID = _random.nextInt(edgeSize);
            _destinationEdge = allEdges.get(destinationID);

            int tupleID = (int)(Math.pow(10,new String(allEdges.size()+"").length())*targetID) + destinationID;
            if(edgesTried.contains(tupleID)){
                endSampling = false;
                continue;
            }

            edgesTried.add(tupleID);
            if(_reticulationDiameterLimit!=-1){
                if(_nodeDistanceMatrix[_node2ID.get(_targetEdge.Item1)][_node2ID.get(_destinationEdge.Item2)]>_reticulationDiameterLimit){
                    endSampling = false;
                    continue;
                }
            }
            if(_moveDiameterLimit!=-1){
                if(_nodeDistanceMatrix[_node2ID.get(_targetEdge.Item2)][_node2ID.get(_destinationEdge.Item2)]>_moveDiameterLimit){
                    endSampling = false;
                    continue;
                }
            }

            if(_destinationEdge.Item2.equals(_targetEdge.Item1) || _destinationEdge.Item1.equals(_targetEdge.Item1)
                    || _destinationEdge.Item2.equals(_targetEdge.Item2) || _destinationEdge.Item1.equals(_targetEdge.Item2)){
                //edgesTried.add(edgeTuple);
                endSampling = false;
                continue;
            }

        }while(!endSampling);

        if(_printDetails){
            System.out.print("Redirect destination of " + printEdge(_targetEdge) + " to " + printEdge(_destinationEdge));
        }
        return true;
    }

    private boolean setParametersForEdgeSourceChange(ArrayList<Tuple<NetNode,NetNode>> allEdges, Set<Integer> edgesTried){
        _sourceEdge = null;
        int allEdgeSize = allEdges.size();

        boolean endSampling;
        do{
            if(edgesTried.size()==allEdgeSize*(allEdgeSize-1)){
                return false;
            }
            endSampling = true;
            int targetID = _random.nextInt(allEdgeSize);
            //int targetID = 5;
            _targetEdge = allEdges.get(targetID);
            //System.out.print("Target: " + printEdge(_targetEdge));

            if(_targetEdge.Item1.isNetworkNode()){
                for(int i=0; i<allEdgeSize; i++){
                    if(i!=targetID)
                        edgesTried.add((int)(Math.pow(10,new String(allEdgeSize+"").length()))*targetID + i);
                }
                endSampling = false;
                continue;
            }

            int destinationID = targetID;
            //int destinationID = 2;
            while(targetID==destinationID){
                destinationID = _random.nextInt(allEdgeSize);
            }
            //System.out.println(" Destination: " + printEdge(_destinationEdge));

            int tupleID = (int)(Math.pow(10,new String(allEdgeSize+"").length()))*targetID + destinationID;
            if(edgesTried.contains(tupleID)){
                endSampling = false;
                continue;
            }
            edgesTried.add(tupleID);
            _destinationEdge = allEdges.get(destinationID);

            if(_moveDiameterLimit!=-1){
                if(_nodeDistanceMatrix[_node2ID.get(_targetEdge.Item2)][_node2ID.get(_destinationEdge.Item2)]>_moveDiameterLimit){
                    endSampling = false;
                    continue;
                }
            }

            if(_reticulationDiameterLimit!=-1 && _targetEdge.Item2.isNetworkNode()){
                if(_nodeDistanceMatrix[_node2ID.get(_targetEdge.Item2)][_node2ID.get(_destinationEdge.Item2)]>_reticulationDiameterLimit){
                    endSampling = false;
                    continue;
                }
            }

            if(_targetEdge.Item1.isRoot()){
                boolean valid = true;
                for(Object anotherChild: _targetEdge.Item1.getChildren()){
                    if(anotherChild!=_targetEdge.Item2 && ((NetNode)anotherChild).isNetworkNode()){
                        valid = false;
                    }
                }
                if(!valid){
                    for(int i=0; i<allEdgeSize; i++){
                        if(i!=targetID)
                            edgesTried.add((int)(Math.pow(10,new String(allEdgeSize+"").length()))*targetID + i);
                    }
                    endSampling = false;
                    continue;
                }
            }

            if(_targetEdge.Item1.equals(_destinationEdge.Item1) || _targetEdge.Item1.equals(_destinationEdge.Item2) ||
                    _targetEdge.Item2.equals(_destinationEdge.Item1) || _targetEdge.Item2.equals(_destinationEdge.Item2)){
                endSampling = false;
                continue;
            }

        }while(!endSampling);
        if(_printDetails){
            System.out.print("Redirect source of " + printEdge(_targetEdge) + " to " + printEdge(_destinationEdge));
        }
        return true;
    }


    private boolean setParametersForReticulationFlip(ArrayList<Tuple<NetNode,NetNode>> allRemovableReticulationEdges, Set<Integer> edgesTried){
        _sourceEdge = null;
        _destinationEdge = null;
        boolean endSampling;
        int size = allRemovableReticulationEdges.size();
        do{
            if(edgesTried.size()==size){
                return false;
            }
            endSampling = true;
            int targetID;
            do{
                targetID = _random.nextInt(size);
            }while(edgesTried.contains(targetID));
            edgesTried.add(targetID);
            _targetEdge = allRemovableReticulationEdges.get(targetID);
            if(_targetEdge.Item1.isRoot()){
                endSampling = false;
            }

        }while(!endSampling);
        if(_printDetails){
            System.out.print("Flip reticulation edge " + printEdge(_targetEdge));
        }
        return true;
    }

    private boolean setParametersForEdgeLengthChange(ArrayList<Tuple<NetNode,NetNode>> allEdges){
        _sourceEdge = null;
        _destinationEdge = null;
        _targetEdge = allEdges.get(_random.nextInt(allEdges.size()));
        if(_printDetails){
            System.out.print("Change length of " + printEdge(_targetEdge));
        }
        return true;
    }


    private boolean setParametersForEdgeInheritanceProbabilityChange(ArrayList<Tuple<NetNode,NetNode>> allReticulationEdges){
        if(allReticulationEdges.size()==0)return false;
        _sourceEdge = null;
        _destinationEdge = null;
        _targetEdge = allReticulationEdges.get(_random.nextInt(allReticulationEdges.size()));
        if(_printDetails){
            System.out.print("Change inheritance probability of " + printEdge(_targetEdge));
        }
        return true;
    }




    public void performRearrangement(Network network, Integer operation, Tuple<NetNode,NetNode> targetEdge, Tuple<NetNode,NetNode> sourceEdge, Tuple<NetNode,NetNode> destinationEdge){
        _networkOperators[operation].setParameters(network, targetEdge, sourceEdge, destinationEdge);
        _networkOperators[operation].performOperation();
    }




}