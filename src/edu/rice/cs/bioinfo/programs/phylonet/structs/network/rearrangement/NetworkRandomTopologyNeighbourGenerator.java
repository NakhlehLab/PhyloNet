package edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement;

import edu.rice.cs.bioinfo.library.programming.Ref;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;

/**
 * Created by Yun Yu
 *
 * This class is a subclass of NetworkNeighbourhoodGenerator.
 * It generates a random neighbor of a given network by changing its topology
 */
public class NetworkRandomTopologyNeighbourGenerator extends NetworkNeighbourhoodGenerator{
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
    private int _maxReticulations;
    private Set<String> _fixedHybrid;



    /**
     * Constructor of this class
     *
     * @param probabilities         the weights of different rearrangement moves
     * @param maxReticulations      the maximum number of reticulations allowed
     * @param moveDiameterLimit     the maximum diameter for random move
     * @param reticulationDiameterLimit           the maximum diameter for a reticulation (the distance between two parents)
     * @param fixedHybrid           the user-specified set of hybrid species
     * @param seed                  the seed for controlling randomness
     */
    public NetworkRandomTopologyNeighbourGenerator(double[] probabilities, int maxReticulations, int moveDiameterLimit, int reticulationDiameterLimit, Set<String> fixedHybrid, Long seed)
    {
        _networkOperators = new NetworkRearrangementOperation[6];
        _networkOperators[0] = new ReticulationEdgeAddition();
        _networkOperators[1] = new ReticulationEdgeDeletion();
        _networkOperators[2] = new ReticulationEdgeDestinationChange();
        _networkOperators[3] = new EdgeSourceChange();
        _networkOperators[4] = new ReticulationFlip();
        _networkOperators[5] = new ReticulationEdgeReplace();
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
        _fixedHybrid = fixedHybrid;
        if(_fixedHybrid!=null && _fixedHybrid.isEmpty()){
            _fixedHybrid = null;
        }
    }


    /**
     * This function is to mutate the current network by changing its topology
     */
    public void mutateNetwork(Network network){
        ArrayList<Tuple<NetNode, NetNode>> allEdges = new ArrayList<>();
        ArrayList<Tuple<NetNode, NetNode>> allRemovableReticulationEdges = new ArrayList<>();
        Ref<Boolean> incrementHybrid = new Ref<>(false);
        Set<String> taxa = new HashSet<>();
        getNetworkInfo(network, taxa, allEdges, allRemovableReticulationEdges, incrementHybrid);

        if(_reticulationDiameterLimit!=-1 || _moveDiameterLimit!=-1){
            _node2ID = new HashMap<NetNode, Integer>();
            int id = 0;
            for(Object node: Networks.postTraversal(network)){
                _node2ID.put((NetNode)node, id++);
            }
            computeNodeDistances(network);
        }
        if(_printDetails){
            System.out.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
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
                System.out.println("Select operation: "+ _networkOperators[_operationID].getClass().getSimpleName());
            }
            while(!successRearrangment && setNextMove(network, allEdges, allRemovableReticulationEdges, previousTriedEdges)){
                if (_networkOperators[_operationID].performOperation()) {

                    if (Networks.hasCycle(network) || !isNetworkValid(network, taxa) || !checkHybrid(network, _fixedHybrid)) {
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
            System.out.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
        }
    }


    /**
     * This function is to undo the last rearrangement move
     */
    public void undo(){
        _networkOperators[_operationID].undoOperation();
    }



    /**
     * This function is to compute the distance between every pair of nodes in the network
     */
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


    /**
     * This function is to generate the next random rearrangement
     *
     * @param incrementHybrid   whether adding reticulations is allowed
     * @param decrementHybrid   whether deleting reticulations is allowed
     */
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


    /**
     * This function is to gather the information about the network
     *
     * @param network       the given species network
     * @param allEdges      all edges in the species network
     * @param removableReticulationEdges      all reticulation edges that can be removed
     *                                        note that a reticulation edge can be removed only if removing it reduces the number of reticulations in the network by one
     * @param increaseReticulations     whether adding reticulations is allowed
     */
    private void getNetworkInfo(Network network, Set<String> taxa, ArrayList<Tuple<NetNode, NetNode>> allEdges, ArrayList<Tuple<NetNode, NetNode>> removableReticulationEdges, Ref<Boolean> increaseReticulations){
        int numReticulations = 0;
        for(Object nodeO: Networks.postTraversal(network)){
            NetNode node = (NetNode)nodeO;
            if(node.isLeaf()){
                taxa.add(node.getName());
            }
            else if(node.isNetworkNode()){
                numReticulations++;
            }
            for(Object childO: node.getChildren()){
                NetNode childNode = (NetNode)childO;
                Tuple<NetNode,NetNode> edge = new Tuple<>(node, childNode);
                allEdges.add(edge);
                if(childNode.isNetworkNode()) {
                    if (node.isTreeNode()) {
                        removableReticulationEdges.add(edge);
                    }
                }

            }
        }
        if(numReticulations>=_maxReticulations){
            increaseReticulations.set(false);
        }
        else{
            increaseReticulations.set(true);
        }
    }


    /**
     * This function is to set the next move
     *
     * @param network       the given species network
     * @param allEdges      all edges in the species network
     * @param allRemovableReticulationEdges      all reticulation edges that can be removed
     *                                        note that a reticulation edge can be removed only if removing it reduces the number of reticulations in the network by one
     * @param previousTriedEdges     all edges that have been tried
     */
    private boolean setNextMove(Network network, ArrayList<Tuple<NetNode,NetNode>> allEdges, ArrayList<Tuple<NetNode,NetNode>> allRemovableReticulationEdges, Set<Integer> previousTriedEdges){
        boolean successMove = true;
        switch(_operationID){
            case 0:
                successMove = setParametersForReticulationEdgeAddition(allEdges, previousTriedEdges);
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
                successMove = setParametersForReticulationEdgeReplace(allRemovableReticulationEdges, allEdges, previousTriedEdges);
                break;
        }
        if(successMove){
            _networkOperators[_operationID].setParameters(network, _targetEdge, _sourceEdge, _destinationEdge);
        }
        return successMove;
    }


    /**
     * This function is to set parameters for adding reticulation edge
     *
     * @param allEdges      all edges in the species network
     * @param edgesTried    all edges that have been tried which results in invalid networks
     */
    private boolean setParametersForReticulationEdgeAddition(ArrayList<Tuple<NetNode,NetNode>> allEdges, Set<Integer> edgesTried){
        _targetEdge = null;
        int size = allEdges.size();
        boolean endSampling;
        do{
            if(edgesTried.size()==size*(size-1)/2){
                return false;
            }
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
        return true;
    }



    /**
     * This function is to set parameters for deleting reticulation edge
     *
     * @param allRemovableReticulationEdges   all reticulation edges in the species network that can be removed
     *                                        note that a reticulation edge can be removed only if removing it reduces the number of reticulations in the network by one
     * @param edgesTried    all edges that have been tried which results in invalid networks
     */
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


    /**
     * This function is to set parameters for changing the tail of a reticulation edge
     *
     * @param allRemovableReticulationEdges   all reticulation edges in the species network that can be edited
     *                                        note that a reticulation edge can be edited only if editing it without deleting does not change the number of reticulations in the network
     * @param allEdges      all edges in the species network
     * @param edgesTried    all edges that have been tried which results in invalid networks
     */
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
                endSampling = false;
                continue;
            }

        }while(!endSampling);

        if(_printDetails){
            System.out.print("Redirect destination of " + printEdge(_targetEdge) + " to " + printEdge(_destinationEdge));
        }
        return true;
    }


    /**
     * This function is to set parameters for changing the head of a reticulation edge
     *
     * @param allEdges      all edges in the species network
     * @param edgesTried    all edges that have been tried which results in invalid networks
     */
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
            _targetEdge = allEdges.get(targetID);

            if(_targetEdge.Item1.isNetworkNode()){
                for(int i=0; i<allEdgeSize; i++){
                    if(i!=targetID)
                        edgesTried.add((int)(Math.pow(10,new String(allEdgeSize+"").length()))*targetID + i);
                }
                endSampling = false;
                continue;
            }

            int destinationID = targetID;
            while(targetID==destinationID){
                destinationID = _random.nextInt(allEdgeSize);
            }

            int tupleID = (int)(Math.pow(10,new String(allEdgeSize+"").length()))*targetID + destinationID;
            if(edgesTried.contains(tupleID)){
                endSampling = false;
                continue;
            }
            edgesTried.add(tupleID);
            _destinationEdge = allEdges.get(destinationID);

            if(_moveDiameterLimit!=-1){
                if(_nodeDistanceMatrix[_node2ID.get(_targetEdge.Item1)][_node2ID.get(_destinationEdge.Item2)]>_moveDiameterLimit){
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



    /**
     * This function is to set parameters for flipping the direction of a reticulation edge
     *
     * @param allRemovableReticulationEdges   all reticulation edges in the species network that can be edited
     *                                        note that a reticulation edge can be edited only if editing it without deleting does not change the number of reticulations in the network
     * @param edgesTried    all edges that have been tried which results in invalid networks
     */
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


    /**
     * This function is to set parameters for replacing a reticulation edge
     *
     * @param allRemovableReticulationEdges   all reticulation edges in the species network that can be replaced
     *                                        note that a reticulation edge can be edited only if replacing it does not change the number of reticulations in the network
     * @param edgesTried    all edges that have been tried which results in invalid networks
     */
    private boolean setParametersForReticulationEdgeReplace(ArrayList<Tuple<NetNode,NetNode>> allRemovableReticulationEdges, ArrayList<Tuple<NetNode,NetNode>> allEdges, Set<Integer> edgesTried){
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

        size = allEdges.size();
        boolean endSampling;
        do{
            endSampling = true;
            int sourceID = _random.nextInt(size);
            _sourceEdge = allEdges.get(sourceID);
            if(_sourceEdge.Item1 == _targetEdge.Item1 || _sourceEdge.Item2 == _targetEdge.Item1){
                endSampling = false;
            }
            else{
                int destinationID = sourceID;
                while(sourceID==destinationID){
                    destinationID = _random.nextInt(size);
                }
                _destinationEdge = allEdges.get(destinationID);
                if(_destinationEdge.Item1 == _targetEdge.Item2 || _destinationEdge.Item2 == _targetEdge.Item2){
                    endSampling = false;
                }
                else{
                    if(_sourceEdge.Item1 == _targetEdge.Item2 || _sourceEdge.Item2 == _targetEdge.Item2){
                        if(_destinationEdge.Item1 == _targetEdge.Item1 || _destinationEdge.Item2 == _targetEdge.Item1) {
                            endSampling = false;

                        }
                    }
                    else{
                        if(_reticulationDiameterLimit!=-1){
                            if(_nodeDistanceMatrix[_node2ID.get(_sourceEdge.Item2)][_node2ID.get(_destinationEdge.Item2)]>_reticulationDiameterLimit){
                                endSampling = false;
                            }
                        }

                        if(_moveDiameterLimit!=-1){
                            if(_nodeDistanceMatrix[_node2ID.get(_targetEdge.Item1)][_node2ID.get(_sourceEdge.Item2)]>_moveDiameterLimit ||
                                    _nodeDistanceMatrix[_node2ID.get(_targetEdge.Item2)][_node2ID.get(_destinationEdge.Item2)]>_moveDiameterLimit){
                                endSampling = false;
                            }
                        }
                    }
                }
            }

        }while(!endSampling);
        if(_printDetails){
            System.out.print(" and add " + printEdge(_sourceEdge) + " to " + printEdge(_destinationEdge) );
        }
        return true;
    }


    /**
     * This function is to perform the rearrangement move
     *
     * @param network       the network to be mutated
     * @param operation     the rearrangement operation that to be performed
     * @param targetEdge    the target edge that is involved in this move
     * @param sourceEdge    the source edge that is involved in this move
     * @param destinationEdge    the destination edge that is involved in this move
     */
    public void performRearrangement(Network network, Integer operation, Tuple<NetNode,NetNode> targetEdge, Tuple<NetNode,NetNode> sourceEdge, Tuple<NetNode,NetNode> destinationEdge){
        _networkOperators[operation].setParameters(network, targetEdge, sourceEdge, destinationEdge);
        _networkOperators[operation].performOperation();
    }


    public void resetList(){}
    public int getOperationID(){
        return -1;
    }

}