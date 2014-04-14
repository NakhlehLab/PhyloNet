package edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.allNeighbours;

import edu.rice.cs.bioinfo.library.phylogenetics.GetDirectSuccessors;
import edu.rice.cs.bioinfo.library.phylogenetics.GetInDegree;
import edu.rice.cs.bioinfo.library.phylogenetics.GetNodesPostOrder;
import edu.rice.cs.bioinfo.library.phylogenetics.Graph;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func3;
import edu.rice.cs.bioinfo.library.programming.Func4;
import edu.rice.cs.bioinfo.library.programming.Tuple;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: yy9
 * Date: 6/7/12
 * Time: 9:26 PM
 * To change this template use File | Settings | File Templates.
 */
public class NetworkNeighbourhoodRandomWalkGenerator<G extends Graph<N,E>,N,E> extends NetworkNeighbourhoodGenerator<G,N,E> {
    private final double[] _operationProbabilities;
    private int _operationID;
    private ArrayList<Set<Integer>> _previousTried;
    private E _operationEdge1;
    private E _operationEdge2;
    private int[][] _nodeDistanceMatrix;
    private int _diameterLimit;
    private Map<N, Integer> _node2ID;
    private Random _random;


    public NetworkNeighbourhoodRandomWalkGenerator(double[] probabilities, Func1<G, N> makeNode, Func3<G, N, N, E> makeEdge, Long seed)
    {
        super(makeNode,makeEdge);
        _operationProbabilities = new double[4];
        for(int i=0; i<4; i++){
            if(i==0){
                _operationProbabilities[i] = probabilities[i];
            }
            else{
                _operationProbabilities[i] = _operationProbabilities[i-1] + probabilities[i];
            }
        }
        _previousTried = new ArrayList<Set<Integer>>();
        for(int i=0; i<4; i++){
            _previousTried.add(new HashSet<Integer>());
        }
        //initializeReticulationNodeSet();
        if(seed!=null) {
            _random = new Random(seed);
        }else{
            _random = new Random();
        }
    }

    /*
    public NetworkNeighbourhoodRandomWalkGenerator(double[] probabilities, Func1<G, N> makeNode, Func3<G, N, N, E> makeEdge)
    {
        super(makeNode,makeEdge);
        _operationProbabilities = new double[4];
        for(int i=0; i<4; i++){
            if(i==0){
                _operationProbabilities[i] = probabilities[i];
            }
            else{
                _operationProbabilities[i] = _operationProbabilities[i-1] + probabilities[i];
            }
        }
        _previousTried = new ArrayList<Set<Integer>>();
        for(int i=0; i<4; i++){
            _previousTried.add(new HashSet<Integer>());
        }
        //initializeReticulationNodeSet();
        _random = new Random();
    }
    */



    public NetworkNeighbourhoodRandomWalkGenerator(Func1<G, N> makeNode, Func3<G, N, N, E> makeEdge)
    {
        super(makeNode,makeEdge);
        _operationProbabilities = null;
        //initializeReticulationNodeSet();
        _previousTried = new ArrayList<Set<Integer>>();
        for(int i=0; i<4; i++){
            _previousTried.add(new HashSet<Integer>());
        }
    }

    /*
    private void initializeReticulationNodeSet(){
        _reticulationNodeSet = new HashSet<N>();
        GetInDegree<N,E> getInDegree = new GetInDegree<N,E>();
        for(N node: _network.getNodes()){
            if(getInDegree.execute(_network,node)==2){
                _reticulationNodeSet.add(node);
            }
        }
    }
    */

    public G computeRandomNeighbour(G network, boolean incrementHybrid, Func4<G,Integer,E,E,Boolean> rearrangementComputed, int diameterLimit){
        _diameterLimit = diameterLimit;
        boolean findValidNetwork;
        ArrayList<E> allEdges = new ArrayList<E>();
        ArrayList<E> allReticulationEdges = new ArrayList<E>();
        ArrayList<E> allRemovableReticulationEdges = new ArrayList<E>();
        ArrayList<E> allTreeEdges = new ArrayList<E>();
        computeNetworkEdges(network, allEdges, allTreeEdges, allReticulationEdges, allRemovableReticulationEdges);
        _operationID = getNextOperationID(incrementHybrid, allRemovableReticulationEdges.size()!=0);

        if(diameterLimit!=0){
            _node2ID = new HashMap<N, Integer>();
            int id = 0;
            for(N node: network.getNodes()){
                _node2ID.put(node, id++);
            }
            computeNodeDistances(network);
        }
        //System.out.println("#reticulationEdges:"+allReticulationEdges.size());
        boolean operationValid;
        do{
            operationValid = true;
            findValidNetwork = true;
            try{
                operationValid = setNextMove(network, allEdges, allTreeEdges, allReticulationEdges,allRemovableReticulationEdges);
                if(operationValid){
                    _networkOperators[_operationID].performOperation();
                    assertValidNetwork(network);
                }
            }catch(IllegalArgumentException e){
                _networkOperators[_operationID].undoOperation();
                findValidNetwork = false;
            }catch(IllegalStateException e){
                findValidNetwork = false;
            }

        }while(!findValidNetwork && operationValid);

        /*
        N reticulationNode = null;
        if(_operationID == 0){
            reticulationNode = (N)((ReticulationEdgeAddition)_networkOperators[_operationID]).getReticulationNode();
            _reticulationNodeSet.add(reticulationNode);
        }else if(_operationID == 1){
            reticulationNode = (N)((ReticulationEdgeDeletion)_networkOperators[_operationID]).getReticulationNode();
            _reticulationNodeSet.remove(reticulationNode);
        }
        */
        _previousTried.get(_operationID).clear();

        if(operationValid){
            rearrangementComputed.execute(network,_operationID,_operationEdge1,_operationEdge2);
            _networkOperators[_operationID].undoOperation();
        }else{
            network = computeRandomNeighbour(network, incrementHybrid, rearrangementComputed);
        }

        return network;
    }

    public G next(G network, int operationID){
        _operationID = operationID;
        boolean findValidNetwork;
        ArrayList<E> allEdges = new ArrayList<E>();
        ArrayList<E> allTreeEdges = new ArrayList<E>();
        ArrayList<E> allReticulationEdges = new ArrayList<E>();
        ArrayList<E> allRemovableReticulationEdges = new ArrayList<E>();
        computeNetworkEdges(network, allEdges, allTreeEdges, allReticulationEdges,allRemovableReticulationEdges);
        do{
            try{
                setNextMove(network, allEdges,allTreeEdges, allReticulationEdges,allRemovableReticulationEdges);
                _networkOperators[_operationID].performOperation();
                assertValidNetwork(network);
                //assertUnallowedCase();
                findValidNetwork = true;
            }catch(IllegalArgumentException e){
                System.out.println("Not valid network!");
                _networkOperators[_operationID].undoOperation();
                findValidNetwork = false;
            }catch(IllegalStateException e){
                System.out.println("Not valid state!");
                findValidNetwork = false;
            }
        }while(!findValidNetwork);

        _previousTried.get(_operationID).clear();
        /*
        N reticulationNode;
        if(_operationID == 0){
            reticulationNode = (N)((ReticulationEdgeAddition)_networkOperators[_operationID]).getReticulationNode();
            _reticulationNodeSet.add(reticulationNode);
        }else if(_operationID == 1){
            reticulationNode = (N)((ReticulationEdgeDeletion)_networkOperators[_operationID]).getReticulationNode();
            _reticulationNodeSet.remove(reticulationNode);
        }
        */
        return network;
    }

    /*
    private void assertUnallowedCase(){
        N reticulationNode = null;
        try{

            if(_operationID == 0){
                reticulationNode = (N)((ReticulationEdgeAddition)_networkOperators[_operationID]).getReticulationNode();
                _reticulationNodeSet.add(reticulationNode);
            }else if(_operationID == 1){
                reticulationNode = (N)((ReticulationEdgeDeletion)_networkOperators[_operationID]).getReticulationNode();
                _reticulationNodeSet.remove(reticulationNode);
            }
            GetDirectPredecessors<N,E> getParent = new GetDirectPredecessors<N,E>();
            for(N node1: _reticulationNodeSet){
                N node1Parent1=null, node1Parent2=null;
                int index = 0;
                for(N n1p: getParent.execute(_network,node1)){
                    if(index==0){
                        node1Parent1 = n1p;
                    }
                    else{
                        node1Parent2 = n1p;
                    }
                    index++;
                }
                

                for(N node2: _reticulationNodeSet){
                    if(node1.equals(node2))continue;
                    N node2Parent1=null, node2Parent2=null;
                    index = 0;
                    for(N n2p: getParent.execute(_network,node2)){
                        if(index==0){
                            node2Parent1 = n2p;
                        }
                        else{
                            node2Parent2 = n2p;
                        }
                        index++;
                    }
                    GetInDegree<N,E> getInDegree = new GetInDegree<N, E>();
                    if(_network.getEdge(node1,node2)!=null){
                        N anotherParent=node2Parent1.equals(node1)?node2Parent2:node2Parent1;
                        if((anotherParent.equals(node1Parent1) || anotherParent.equals(node1Parent2)) && (_network.getEdge(node1Parent1, node1Parent2)!=null || _network.getEdge(node1Parent2, node1Parent1)!=null)){
                            System.out.println(1);
                            throw new IllegalArgumentException();
                        }

                        
                        if((_network.getEdge(node1Parent1, anotherParent) != null || _network.getEdge(node1Parent2, anotherParent) != null) && getInDegree.execute(_network,anotherParent)==1){
                            System.out.println(2);
                            throw new IllegalArgumentException();
                        }

                        else if(_network.getEdge(anotherParent, node1Parent1) != null && getInDegree.execute(_network,node1Parent1)==1){
                            System.out.println(3);
                            throw new IllegalArgumentException();
                        }
                        else if(_network.getEdge(anotherParent, node1Parent2) != null && getInDegree.execute(_network,node1Parent2)==1){
                            System.out.println(4);
                            throw new IllegalArgumentException();
                        }
                    }
                    else{

                        if((node1Parent1.equals(node2Parent1)||node1Parent2.equals(node2Parent1)) && _network.getEdge(node1,node2Parent2)!=null && getInDegree.execute(_network,node2Parent2)==1){
                            System.out.println(5);
                            throw new IllegalArgumentException();
                        }
                        if((node1Parent1.equals(node2Parent2)||node1Parent2.equals(node2Parent2)) && _network.getEdge(node1,node2Parent1)!=null && getInDegree.execute(_network,node2Parent1)==1){
                            System.out.println(6);
                            throw new IllegalArgumentException();
                        }
                        
                    }
                }
            }
        }catch(IllegalArgumentException e){
            System.out.println("Caused by " + _operationID);
            throw new IllegalArgumentException();
        }finally {
            if(_operationID == 0){
                _reticulationNodeSet.remove(reticulationNode);
            }else if(_operationID == 1){
                _reticulationNodeSet.add(reticulationNode);
            }
        }
    }
    */



    public G undo(G network){
        _networkOperators[_operationID].undoOperation();
        return network;
    }


    private void computeNodeDistances(G network){
        _nodeDistanceMatrix = new int[_node2ID.size()][_node2ID.size()];
        GetNodesPostOrder<N, E> postOrderTraversal = new GetNodesPostOrder<N, E>();
        HashMap<N, List<Tuple<N,Integer>>> node2children = new HashMap<N, List<Tuple<N,Integer>>>();
        GetDirectSuccessors<N, E> getChildren = new GetDirectSuccessors<N, E>();
        for(N node: postOrderTraversal.execute(network)){
            List<Tuple<N, Integer>> child2depth = new ArrayList<Tuple<N, Integer>>();
            child2depth.add(new Tuple<N,Integer>(node, 0));
            int nodeID = _node2ID.get(node);
            List<N> children = new ArrayList<N>();
            for(N child: getChildren.execute(network, node)){
                children.add(child);
                for(Tuple<N, Integer> tuple: node2children.get(child)){
                    int depth = tuple.Item2+1;
                    child2depth.add(new Tuple<N, Integer>(tuple.Item1, depth));
                    int childID = _node2ID.get(tuple.Item1);
                    int dist = _nodeDistanceMatrix[nodeID][childID];
                    if(dist==0 || dist>depth){
                        _nodeDistanceMatrix[nodeID][childID] = depth - 1;
                        _nodeDistanceMatrix[childID][nodeID] = depth - 1;
                    }
                }
            }
            if(children.size()==2){
                for(Tuple<N, Integer> tuple1: node2children.get(children.get(0))){
                    int child1ID = _node2ID.get(tuple1.Item1);
                    int depth1 = tuple1.Item2;
                    for(Tuple<N, Integer> tuple2: node2children.get(children.get(1))){
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
            for(int i=0; i<4; i++){
                if(random<_operationProbabilities[i]){
                    operationID = i;
                    break;
                }
            }

            if((operationID == 0 && !incrementHybrid) || ((operationID == 1||operationID == 2) && !decrementHybrid)){
                stop = false;
            }

        }while(!stop);
        //System.out.println();
        //System.out.println("Next Move :" + operationID);
        return operationID;
    }

    private void computeNetworkEdges(G network, ArrayList<E> allEdges, ArrayList<E> allTreeEdges, ArrayList<E> allReticulationEdges, ArrayList<E> allRemovableReticulationEdges){
        allEdges.addAll((Collection)network.getEdges());

        for (E edge : allEdges) {
            Tuple<N, N> nodesOfEdge = network.getNodesOfEdge(edge);
            GetInDegree<N, E> getInDegree = new GetInDegree<N, E>();
            if (getInDegree.execute(network, nodesOfEdge.Item2) == 2) {
                if (getInDegree.execute(network, nodesOfEdge.Item1) != 2) {
                    allRemovableReticulationEdges.add(edge);
                }
                allReticulationEdges.add(edge);
            } else {
                allTreeEdges.add(edge);
            }
        }

    }

    private boolean setNextMove(G network, ArrayList<E> allEdges, ArrayList<E> allTreeEdges, ArrayList<E> allReticulationEdges,ArrayList<E> allRemovableReticulationEdges){
        switch(_operationID){
            case 0: //setParametersForReticulationEdgeAddition(allEdges,allReticulationEdges);
                    setParametersForReticulationEdgeAddition(network, allTreeEdges);
                    break;
            case 1:
                    if(_previousTried.get(_operationID).size()==allRemovableReticulationEdges.size()){
                        return false;
                    }

                    setParametersForReticulationEdgeDeletion(network, allRemovableReticulationEdges);
                    break;
            case 2: //setParametersForReticulationEdgeDestinationChange(allRemovableReticulationEdges,allEdges);

                    if(_previousTried.get(_operationID).size()==allRemovableReticulationEdges.size()*allEdges.size()){
                        return false;
                    }
                    setParametersForReticulationEdgeDestinationChange(network, allRemovableReticulationEdges,allEdges);
                    break;
            case 3: //setParametersForEdgeSourceChange(allEdges, allEdges);

                    //int numTreeEdge = allTreeEdges.size();
                    int upper = allEdges.size()*(allEdges.size()-1);
                    if(_previousTried.get(_operationID).size()==upper){
                        return false;
                        //System.out.println(_previousTried.get(_operationID));
                        //throw new RuntimeException("tried all EdgeSourceChange");
                    }
                    setParametersForEdgeSourceChange(network, allEdges, allEdges);
                    break;
        }
        return true;
    }

    private void setParametersForReticulationEdgeAddition(G network, ArrayList<E> treeEdges){
        E sourceEdge, destinationEdge;
        int size = treeEdges.size();
        //boolean endSampling;
        //do{
            //endSampling = true;
            int source = _random.nextInt(size);
            int destination = source;
            while(source==destination){
                destination = _random.nextInt(size);
            }
            sourceEdge = treeEdges.get(source);
            destinationEdge = treeEdges.get(destination);

            /*
            PhyloEdge<String> te = new PhyloEdge("16653450","77758997");
            PhyloEdge<String> de = new PhyloEdge("6973840","2184769");
            sourceEdge = (E)te;
            destinationEdge = (E)de;
            */

            /*
            Tuple<N, N> nodesOfSourceEdge = _network.getNodesOfEdge(sourceEdge);
            Tuple<N, N> nodesOfDestinationEdge = _network.getNodesOfEdge(destinationEdge);
            if(_network.getEdge(nodesOfSourceEdge.Item1, nodesOfDestinationEdge.Item1)!=null
                    || _network.getEdge(nodesOfDestinationEdge.Item1, nodesOfSourceEdge.Item1)!=null
                    || _network.getEdge(nodesOfSourceEdge.Item1, nodesOfDestinationEdge.Item2)!=null
                    || _network.getEdge(nodesOfSourceEdge.Item2, nodesOfDestinationEdge.Item1)!=null){
                endSampling = false;
                break;
            }
            */

        //}while(!endSampling);
        //System.out.println("Add " + sourceEdge + " to " + destinationEdge);
        _networkOperators[0].setParameters(network,null, sourceEdge, destinationEdge);
        _operationEdge1 = sourceEdge;
        _operationEdge2 = destinationEdge;

    }

    private void setParametersForReticulationEdgeDeletion(G network, ArrayList<E> allRemovableReticulationEdges){
        Set<Integer> edgesTried = _previousTried.get(_operationID);
        Tuple<E,E> edgeTuple;
        int randomID;
        do{
            randomID = _random.nextInt(allRemovableReticulationEdges.size());
            E targetEdge = allRemovableReticulationEdges.get(randomID);
            edgeTuple = new Tuple<E, E>(targetEdge, targetEdge);
        }while(edgesTried.contains(randomID));
        edgesTried.add(randomID);
        //PhyloEdge<String> te = new PhyloEdge("22810764","54816791");
        //targetEdge = (E)te;

        //System.out.println("Remove " + edgeTuple.Item1);
        _networkOperators[1].setParameters(network,edgeTuple.Item1,null,null);
        _operationEdge1 = edgeTuple.Item1;
    }


    private void setParametersForReticulationEdgeDestinationChange(G network, ArrayList<E> allRemovableReticulationEdges, ArrayList<E> allEdges){

        Set<Integer> edgesTried = _previousTried.get(_operationID);
        Tuple<E,E> edgeTuple;
        int reticulationEdgeSize = allRemovableReticulationEdges.size();
        int edgeSize = allEdges.size();

        boolean endSampling;
        do{
            endSampling = true;
            int targetID = _random.nextInt(reticulationEdgeSize);
            E targetEdge = allRemovableReticulationEdges.get(targetID);
            Tuple<N,N> nodesOfTargetEdge = network.getNodesOfEdge(targetEdge);
            int node1ID = _diameterLimit==0? -1:_node2ID.get(nodesOfTargetEdge.Item2);
            /*
            N sourceEdgeItem1=null;            
            GetDirectPredecessors<N, E> findParents = new GetDirectPredecessors<N, E>();
            int count = 0;
            for (N node : findParents.execute(_network, nodesOfTargetEdge.Item1)) {
                sourceEdgeItem1 = node;
                count++;
            }
            if(count!=1 && count!=0){
                throw new RuntimeException(nodesOfTargetEdge.Item1 + " should have zero or one parent");
            }

            N sourceEdgeItem2=null;
            GetDirectSuccessors<N, E> findChildren = new GetDirectSuccessors<N, E>();
            count = 0;
            for (N node : findChildren.execute(_network, nodesOfTargetEdge.Item1)) {
                if(!node.equals(nodesOfTargetEdge.Item2)){
                    sourceEdgeItem2 = node;
                }
                count++;
            }
            if(count!=2){
                throw new RuntimeException(nodesOfTargetEdge.Item1+" should have two children");
            }
            */
            int destinationID = _random.nextInt(edgeSize);
            E destinationEdge = allEdges.get(destinationID);

            edgeTuple = new Tuple<E, E>(targetEdge, destinationEdge);
            int tupleID = (int)(Math.pow(10,new String(allEdges.size()+"").length())*targetID) + destinationID;
            if(edgesTried.contains(tupleID)){
                endSampling = false;
                continue;
            }

            Tuple<N,N> nodesOfDestinationEdge = network.getNodesOfEdge(destinationEdge);
            if(_diameterLimit!=0 && _nodeDistanceMatrix[node1ID][_node2ID.get(nodesOfDestinationEdge.Item2)]>_diameterLimit){
                edgesTried.add(tupleID);
                endSampling = false;
                continue;
            }

            edgesTried.add(tupleID);
            //System.out.println("trying redirect " + targetEdge + " to " + destinationEdge);


            if(nodesOfDestinationEdge.Item2.equals(nodesOfTargetEdge.Item1) || nodesOfDestinationEdge.Item1.equals(nodesOfTargetEdge.Item1)
                    || nodesOfDestinationEdge.Item2.equals(nodesOfTargetEdge.Item2) || nodesOfDestinationEdge.Item1.equals(nodesOfTargetEdge.Item2)){
                //edgesTried.add(edgeTuple);
                endSampling = false;
                continue;
            }

            /*
            if(new GetInDegree<N,E>().execute(_network,nodesOfDestinationEdge.Item1)==2 && _network.getEdge(nodesOfTargetEdge.Item1,nodesOfDestinationEdge.Item1)!=null){
                endSampling = false;
                continue;
            }

            if(endSampling && sourceEdgeItem1!=null){
                if(_network.getEdge(sourceEdgeItem1, nodesOfDestinationEdge.Item1)!=null
                        || _network.getEdge(nodesOfDestinationEdge.Item1, sourceEdgeItem1)!=null
                        || _network.getEdge(sourceEdgeItem1, nodesOfDestinationEdge.Item2)!=null
                        || _network.getEdge(sourceEdgeItem2, nodesOfDestinationEdge.Item1)!=null
                        || _network.getEdge(sourceEdgeItem2, nodesOfDestinationEdge.Item2)!=null
                        || _network.getEdge(nodesOfDestinationEdge.Item2, sourceEdgeItem2)!=null){
                    endSampling = false;
                }
            }
            */

        }while(!endSampling);

        //System.out.println("Redirect " + edgeTuple.Item1 + " to " + edgeTuple.Item2);
        /*
        Tuple<N,N> nodesofTarget = _network.getNodesOfEdge(targetEdge);
        GetDirectPredecessors<N,E> findParents = new GetDirectPredecessors<N,E>();
        for(N node : findParents.execute(_network, nodesofTarget.Item2)){
            System.out.println(node);
        }
        GetInDegree<N,E> getInDegree = new GetInDegree<N, E>();
        System.out.println(getInDegree.execute(_network,nodesofTarget.Item2));
        */
        _networkOperators[2].setParameters(network, edgeTuple.Item1, null, edgeTuple.Item2);
        _operationEdge1 = edgeTuple.Item1;
        _operationEdge2 = edgeTuple.Item2;

    }

    private void setParametersForEdgeSourceChange(G network, ArrayList<E> allEdges1, ArrayList<E> allEdges2){
        int allEdgeSize = allEdges1.size();
        //int treeEdgeSize = allTreeEdges.size();
        Set<Integer> edgesTried = _previousTried.get(_operationID);
        Tuple<E,E> edgeTuple;

        /*
        PhyloEdge<String> te = new PhyloEdge("1509113","22283294");
        PhyloEdge<String> de = new PhyloEdge("36951054","E");
        targetEdge = (E)te;
        destinationEdge = (E)de;
        */

        boolean endSampling;
        do{
            endSampling = true;
            int target = _random.nextInt(allEdgeSize);
            E targetEdge = allEdges1.get(target);
            Tuple<N,N> nodesOfTargetEdge = network.getNodesOfEdge(targetEdge);
            int node1ID = _diameterLimit==0? -1:_node2ID.get(nodesOfTargetEdge.Item2);
            GetInDegree<N,E> getInDegree = new GetInDegree<N, E>();
            //boolean targetIsReticulation = getInDegree.execute(network,nodesOfTargetEdge.Item2)==2;
            //int size = targetIsReticulation ? treeEdgeSize:allEdgeSize;
            int size = allEdgeSize;
            int destination = target;
            while(target==destination){
                destination = _random.nextInt(size);
            }
            //E destinationEdge = targetIsReticulation? allTreeEdges.get(destination):allEdges.get(destination);
            E destinationEdge = allEdges1.get(destination);
            edgeTuple = new Tuple<E, E>(targetEdge, destinationEdge);
            //System.out.println(new String(allEdges1.size()+"").length());
            int tupleID = (int)(Math.pow(10,new String(allEdgeSize+"").length()))*target + destination;
            if(edgesTried.contains(tupleID)){
                endSampling = false;
                continue;
            }

            Tuple<N,N> nodesOfDestinationEdge = network.getNodesOfEdge(destinationEdge);
            if(_diameterLimit!=0 && _nodeDistanceMatrix[node1ID][_node2ID.get(nodesOfDestinationEdge.Item2)]>_diameterLimit){
                edgesTried.add(tupleID);
                endSampling = false;
                continue;
            }

            edgesTried.add(tupleID);
            //System.out.println("add trying1 "+edgeTuple);

            //
            int targetInDegree = getInDegree.execute(network,nodesOfTargetEdge.Item1);
            if(targetInDegree==2){
                //int targetOutDegree = getInDegree.execute(network,nodesOfTargetEdge.Item2);
                int index = 0;
                for(E edge: allEdges2){
                    if(!edge.equals(targetEdge))
                        edgesTried.add((int)(Math.pow(10,new String(allEdgeSize+"").length()))*target + index++);
                    //System.out.println("add trying3 "+new Tuple<E, E>(targetEdge,edge));
                }
                /*
                if(targetOutDegree==1){
                    for(E edge: allEdges1){
                        if(!edge.equals(targetEdge))
                            edgesTried.add(new Tuple<E, E>(targetEdge,edge));
                            System.out.println("add trying2 "+new Tuple<E, E>(targetEdge,edge));
                    }
                }
                else{
                    for(E edge: allEdges2){
                        if(!edge.equals(targetEdge))
                            edgesTried.add(new Tuple<E, E>(targetEdge,edge));
                        System.out.println("add trying3 "+new Tuple<E, E>(targetEdge,edge));
                    }
                }
                */
                endSampling = false;
                continue;
            }


            if(targetInDegree == 0){
                N targetEdgeSibling=null;
                int count = 0;
                for (N node : new GetDirectSuccessors<N, E>().execute(network, nodesOfTargetEdge.Item1)) {
                    if (!node.equals(nodesOfTargetEdge.Item2)) {
                        targetEdgeSibling = node;
                    }
                    count++;
                }
                if(count!=2){
                    throw new RuntimeException(nodesOfTargetEdge.Item1+" should have two children!");
                }
                if(getInDegree.execute(network,targetEdgeSibling)==2){
                    int index = 0;
                    for(E edge: allEdges2){
                        if(!edge.equals(targetEdge))
                            edgesTried.add((int)(Math.pow(10,new String(allEdgeSize+"").length()))*target + index++);
                            //System.out.println("add trying4 "+new Tuple<E, E>(targetEdge,edge));
                    }
                    endSampling = false;
                    continue;
                }
            }

            if(nodesOfTargetEdge.Item1.equals(nodesOfDestinationEdge.Item1) || nodesOfTargetEdge.Item1.equals(nodesOfDestinationEdge.Item2) ||
                    nodesOfTargetEdge.Item2.equals(nodesOfDestinationEdge.Item1) || nodesOfTargetEdge.Item2.equals(nodesOfDestinationEdge.Item2)){
                endSampling = false;
                continue;
            }

            /*
            if(targetIsReticulation){
                N sourceEdgeItem1=nodesOfDestinationEdge.Item1;
                N sourceEdgeItem2=nodesOfDestinationEdge.Item2;
                N destinationEdgeItem1 = null;
                GetDirectPredecessors<N, E> findParents = new GetDirectPredecessors<N, E>();
                int count = 0;
                for (N node : findParents.execute(_network, nodesOfTargetEdge.Item2)) {
                    if(node.equals(nodesOfTargetEdge.Item1)){
                        destinationEdgeItem1 = node;
                    }
                    count++;
                }
                if(count!=2){
                    throw new RuntimeException(nodesOfTargetEdge.Item2 + " should two parents");
                }

                N destinationEdgeItem2 = null;
                GetDirectSuccessors<N, E> findChildren = new GetDirectSuccessors<N, E>();
                count = 0;
                for (N node : findChildren.execute(_network, nodesOfTargetEdge.Item2)) {
                    destinationEdgeItem2 = node;
                    count++;
                }
                if(count!=1){
                    throw new RuntimeException(nodesOfTargetEdge.Item1+" should have two children");
                }

                if(endSampling && sourceEdgeItem1!=null){
                    if(_network.getEdge(sourceEdgeItem1, destinationEdgeItem1)!=null
                            || _network.getEdge(destinationEdgeItem1, sourceEdgeItem1)!=null
                            || _network.getEdge(sourceEdgeItem1, destinationEdgeItem2)!=null
                            || _network.getEdge(sourceEdgeItem2, destinationEdgeItem1)!=null

                            || _network.getEdge(sourceEdgeItem2, destinationEdgeItem2)!=null
                            || _network.getEdge(destinationEdgeItem2, sourceEdgeItem2)!=null){
                        endSampling = false;
                        break;
                    }
                }
            }
            */
        }while(!endSampling);

        //System.out.println("Redirect " + edgeTuple.Item1 + " to " + edgeTuple.Item2);

        _networkOperators[3].setParameters(network,edgeTuple.Item1,null,edgeTuple.Item2);
        _operationEdge1 = edgeTuple.Item1;
        _operationEdge2 = edgeTuple.Item2;
    }



}