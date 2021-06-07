package edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.allNeighbours;

import edu.rice.cs.bioinfo.library.phylogenetics.GetDirectSuccessors;
import edu.rice.cs.bioinfo.library.phylogenetics.GetInDegree;
import edu.rice.cs.bioinfo.library.phylogenetics.GetNodesPostOrder;
import edu.rice.cs.bioinfo.library.phylogenetics.Graph;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func3;
import edu.rice.cs.bioinfo.library.programming.Func4;
import edu.rice.cs.bioinfo.library.programming.Tuple;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: yy9
 * Date: 7/2/12
 * Time: 2:04 PM
 * To change this template use File | Settings | File Templates.
 */
public class NetworkWholeNeighbourhoodGenerator<G extends Graph<N,E>,N,E> extends NetworkNeighbourhoodGenerator<G,N,E> {

    public NetworkWholeNeighbourhoodGenerator(Func1<G, N> makeNode, Func3<G, N, N, E> makeEdge)
    {
        super(makeNode, makeEdge);

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

    private int[][] computeNodeDistances(G network, HashMap<N, Integer> node2ID){
        int[][] nodeDistanceMatrix = new int[node2ID.size()][node2ID.size()];
        GetNodesPostOrder<N, E> postOrderTraversal = new GetNodesPostOrder<N, E>();
        HashMap<N, List<Tuple<N,Integer>>> node2children = new HashMap<N, List<Tuple<N,Integer>>>();
        GetDirectSuccessors<N, E> getChildren = new GetDirectSuccessors<N, E>();
        for(N node: postOrderTraversal.execute(network)){
            List<Tuple<N, Integer>> child2depth = new ArrayList<Tuple<N, Integer>>();
            child2depth.add(new Tuple<N,Integer>(node, 0));
            int nodeID = node2ID.get(node);
            List<N> children = new ArrayList<N>();
            for(N child: getChildren.execute(network, node)){
                children.add(child);
                for(Tuple<N, Integer> tuple: node2children.get(child)){
                    int depth = tuple.Item2+1;
                    child2depth.add(new Tuple<N, Integer>(tuple.Item1, depth));
                    int childID = node2ID.get(tuple.Item1);
                    int dist = nodeDistanceMatrix[nodeID][childID];
                    if(dist==0 || dist>depth){
                        nodeDistanceMatrix[nodeID][childID] = depth - 1;
                        nodeDistanceMatrix[childID][nodeID] = depth - 1;
                    }
                }
            }
            if(children.size()==2){
                for(Tuple<N, Integer> tuple1: node2children.get(children.get(0))){
                    int child1ID = node2ID.get(tuple1.Item1);
                    int depth1 = tuple1.Item2;
                    for(Tuple<N, Integer> tuple2: node2children.get(children.get(1))){
                        int child2ID = node2ID.get(tuple2.Item1);
                        int depth2 = tuple2.Item2;

                        if(child1ID == child2ID){
                            continue;
                        }
                        int dist = nodeDistanceMatrix[child1ID][child2ID];
                        int newDist = depth1 + depth2;
                        if(dist==0 || dist>newDist){
                            nodeDistanceMatrix[child1ID][child2ID] = newDist;
                            nodeDistanceMatrix[child2ID][child1ID] = newDist;
                        }
                    }
                }
            }
            node2children.put(node, child2depth);
        }
        return nodeDistanceMatrix;
    }

    public void generateHorizontalNeighbours(G network, Func4<G,Integer,E,E,Boolean> rearrangementComputed){
        generateHorizontalNeighbours(network, rearrangementComputed, 0);
    }

    public void generateHorizontalNeighbours(G network, Func4<G,Integer,E,E,Boolean> rearrangementComputed, int diameterLimit){
        //System.out.println("in Horizontal");

        ArrayList<E> allEdges = new ArrayList<E>();
        ArrayList<E> allReticulationEdges = new ArrayList<E>();
        ArrayList<E> allRemovableReticulationEdges = new ArrayList<E>();
        ArrayList<E> allTreeEdges = new ArrayList<E>();
        computeNetworkEdges(network, allEdges, allTreeEdges, allReticulationEdges, allRemovableReticulationEdges);

        int[][] nodeDistanceMatrix = null;
        HashMap<N, Integer> node2ID = null;
        if(diameterLimit!=0){
            node2ID = new HashMap<N, Integer>();
            int id = 0;
            for(N node: network.getNodes()){
                node2ID.put(node, id++);
            }
            nodeDistanceMatrix = computeNodeDistances(network, node2ID);
        }

        for(E targetEdge: allReticulationEdges){
            Tuple<N,N> nodesOfTargetEdge = network.getNodesOfEdge(targetEdge);
            int node1ID = diameterLimit==0? -1:node2ID.get(nodesOfTargetEdge.Item2);
            for(E destinationEdge: allEdges){
                if(targetEdge.equals(destinationEdge))continue;
                Tuple<N,N> nodesOfDestinationEdge = network.getNodesOfEdge(destinationEdge);
                if(diameterLimit!=0 && nodeDistanceMatrix[node1ID][node2ID.get(nodesOfDestinationEdge.Item2)]>diameterLimit){
                    continue;
                }
                if(nodesOfDestinationEdge.Item2.equals(nodesOfTargetEdge.Item1) || nodesOfDestinationEdge.Item1.equals(nodesOfTargetEdge.Item1)
                        || nodesOfDestinationEdge.Item2.equals(nodesOfTargetEdge.Item2) || nodesOfDestinationEdge.Item1.equals(nodesOfTargetEdge.Item2)){
                    continue;
                }
                _networkOperators[2].setParameters(network, targetEdge, null, destinationEdge);
                boolean findValidNetwork = true;
                try{
                    _networkOperators[2].performOperation();
                    assertValidNetwork(network);
                }catch(IllegalArgumentException e){
                    _networkOperators[2].undoOperation();
                    findValidNetwork = false;
                }catch(IllegalStateException e){
                    findValidNetwork = false;
                }
                if(!findValidNetwork){
                    continue;
                }
                boolean ifContinue = rearrangementComputed.execute(network,2,targetEdge,destinationEdge);
                _networkOperators[2].undoOperation();
                if(!ifContinue){
                    return;
                }
            }
        }

        GetInDegree<N,E> getInDegree = new GetInDegree<N, E>();
        for(E targetEdge: allEdges){
            Tuple<N,N> nodesOfTargetEdge = network.getNodesOfEdge(targetEdge);
            int targetInDegree = getInDegree.execute(network,nodesOfTargetEdge.Item1);

            if(targetInDegree==2){
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
                    continue;
                }
            }
            int node1ID = diameterLimit==0? -1:node2ID.get(nodesOfTargetEdge.Item1);
            for(E destinationEdge: allEdges){
                Tuple<N,N> nodesOfDestinationEdge = network.getNodesOfEdge(destinationEdge);
                if(diameterLimit!=0 && nodeDistanceMatrix[node1ID][node2ID.get(nodesOfDestinationEdge.Item2)]>diameterLimit){
                    continue;
                }
                if(nodesOfTargetEdge.Item1.equals(nodesOfDestinationEdge.Item1) || nodesOfTargetEdge.Item1.equals(nodesOfDestinationEdge.Item2) ||
                        nodesOfTargetEdge.Item2.equals(nodesOfDestinationEdge.Item1) || nodesOfTargetEdge.Item2.equals(nodesOfDestinationEdge.Item2)){
                    continue;
                }
                _networkOperators[3].setParameters(network, targetEdge, null, destinationEdge);
                boolean findValidNetwork = true;
                try{
                    _networkOperators[3].performOperation();
                    assertValidNetwork(network);
                }catch(IllegalArgumentException e){
                    _networkOperators[3].undoOperation();
                    findValidNetwork = false;
                }catch(IllegalStateException e){
                    findValidNetwork = false;
                }
                if(!findValidNetwork){
                    continue;
                }
                //System.out.println(targetEdge + " to " + destinationEdge);
                //boolean ifContinue = true;
                boolean ifContinue = rearrangementComputed.execute(network,3,targetEdge,destinationEdge);
                _networkOperators[3].undoOperation();
                //System.out.println(network);
                if(!ifContinue){
                    return;
                }
            }
        }
    }

    public void generateVerticalNeighbours(G network, Func4<G,Integer,E,E,Boolean> rearrangementComputed){
        ArrayList<E> allEdges = new ArrayList<E>();
        ArrayList<E> allReticulationEdges = new ArrayList<E>();
        ArrayList<E> allRemovableReticulationEdges = new ArrayList<E>();
        ArrayList<E> allTreeEdges = new ArrayList<E>();
        computeNetworkEdges(network, allEdges, allTreeEdges, allReticulationEdges, allRemovableReticulationEdges);

        for(E sourceEdge: allEdges){
            for(E destinationEdge: allEdges){
            //for(int i=index+1; i<allEdges.size(); i++){
                //E destinationEdge = allEdges.get(i);
                if(sourceEdge.equals(destinationEdge)){
                    continue;
                }
                _networkOperators[0].setParameters(network, null, sourceEdge, destinationEdge);
                boolean findValidNetwork = true;
                try{
                    _networkOperators[0].performOperation();
                    assertValidNetwork(network);
                }catch(IllegalArgumentException e){
                    _networkOperators[0].undoOperation();
                    findValidNetwork = false;
                }catch(IllegalStateException e){
                    findValidNetwork = false;
                }
                if(!findValidNetwork){
                    continue;
                }
                boolean ifContinue = rearrangementComputed.execute(network,0,sourceEdge,destinationEdge);
                _networkOperators[0].undoOperation();
                if(!ifContinue){
                    return;
                }
            }
        }
    }

}