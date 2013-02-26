package edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.allNeighbours;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func3;

/**
 * Created by IntelliJ IDEA.
 * User: yy9
 * Date: 6/7/12
 * Time: 10:43 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class NetworkRearrangementOperation<G,N,E>{
    protected G _network;
    protected E _targetEdge;
    protected E _sourceEdge;
    protected E _destinationEdge;
    protected Func3<G, N, N, E> _makeEdge;
    protected Func1<G, N> _makeNode;


    
    public NetworkRearrangementOperation(Func1<G, N> makeNode, Func3<G, N, N, E> makeEdge){
        _makeNode = makeNode;
        _makeEdge = makeEdge;
    }

    //abstract public void updateNode2Ancestors(Map<N,Set<N>> nodeToAncestors);

    public void setParameters(G network, E targetEdge, E sourceEdge, E destinationEdge){
        _network = network;
        _targetEdge = targetEdge;
        _sourceEdge = sourceEdge;
        _destinationEdge = destinationEdge;
    }

    abstract public G performOperation();

    abstract public G undoOperation();

}