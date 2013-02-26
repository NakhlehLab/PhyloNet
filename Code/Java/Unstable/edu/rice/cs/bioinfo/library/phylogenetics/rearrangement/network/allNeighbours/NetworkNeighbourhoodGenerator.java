package edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.allNeighbours;


import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.NetworkValidatorBase;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func3;
import edu.rice.cs.bioinfo.library.programming.Func4;

/**
 * Created by IntelliJ IDEA.
 * User: yy9
 * Date: 6/7/12
 * Time: 9:26 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class NetworkNeighbourhoodGenerator<G,N,E> extends NetworkValidatorBase<N,E> {
    protected final NetworkRearrangementOperation[] _networkOperators;


    public NetworkNeighbourhoodGenerator(Func1<G, N> makeNode, Func3<G, N, N, E> makeEdge)
    {
        _networkOperators = new NetworkRearrangementOperation[4];
        _networkOperators[0] = new ReticulationEdgeAddition(makeNode, makeEdge);
        _networkOperators[1] = new ReticulationEdgeDeletion(makeNode, makeEdge);
        _networkOperators[2] = new ReticulationEdgeDestinationChange(makeNode, makeEdge);
        _networkOperators[3] = new EdgeSourceChange(makeNode, makeEdge);
        
    }


   public G performRearrangement(G network, Integer operation, E edge1, E edge2){
       if(operation==0){
           _networkOperators[operation].setParameters(network, null, edge1, edge2);
       }
       else if(operation==1){
           _networkOperators[operation].setParameters(network, edge1, null, null);
       }
       else if(operation==2 || operation==3){
           _networkOperators[operation].setParameters(network, edge1, null, edge2);
       }
       _networkOperators[operation].performOperation();
       return network;
   }

    public G computeRandomNeighbour(G network, boolean incrementHybrid, Func4<G,Integer,E,E,Boolean> rearrangementComputed){
        return null;
    }

}