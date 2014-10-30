package edu.rice.cs.bioinfo.programs.phylonet.structs.network.rearrangement;



import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;

/**
 * Created by IntelliJ IDEA.
 * User: yy9
 * Date: 6/7/12
 * Time: 9:26 PM
 * To change this template use File | Settings | File Templates.
 */
public interface NetworkNeighbourhoodGenerator{

   public abstract void performRearrangement(Network network, Integer operation, Tuple<NetNode,NetNode> targetEdge, Tuple<NetNode,NetNode> sourceEdge, Tuple<NetNode,NetNode> destinationEdge);

   public abstract void undo();


   public abstract void computeRandomNeighbour(Network network);

}