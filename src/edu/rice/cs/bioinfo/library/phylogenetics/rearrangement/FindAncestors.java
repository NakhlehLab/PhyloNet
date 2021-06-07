//package edu.rice.cs.bioinfo.library.phylogenetics.rearrangement;
//
//import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.tree.nni.NearestNeighborInterchange;
//import edu.rice.cs.bioinfo.library.programming.Func2;
//import edu.rice.cs.bioinfo.library.programming.Tuple;
//
//import javax.swing.text.Element;
//import java.util.HashSet;
//import java.util.LinkedList;
//import java.util.Queue;
//import java.util.Set;
//
///**
// * Created by IntelliJ IDEA.
// * User: Matt
// * Date: 5/7/12
// * Time: 6:25 PM
// * To change this template use File | Settings | File Templates.
// */
//public class FindAncestors<N,E> implements Func2<GraphReadOnly<N,E>,N,Set<N>>
//{
//    public Set<N> execute(GraphReadOnly<N, E> network, N node)
//    {
//        if(!network.isRooted())
//        {
//            throw new IllegalArgumentException("Passed graph must be rooted.");
//        }
//
//        LinkedList<N> toExpand = new LinkedList<N>();
//
//        for(N directPredecessor : getDirectPredecessors(network, node))
//        {
//            toExpand.addLast(directPredecessor);
//        }
//
//        HashSet<N> ancestors = new HashSet<N>();
//        while(toExpand.size() > 0)
//        {
//            N ancestor = toExpand.removeFirst();
//
//            if(!ancestors.contains(ancestor))
//            {
//                ancestors.add(ancestor);
//                for(N directPredecessor : getDirectPredecessors(network, ancestor))
//                {
//                    toExpand.addFirst(directPredecessor);
//                }
//            }
//        }
//
//        return ancestors;
//    }
//
//    private Iterable<N> getDirectPredecessors(GraphReadOnly<N,E> network, N node)
//    {
//        LinkedList<N> directPredecessors = new LinkedList<N>();
//
//        for(E edge : network.getIncidentEdges(node))
//        {
//            Tuple<N,N> nodesOfEdge = network.getNodesOfEdge(edge);
//            if(nodesOfEdge.Item2.equals(node))
//            {
//                directPredecessors.add(nodesOfEdge.Item1);
//            }
//        }
//
//        return directPredecessors;
//    }
//}
