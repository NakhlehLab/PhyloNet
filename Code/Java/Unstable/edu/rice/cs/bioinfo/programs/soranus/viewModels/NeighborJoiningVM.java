package edu.rice.cs.bioinfo.programs.soranus.viewModels;

import edu.rice.cs.bioinfo.library.programming.Tuple;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/5/13
 * Time: 6:07 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class NeighborJoiningVM<N,E> implements DocumentVM
{

    public final Set<E> Edges;

    public final Set<N> Nodes;

    private final Map<N,Integer> _nodeDegree = new HashMap<N,Integer>();

    public NeighborJoiningVM(Set<N> nodes, Set<E> edges) {
        Edges = edges;
        Nodes = nodes;

        for(N node : Nodes)
        {
            _nodeDegree.put(node, 0);
        }

        for(E edge : Edges)
        {
            Tuple<N,N> nodesOfEdge = getNodesOfEdge(edge);

            Integer node1Degree = _nodeDegree.get(nodesOfEdge.Item1);
            Integer node2Degree = _nodeDegree.get(nodesOfEdge.Item2);

            _nodeDegree.put(nodesOfEdge.Item1, node1Degree+1);
            _nodeDegree.put(nodesOfEdge.Item2, node2Degree+1);

        }
    }

    public abstract Tuple<N,N> getNodesOfEdge(E edge);

    public abstract String getNodeLabel(N node);

    public <R,E extends Exception> R execute(DocumentVMAlgo<R,E> algo) throws E
    {
        return algo.forNeighborJoiningVM(this);
    }


    public boolean isLeaf(N node)
    {
        return  _nodeDegree.get(node).equals(1);
    }

    public abstract String getEdgeLabel(E edge);
}
