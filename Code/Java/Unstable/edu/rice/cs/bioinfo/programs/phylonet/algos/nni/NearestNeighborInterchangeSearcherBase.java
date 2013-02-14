package edu.rice.cs.bioinfo.programs.phylonet.algos.nni;

import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.UnrootedTreespaceHeuristicSearcherBase;

import java.util.Iterator;

public abstract class NearestNeighborInterchangeSearcherBase<T, N, E, S> extends UnrootedTreespaceHeuristicSearcherBase<T, N, E, S>
{


    public NearestNeighborInterchangeSearcherBase(Predicate1<T> isRooted, Func3<T, N, E, Boolean> isDestinationNode, Func1<T, Iterable<N>> getNodes, Func1<T, Iterable<E>> getEdges, Func2<T, N, Iterable<E>> getIncidentEdges, Func2<T, E, Tuple<N, N>> getNodesOfEdge)
    {
    	super(isRooted, isDestinationNode, getNodes, getEdges, getIncidentEdges, getNodesOfEdge);
	
    }

    protected /*override*/ abstract T searchWithoutTreeValidation(T tree, Func1<T, S> getTreeScore, Func2<S, S, Boolean> isBetterScore, long maxIterations, S givenTreeScore);

    @Override
    public /*override*/ void assertValidTree(T tree, Predicate1<T> isRooted, Func1<T, Iterable<N>> getNodes, Func1<T, Iterable<E>> getEdges, Func2<T, N, Iterable<E>> getIncidentEdges, Func2<T, E, Tuple<N, N>> getNodesOfEdge)
    {
        super.assertValidTree(tree, isRooted, getNodes, getEdges, getIncidentEdges, getNodesOfEdge);

        for (N node : getNodes.execute(tree))
        {
            Iterator<E> incidentEdges = getIncidentEdges.execute(tree, node).iterator();
            incidentEdges.next();
            if (incidentEdges.hasNext())
            {
                return;
            }
        }

        throw new IllegalArgumentException("Given tree does not contain an internal edge.");
    }
}
