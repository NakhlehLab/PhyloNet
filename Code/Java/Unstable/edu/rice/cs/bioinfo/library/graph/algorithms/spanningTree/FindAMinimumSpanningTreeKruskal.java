package edu.rice.cs.bioinfo.library.graph.algorithms.spanningTree;

import edu.rice.cs.bioinfo.library.programming.Tuple;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/12/13
 * Time: 1:44 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class FindAMinimumSpanningTreeKruskal<G,E,W extends Comparable<W>> implements FindAMinimumSpanningTree<G,E,W>
{
    public Tuple<Set<E>,W> execute(final G graph) throws GraphDisconnectedException
    {
        Set<?> graphNodes = getNodes(graph);

        if(graphNodes.size() == 0)
            return new Tuple<Set<E>, W>(new HashSet<E>(), makeZero());

        Map<Object,Set<Object>> nodeToSpanSet = new HashMap<Object,Set<Object>>();


        for(Object obj : graphNodes)
        {
            HashSet<Object> spanSet = new HashSet<Object>();
            spanSet.add(obj);
            nodeToSpanSet.put(obj, spanSet);
        }

        int distinctSpanSetsCount = graphNodes.size();

        LinkedList<E> edgesByWeightAscending = new LinkedList<E>(getEdges(graph));
        Collections.sort(edgesByWeightAscending, new Comparator<E>()
        {
            public int compare(E edge1, E edge2)
            {
                W weight1 = getWeight(edge1, graph);
                W weight2 = getWeight(edge2, graph);
                return weight1.compareTo(weight2);
            }
        } );

        if(edgesByWeightAscending.size() == 0)
            throw new GraphDisconnectedException();

        Set<E> spanTreeEdges = new HashSet<E>();

        Iterator<E> edgesByWeightAscendingElements = edgesByWeightAscending.iterator();
        while(distinctSpanSetsCount > 1)
        {
           E edge = edgesByWeightAscendingElements.next();

            Tuple<?,?> nodesOfEdge = getNodesOfEdge(edge, graph);
            Set<Object> setOfNode1 = nodeToSpanSet.get(nodesOfEdge.Item1);
            Set<Object> setOfNode2 = nodeToSpanSet.get(nodesOfEdge.Item2);

            if(setOfNode1 == setOfNode2)
                continue;

            spanTreeEdges.add(edge);

            setOfNode1.addAll(setOfNode2);

            for(Object node : setOfNode2)
            {
                nodeToSpanSet.put(node, setOfNode1);
            }

            distinctSpanSetsCount--;
        }

        if(nodeToSpanSet.get(graphNodes.iterator().next()).size() != graphNodes.size())
            throw new GraphDisconnectedException();

        W totalWeight = sumWeights(spanTreeEdges, graph);
        return new Tuple<Set<E>, W>(spanTreeEdges, totalWeight);

    }

    private W sumWeights(Set<E> edges, G graph)
    {
        W accum = makeZero();

        for(E edge : edges)
        {
            accum = add(accum, getWeight(edge, graph));
        }

        return accum;
    }

    protected abstract W add(W term1, W term2);

    protected abstract W makeZero();

    protected abstract Tuple<?,?> getNodesOfEdge(E edge, G graph);

    protected abstract W getWeight(E edge, G graph);

    protected abstract Set<? extends E> getEdges(G graph);

    protected abstract Set<?> getNodes(G graph);
}
