package edu.rice.cs.bioinfo.library.graph.algorithms.spanningTree;

import edu.rice.cs.bioinfo.library.graph.algorithms.cycleDetection.CycleDetectorFromEdgesSets;
import edu.rice.cs.bioinfo.library.math.discrete.Combinations;
import edu.rice.cs.bioinfo.library.programming.Proc2;
import edu.rice.cs.bioinfo.library.programming.Tuple;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/15/13
 * Time: 1:07 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class FindAllMinimumSpanningTreesCombinationsAsc<G,N,E,W extends Comparable<W>>
        implements Iterable<Set<E>>
{
    class Edge<E,W extends Comparable<W>> implements Comparable<Edge<E,W>>
    {
        public final E Edge;

        public final W Weight;

        public Edge(E edge, W weight)
        {
            Edge = edge;
            Weight = weight;
        }

        public int compareTo(Edge<E, W> o)
        {
            return this.Weight.compareTo(o.Weight);
        }
    }

    private final G _graph;

    private W _minSpanWeight;

    private List<Edge<E,W>> _edgesByWeightAsc;

    private Set<N> _nodesOfGraph;


    public FindAllMinimumSpanningTreesCombinationsAsc(G graph) throws GraphDisconnectedException
    {
        _graph = graph;
    }

    protected void setup() throws GraphDisconnectedException
    {
        Tuple<Set<E>,W> minSpanTree = new FindAMinimumSpanningTreeKruskal<G,E,W>()
        {
            @Override
            protected W add(W term1, W term2)
            {
                return FindAllMinimumSpanningTreesCombinationsAsc.this.add(term1, term2);  //To change body of implemented methods use File | Settings | File Templates.
            }

            @Override
            protected W makeZero()
            {
                return FindAllMinimumSpanningTreesCombinationsAsc.this.makeZero();  //To change body of implemented methods use File | Settings | File Templates.
            }

            @Override
            protected Tuple<?, ?> getNodesOfEdge(E edge, G graph)
            {
                return FindAllMinimumSpanningTreesCombinationsAsc.this.getNodesOfEdge(edge, graph);  //To change body of implemented methods use File | Settings | File Templates.
            }

            @Override
            protected W getWeight(E edge, G graph)
            {
                return FindAllMinimumSpanningTreesCombinationsAsc.this.getWeight(edge);  //To change body of implemented methods use File | Settings | File Templates.
            }

            @Override
            protected Set<? extends E> getEdges(G graph)
            {
                return FindAllMinimumSpanningTreesCombinationsAsc.this.getEdgesOfGraph(graph);  //To change body of implemented methods use File | Settings | File Templates.
            }

            @Override
            protected Set<?> getNodes(G graph)
            {
                return FindAllMinimumSpanningTreesCombinationsAsc.this.getNodesOfGraph(graph);  //To change body of implemented methods use File | Settings | File Templates.
            }
        }.execute(_graph);

        _minSpanWeight = minSpanTree.Item2;

        _edgesByWeightAsc = new ArrayList();
        for(E edge : getEdgesOfGraph(_graph))
        {
            _edgesByWeightAsc.add(new Edge<E,W>(edge, getWeight(edge)));
        }

        Collections.sort(_edgesByWeightAsc);

        _nodesOfGraph = getNodesOfGraph(_graph);
    }

    public Iterator<Set<E>> iterator()
    {
        final Combinations<Edge<E, W>> combinations =
                new Combinations<Edge<E, W>>(_edgesByWeightAsc, _nodesOfGraph.size() -1);

        final CycleDetectorFromEdgesSets<Edge<E,W>> cycleDetector = new CycleDetectorFromEdgesSets<Edge<E, W>>()
        {
            @Override
            protected Tuple<?, ?> getNodesOfEdge(Edge<E, W> edge)
            {
                return FindAllMinimumSpanningTreesCombinationsAsc.this.getNodesOfEdge(edge.Edge, _graph);
            }
        };

        return new Iterator<Set<E>>()
        {
            private Set<E> _nextMst = null;

            private boolean _moreMstsPossible = true;

            private Iterator<Set<Edge<E, W>>> edgeCombinations = combinations.iterator();

            {
                combinations.addElementExhaustedListener(new Proc2<Integer,Set<Edge<E, W>>>()
                {
                    public void execute(Integer exhaustedIndex, Set<Edge<E, W>> edges)
                    {
                        W totalWeight = makeZero();
                        for(Edge<E,W> edge : edges)
                            totalWeight = add(totalWeight, edge.Weight);

                        if(totalWeight.compareTo(_minSpanWeight) > 0)
                            _moreMstsPossible = false;

                    }
                });
            }

            public boolean hasNext()
            {
                while(edgeCombinations.hasNext() && _moreMstsPossible)
                {
                    Set<Edge<E,W>> combination = edgeCombinations.next();
                    W totalWeight = makeZero();
                    for(Edge<E,W> edge : combination)
                        totalWeight = add(totalWeight, edge.Weight);

                    if(totalWeight.compareTo(_minSpanWeight) == 0)
                    {

                        Set<N> coveredNodes = new HashSet<N>();
                        for(Edge<E,W> edge : combination)
                        {
                            Tuple<? extends N, ? extends N> nodesOfEdge = getNodesOfEdge(edge.Edge, _graph);
                            coveredNodes.add(nodesOfEdge.Item1);
                            coveredNodes.add(nodesOfEdge.Item2);

                        }

                        if(coveredNodes.size() == _nodesOfGraph.size())
                        {
                            if(!cycleDetector.containsCycle(combination))
                            {
                                _nextMst = new HashSet<E>();
                                for(Edge<E,W> mstEdge : combination)
                                    _nextMst.add(mstEdge.Edge);
                                return true;
                            }
                        }
                    }
                }

                _nextMst = null;
                return false;
            }

            public Set<E> next()
            {
                if(_nextMst == null)
                    if(hasNext())
                        return _nextMst;
                    else
                        throw new NoSuchElementException();
                else
                {
                    return _nextMst;
                }
            }

            public void remove()
            {
                throw new UnsupportedOperationException();
            }
        };
    }


    protected abstract Set<N> getNodesOfGraph(G graph);

    protected abstract Set<E> getEdgesOfGraph(G graph);

    protected abstract Iterable<? extends E> getIncidentEdges(N node, G graph);

    protected abstract Tuple<? extends N, ? extends N> getNodesOfEdge(E edge, G graph);

    protected abstract W add(W term1, W term2);

    protected abstract W makeZero();

    protected abstract W getWeight(E edge);

}
