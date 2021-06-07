package edu.rice.cs.bioinfo.library.graph.algorithms.spanningTree;

import edu.rice.cs.bioinfo.library.programming.Proc1;
import edu.rice.cs.bioinfo.library.programming.Tuple;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/12/13
 * Time: 2:29 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class FindAllMinumumSpanningTreesBacktrack<G,N,E,W extends Comparable<W>>
{
    private Set<Proc1<Set<E>>> _minSpanTreeListeners = new HashSet<Proc1<Set<E>>>();

    public void addMinSpanTreeFoundListener(Proc1<Set<E>> listener)
    {
        _minSpanTreeListeners.add(listener);
    }


    public void execute(G graph) throws GraphDisconnectedException
    {
        Tuple<Set<E>,W> minSpanTree = new FindAMinimumSpanningTreeKruskal<G,E,W>()
        {
            @Override
            protected W add(W term1, W term2)
            {
                return FindAllMinumumSpanningTreesBacktrack.this.add(term1, term2);  //To change body of implemented methods use File | Settings | File Templates.
            }

            @Override
            protected W makeZero()
            {
                return FindAllMinumumSpanningTreesBacktrack.this.makeZero();  //To change body of implemented methods use File | Settings | File Templates.
            }

            @Override
            protected Tuple<?, ?> getNodesOfEdge(E edge, G graph)
            {
                return FindAllMinumumSpanningTreesBacktrack.this.getNodesOfEdge(edge, graph);  //To change body of implemented methods use File | Settings | File Templates.
            }

            @Override
            protected W getWeight(E edge, G graph)
            {
                return FindAllMinumumSpanningTreesBacktrack.this.getWeight(edge);  //To change body of implemented methods use File | Settings | File Templates.
            }

            @Override
            protected Set<? extends E> getEdges(G graph)
            {
                return FindAllMinumumSpanningTreesBacktrack.this.getEdgesOfGraph(graph);  //To change body of implemented methods use File | Settings | File Templates.
            }

            @Override
            protected Set<?> getNodes(G graph)
            {
                return FindAllMinumumSpanningTreesBacktrack.this.getNodesOfGraph(graph);  //To change body of implemented methods use File | Settings | File Templates.
            }
        }.execute(graph);

        final W minSpanWeight = minSpanTree.Item2;

        BacktrackSpanningTrees<G,N,E> backtracking = new BacktrackSpanningTrees<G,N,E>()
        {
            @Override
            protected boolean exploreSubtree(Set<E> partialSpanTree)
            {
                W treeWeightAccum = makeZero();
                for(E edge : partialSpanTree)
                {
                    W edgeWeight = getWeight(edge);
                    treeWeightAccum = add(treeWeightAccum, edgeWeight);
                }

                return treeWeightAccum.compareTo(minSpanWeight) <= 0;
            }

            @Override
            protected Set<N> getNodesOfGraph(G graph)
            {
                return FindAllMinumumSpanningTreesBacktrack.this.getNodesOfGraph(graph);  //To change body of implemented methods use File | Settings | File Templates.
            }

            @Override
            protected Set<E> performGetEdgesOfGraph(G graph)
            {
                TreeSet<E> edgeSet = new TreeSet<E>(new Comparator<E>()
                {
                    public int compare(E edge1, E edge2)
                    {
                        W weight1 = getWeight(edge1);
                        W weight2 = getWeight(edge2);
                        return weight1.compareTo(weight2);
                    }
                });
                return new HashSet<E>(getEdgesOfGraph(graph))
                {
                    private HashSet<E> self = this;

                    @Override
                    public Iterator<E> iterator()
                    {

                        final LinkedList<E> order = new LinkedList<E>();
                        Iterator<E> elements = super.iterator();
                        while(elements.hasNext())
                            order.add(elements.next());

                        Collections.sort(order, new Comparator<E>()
                        {
                            public int compare(E o1, E o2)
                            {
                                return getWeight(o1).compareTo(getWeight(o2));
                            }
                        });
                        return new Iterator<E>()
                        {
                            Iterator<E> inner = order.iterator();

                            private E _lastReturned = null;

                            public boolean hasNext()
                            {
                                return inner.hasNext();  //To change body of implemented methods use File | Settings | File Templates.
                            }

                            public E next()
                            {
                                _lastReturned = inner.next();  //To change body of implemented methods use File | Settings | File Templates.
                                return _lastReturned;
                            }

                            public void remove()
                            {
                                inner.remove();
                                self.remove(_lastReturned);
                            }
                        };

                    }
                };

            }

            @Override
            protected Set<E> getEdgesOfGraph(G graph)
            {
                return FindAllMinumumSpanningTreesBacktrack.this.getEdgesOfGraph(graph);  //To change body of implemented methods use File | Settings | File Templates.
            }

            @Override
            protected Iterable<? extends E> getIncidentEdges(N node, G graph)
            {
                return FindAllMinumumSpanningTreesBacktrack.this.getIncidentEdges(node, graph);  //To change body of implemented methods use File | Settings | File Templates.
            }

            @Override
            protected Tuple<? extends N, ? extends N> getNodesOfEdge(E edge, G graph)
            {
                return FindAllMinumumSpanningTreesBacktrack.this.getNodesOfEdge(edge, graph);  //To change body of implemented methods use File | Settings | File Templates.
            }
        };

        backtracking.addSpanningTreeFoundListener(new Proc1<Set<E>>()
        {
            public void execute(Set<E> minSpanTree)
            {
                fireMinSpanTreeFound(minSpanTree);
            }
        });
        backtracking.execute(graph);
    }

    private void fireMinSpanTreeFound(Set<E> minSpanTree)
    {
        for(Proc1<Set<E>> listener : _minSpanTreeListeners)
        {
            listener.execute(minSpanTree);
        }
    }

    protected abstract Set<N> getNodesOfGraph(G graph);

    protected abstract Set<E> getEdgesOfGraph(G graph);

    protected abstract Iterable<? extends E> getIncidentEdges(N node, G graph);

    protected abstract Tuple<? extends N, ? extends N> getNodesOfEdge(E edge, G graph);

    protected abstract W add(W term1, W term2);

    protected abstract W makeZero();

    protected abstract W getWeight(E edge1);
}
