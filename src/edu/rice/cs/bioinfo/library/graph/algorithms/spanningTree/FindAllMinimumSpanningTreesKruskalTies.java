package edu.rice.cs.bioinfo.library.graph.algorithms.spanningTree;

import com.google.common.collect.Collections2;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/30/13
 * Time: 4:14 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class FindAllMinimumSpanningTreesKruskalTies<G,N,E,W extends Comparable<W>> implements
        Iterable<Set<E>>
{
    class EdgeAndWeight implements Comparable<EdgeAndWeight>
    {
        public final E Edge;

        public final W Weight;

        public EdgeAndWeight(E Edge, W Weight)
        {
            this.Edge = Edge;
            this.Weight = Weight;
        }

        public int compareTo(EdgeAndWeight o)
        {
            return Weight.compareTo(o.Weight);
        }
    }

    class FindMst extends FindAMinimumSpanningTreeKruskal<G,EdgeAndWeight,W>
    {

        private List<EdgeAndWeight> _edges;

        FindMst(List<EdgeAndWeight> edges)
        {
            _edges = edges;
        }

        @Override
        protected W add(W term1, W term2)
        {
            return FindAllMinimumSpanningTreesKruskalTies.this.add(term1, term2);
        }

        @Override
        protected W makeZero()
        {
            return FindAllMinimumSpanningTreesKruskalTies.this.makeZero();
        }

        @Override
        protected Tuple<?, ?> getNodesOfEdge(EdgeAndWeight edge, G graph)
        {
            return FindAllMinimumSpanningTreesKruskalTies.this.getNodesOfEdge(edge.Edge, graph);
        }

        @Override
        protected W getWeight(EdgeAndWeight edge, G graph)
        {
            return edge.Weight;
        }

        @Override
        protected Set<? extends EdgeAndWeight> getEdges(G graph)
        {
            return new HashSet<EdgeAndWeight>(_edges);
        }

        @Override
        protected Set<?> getNodes(G graph)
        {
            return FindAllMinimumSpanningTreesKruskalTies.this.getNodesOfGraph(graph);
        }

        @Override
        protected List<EdgeAndWeight> sortEdgesAscending(final G graph)
        {
            return _edges;
        }
    }

    private final G _graph;

    public FindAllMinimumSpanningTreesKruskalTies(G graph)
    {
        _graph = graph;
    }

    public Iterator<Set<E>> iterator()
    {
        Map<N,Set<N>> nodeToConnectionSet = new HashMap<N,Set<N>>();

        for(N node : getNodesOfGraph(_graph))
        {
            nodeToConnectionSet.put(node, new HashSet<N>(Arrays.asList(node)));
        }

        Iterable<EdgeAndWeight> weightedEdges = IterableHelp.map(getEdgesOfGraph(_graph), new Func1<E, EdgeAndWeight>()
        {
            public EdgeAndWeight execute(E edge)
            {
                return new EdgeAndWeight(edge, getWeightOfEdge(edge));
            }
        });


        List<EdgeAndWeight> edgesByWeightAscending = new ArrayList<EdgeAndWeight>(IterableHelp.toList(weightedEdges));
        Collections.sort(edgesByWeightAscending);

        final List<List<EdgeAndWeight>> edgesByWeightAscendingGroupByWeight = new LinkedList<List<EdgeAndWeight>>();
        List<EdgeAndWeight> fillList = new ArrayList<EdgeAndWeight>();
        fillList.add(edgesByWeightAscending.get(0));

        for(int i = 1; i<edgesByWeightAscending.size(); i++)
        {
            EdgeAndWeight ithEdge = edgesByWeightAscending.get(i);

            if(ithEdge.Weight.compareTo(fillList.get(0).Weight) == 0)
                fillList.add(ithEdge);
            else
            {
                edgesByWeightAscendingGroupByWeight.add(fillList);
                fillList = new ArrayList<EdgeAndWeight>();
                fillList.add(ithEdge);
            }
        }
        edgesByWeightAscendingGroupByWeight.add(fillList);

        HashSet<Set<EdgeAndWeight>> found = new HashSet<Set<EdgeAndWeight>>();

        final List<Iterator<List<EdgeAndWeight>>> permutationGroups = new LinkedList<Iterator<List<EdgeAndWeight>>>();
        permutationGroups.add(Collections2.permutations(edgesByWeightAscendingGroupByWeight.get(0)).iterator());

        return new Iterator<Set<E>>()
        {
            private boolean _hasNext;

            private Set<E> _nextMst;

            private Set<Set<E>> _seenMsts = new HashSet<Set<E>>();

            {
                update();
            }

            public boolean hasNext()
            {
                return _hasNext;
            }

            public Set<E> next()
            {
                if(_hasNext)
                {
                    Set<E> tbr = _nextMst;
                    update();
                    return tbr;
                }
                else
                    throw new IllegalStateException();

            }

            private void update()
            {
                while(true)
                {
                    while(permutationGroups.size() > 0 &&
                            !permutationGroups.get(permutationGroups.size()-1).hasNext())
                    {
                        permutationGroups.remove(permutationGroups.size()-1);
                    }

                    if(permutationGroups.size() == 0)
                    {
                        _hasNext = false;
                        _nextMst = null;
                        return;
                    }


                    for(int i = permutationGroups.size(); i < edgesByWeightAscendingGroupByWeight.size(); i++)
                    {
                        permutationGroups.add(
                                Collections2.permutations(edgesByWeightAscendingGroupByWeight.get(i)).iterator());
                    }

                    List<EdgeAndWeight> edgesByWeightAscending = new LinkedList<EdgeAndWeight>();
                    for(Iterator<List<EdgeAndWeight>> group : permutationGroups)
                    {
                        while(group.hasNext())
                        {
                            edgesByWeightAscending.addAll(group.next());
                        }
                    }

                    try
                    {
                        Set<EdgeAndWeight> found = new FindMst(edgesByWeightAscending).execute(_graph).Item1;
                        Set<E> mstReturn = new HashSet<E>();
                        for(EdgeAndWeight edge : found)
                            mstReturn.add(edge.Edge);

                        if(!_seenMsts.contains(mstReturn))
                        {
                            _seenMsts.add(mstReturn);
                            _nextMst = mstReturn;
                            _hasNext = permutationGroups.get(0).hasNext();
                            return;
                        }

                    }
                    catch (GraphDisconnectedException e)
                    {
                        throw new RuntimeException(e);
                    }
                }
            }

            public void remove()
            {
                throw new UnsupportedOperationException();
            }
        };
    }



    private void search(List<List<EdgeAndWeight>> searchInstance, List<List<EdgeAndWeight>> edgesByWeightAscendingGroupByWeight,
                        int permutationGroupIndex, G graph, Set<Set<EdgeAndWeight>> foundMsts) throws GraphDisconnectedException
    {
        if(searchInstance.size() == edgesByWeightAscendingGroupByWeight.size())
        {
            List<EdgeAndWeight> edgesByWeightAscending = new LinkedList<EdgeAndWeight>();
            for(List<EdgeAndWeight> group : searchInstance)
            {
                for(EdgeAndWeight groupElement : group)
                {
                    edgesByWeightAscending.add(groupElement);
                }
            }

            foundMsts.add(new FindMst(edgesByWeightAscending).execute(graph).Item1);


        }
        else
        {
            List<EdgeAndWeight> permutationGroup = edgesByWeightAscendingGroupByWeight.get(permutationGroupIndex);

            for(List<EdgeAndWeight> permutation : Collections2.permutations(permutationGroup))
            {
                searchInstance.add(permutation);
                search(searchInstance, edgesByWeightAscendingGroupByWeight, permutationGroupIndex + 1, graph, foundMsts);
                searchInstance.remove(searchInstance.size()-1);
            }
        }

    }

    private void findMsts(G graph, Map<N, Set<N>> nodeToConnectionSet, List<EdgeAndWeight> edgesByWeightAscending, int edgeIndex,
                          int numNodesInGraph, Set<EdgeAndWeight> partialMstEdges)
    {
        W edgeIndexWeight = edgesByWeightAscending.get(edgeIndex).Weight;

        for(int i = edgeIndex; i<edgesByWeightAscending.size(); i++)
        {
            EdgeAndWeight iEdge = edgesByWeightAscending.get(i);
            if(iEdge.Weight.compareTo(edgeIndexWeight) != 0)
                break;

            tryIncludeEdge(iEdge, graph, clone(nodeToConnectionSet), edgesByWeightAscending, i, numNodesInGraph,
                    new HashSet<EdgeAndWeight>(partialMstEdges));


        }
    }

    private void tryIncludeEdge(EdgeAndWeight edge, G graph, Map<N, Set<N>> nodeToConnectionSet, List<EdgeAndWeight> edgesByWeightAscending,
                                int edgeIndex, int numNodesInGraph, Set<EdgeAndWeight> partialMstEdges)
    {
        if(partialMstEdges.size() == numNodesInGraph -1)
            return;

        Tuple<? extends N, ? extends N> nodesOfEdge = getNodesOfEdge(edge.Edge, graph);
        Set<N> connectedSetOfNode1 = nodeToConnectionSet.get(nodesOfEdge.Item1);
        Set<N> connectedSetOfNode2 = nodeToConnectionSet.get(nodesOfEdge.Item2);

        if(connectedSetOfNode1 == connectedSetOfNode2)
            return;

        partialMstEdges.add(edge);
        connectedSetOfNode1.addAll(connectedSetOfNode2);

        for(N nodeOf2 : connectedSetOfNode2)
        {
            nodeToConnectionSet.put(nodeOf2, connectedSetOfNode1);
        }

        findMsts(graph, nodeToConnectionSet,edgesByWeightAscending, edgeIndex+1,numNodesInGraph, partialMstEdges);



    }

    private Map<N, Set<N>> clone(Map<N, Set<N>> nodeToConnectionSet)
    {
        Map<N, Set<N>> cloneMap = new HashMap<N, Set<N>>();

        Map<Set<N>, Set<N>> givenSetToClone = new HashMap<Set<N>, Set<N>>();

        for(N key : nodeToConnectionSet.keySet())
        {
            Set<N> value = nodeToConnectionSet.get(key);

            Set<N> clone;
            if(givenSetToClone.containsKey(value))
                clone = givenSetToClone.get(value);
            else
            {
                clone = new HashSet<N>(value);
                givenSetToClone.put(value, clone);
            }

            cloneMap.put(key, clone);

        }

        return cloneMap;
    }

    protected abstract W add(W term1, W term2);

    protected abstract W makeZero();

    protected abstract Set<N> getNodesOfGraph(G graph);

    protected abstract Set<E> getEdgesOfGraph(G graph);

    protected abstract Iterable<? extends E> getIncidentEdges(N node, G graph);

    protected abstract Tuple<? extends N, ? extends N> getNodesOfEdge(E edge, G graph);

    protected abstract W getWeightOfEdge(E edge);
}
