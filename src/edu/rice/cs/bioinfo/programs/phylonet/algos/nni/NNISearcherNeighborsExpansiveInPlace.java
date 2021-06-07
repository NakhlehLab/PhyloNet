package edu.rice.cs.bioinfo.programs.phylonet.algos.nni;

import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;

import java.util.ArrayList;
import java.util.LinkedList;


public class NNISearcherNeighborsExpansiveInPlace<T, N, E, S> extends NearestNeighborInterchangeSearcherBase<T, N, E, S>
{
    private /*final*/ Func3<T, N, N, E> _makeEdge;
    
    private /*final*/ Proc2<T,E> _addEdge;

    private /*final*/ Proc2<T, E> _removeEdge;

    private Proc1<T> BeginningToSearchNeighborhoodDelegate;



    public NNISearcherNeighborsExpansiveInPlace(Predicate1<T> isRooted, Func3<T, N, E, Boolean> isDestinationNode, Func1<T, Iterable<N>> getNodes, Func1<T, Iterable<E>> getEdges, Func2<T, N, Iterable<E>> getIncidentEdges,
                                                Func2<T, E, Tuple<N, N>> getNodesOfEdge, Func3<T, N, N, E> makeEdge, Proc2<T,E> addEdge, Proc2<T,E> removeEdge)
    {
    	super(isRooted, isDestinationNode, getNodes, getEdges, getIncidentEdges, getNodesOfEdge);
	
        _makeEdge = makeEdge;
        _addEdge = addEdge;
        _removeEdge = removeEdge;
    }

    protected /*override*/ T searchWithoutTreeValidation(T tree, Func1<T, S> getTreeScore, Func2<S, S, Boolean> isBetterScore, long maxIterations, S givenTreeScore)
    {
        S bestSeenScore = givenTreeScore;
        N bestSeenSwapNodeA = null /* default(N)*/ ;
        N bestSeenSwapNodeB = null /* default(N)*/ ;
        N bestSeenInternalEdgeNodeA = null /* default(N)*/ ;
        N bestSeenInternalEdgeNodeB = null /* default(N)*/ ;
        E bestSeenSwapEdgeA = null /* default(E)*/ ;
        E bestSeenSwapEdgeB = null /* default(E)*/ ;

        for (int i = 0; i < maxIterations; i++)
        {
            boolean newBestSeen = false;


            LinkedList<E> edges = new LinkedList<E>(); // get a list to prevent concurrent modification exception
            for(E edge : getEdges.execute(tree))
            {
                edges.add(edge);
            }

            for (E edge : edges)
            {
                Tuple<N,N> nodesOfEdge = getNodesOfEdge.execute(tree, edge);
                int node1Degree = IterableHelp.countInt(getIncidentEdges.execute(tree, nodesOfEdge.Item1));
                int node2Degree = IterableHelp.countInt(getIncidentEdges.execute(tree, nodesOfEdge.Item2));
                if (node1Degree > 1 && node2Degree > 1) // is internal edge
                {
                    if (node1Degree == 3 && node2Degree == 3)
                    {
                        Ref<S> referenceToBestSeenScore = new Ref<S>(bestSeenScore);
                        Ref<N> referenceToBestSeenInternalEdgeNodeA = new Ref<N>(bestSeenInternalEdgeNodeA );
                        Ref<N> referenceToBestSeenInternalEdgeNodeB = new Ref<N>(bestSeenInternalEdgeNodeB);
                        Ref<N> referenceToBestSeenSwapNodeA = new Ref<N>(bestSeenSwapNodeA);
                        Ref<E> referenceToBestSeenSwapEdgeA = new Ref<E>(bestSeenSwapEdgeA);
                        Ref<N> referenceToBestSeenSwapNodeB = new Ref<N>(bestSeenSwapNodeB);
                        Ref<E> referenceToBestSeenSwapEdgeB = new Ref<E>(bestSeenSwapEdgeB);
                        Ref<Boolean> referenceToNewBestSeen = new Ref<Boolean>(newBestSeen);
                        swapDoubleDegreeThree(tree, getTreeScore, isBetterScore, givenTreeScore, edge, nodesOfEdge.Item1, nodesOfEdge.Item2, /*ref*/ referenceToBestSeenScore,
                                /*ref*/ referenceToBestSeenInternalEdgeNodeA, /*ref*/ referenceToBestSeenInternalEdgeNodeB, /*ref*/ referenceToBestSeenSwapNodeA, /*ref*/ referenceToBestSeenSwapEdgeA, /*ref*/ referenceToBestSeenSwapNodeB, /*ref*/ referenceToBestSeenSwapEdgeB,
                                /*ref*/ referenceToNewBestSeen);
                        bestSeenScore = referenceToBestSeenScore.get();
                        bestSeenInternalEdgeNodeA = referenceToBestSeenInternalEdgeNodeA.get();
                        bestSeenInternalEdgeNodeB = referenceToBestSeenInternalEdgeNodeB.get();
                        bestSeenSwapNodeA = referenceToBestSeenSwapNodeA.get();
                        bestSeenSwapEdgeA = referenceToBestSeenSwapEdgeA.get();
                        bestSeenSwapNodeB = referenceToBestSeenSwapNodeB.get();
                        bestSeenSwapEdgeB = referenceToBestSeenSwapEdgeB.get();
                        newBestSeen = referenceToNewBestSeen.get();
                    }
                    else
                    {
                        Ref<S> referenceToBestSeenScore = new Ref<S>(bestSeenScore);
                        Ref<N> referenceToBestSeenInternalEdgeNodeA = new Ref<N>(bestSeenInternalEdgeNodeA );
                        Ref<N> referenceToBestSeenInternalEdgeNodeB = new Ref<N>(bestSeenInternalEdgeNodeB);
                        Ref<N> referenceToBestSeenSwapNodeA = new Ref<N>(bestSeenSwapNodeA);
                        Ref<E> referenceToBestSeenSwapEdgeA = new Ref<E>(bestSeenSwapEdgeA);
                        Ref<N> referenceToBestSeenSwapNodeB = new Ref<N>(bestSeenSwapNodeB);
                        Ref<E> referenceToBestSeenSwapEdgeB = new Ref<E>(bestSeenSwapEdgeB);
                        Ref<Boolean> referenceToNewBestSeen = new Ref<Boolean>(newBestSeen);
                        swapGeneralCase(tree, getTreeScore, isBetterScore, givenTreeScore, edge, nodesOfEdge.Item1, nodesOfEdge.Item2, /*ref*/ referenceToBestSeenScore,
                                /*ref*/ referenceToBestSeenInternalEdgeNodeA, /*ref*/ referenceToBestSeenInternalEdgeNodeB, /*ref*/ referenceToBestSeenSwapNodeA, /*ref*/ referenceToBestSeenSwapEdgeA, /*ref*/ referenceToBestSeenSwapNodeB, /*ref*/ referenceToBestSeenSwapEdgeB,
                                /*ref*/ referenceToNewBestSeen);
                        bestSeenScore = referenceToBestSeenScore.get();
                        bestSeenInternalEdgeNodeA = referenceToBestSeenInternalEdgeNodeA.get();
                        bestSeenInternalEdgeNodeB = referenceToBestSeenInternalEdgeNodeB.get();
                        bestSeenSwapNodeA = referenceToBestSeenSwapNodeA.get();
                        bestSeenSwapEdgeA = referenceToBestSeenSwapEdgeA.get();
                        bestSeenSwapNodeB = referenceToBestSeenSwapNodeB.get();
                        bestSeenSwapEdgeB = referenceToBestSeenSwapEdgeB.get();
                        newBestSeen = referenceToNewBestSeen.get();
                    }
                }
            }

            if (!newBestSeen)
            {
                return tree;
            }
            else
            {
               
                _removeEdge.execute(tree, bestSeenSwapEdgeA);
                E newEdge1 = _makeEdge.execute(tree, bestSeenInternalEdgeNodeB, bestSeenSwapNodeA);
                _addEdge.execute(tree, newEdge1);

           
                _removeEdge.execute(tree, bestSeenSwapEdgeB);
                E newEdge2 = _makeEdge.execute(tree, bestSeenInternalEdgeNodeA, bestSeenSwapNodeB);
                _addEdge.execute(tree, newEdge2);

            }
            
        }

        return tree;
    }

    private void swapDoubleDegreeThree(T tree, Func1<T, S> getTreeScore, Func2<S, S, Boolean> isBetterScore, S givenTreeScore, E internalEdge,
                          N internalEdgeNodeA, N internalEdgeNodeB, Ref<S> bestSeenScore,
                          Ref<N> bestSeenInternalEdgeNodeA, Ref<N> bestSeenInternalEdgeNodeB, Ref<N> bestSeenSwapNodeA, Ref<E> bestSeenSwapEdgeA, Ref<N> bestSeenSwapNodeB, Ref<E> bestSeenSwapEdgeB,
                          Ref<Boolean> newBestSeen)
    {
        boolean treeRooted = isRooted.execute(tree);
        ArrayList<E> swapEdgesA = new ArrayList<E>();
        for(E edge : getIncidentEdges.execute(tree, internalEdgeNodeA))
        {
            if(!edge.equals(internalEdge))
            {
                if(treeRooted)
                {
                    if(!isDestinationNode.execute(tree, internalEdgeNodeA, edge))
                    {
                        swapEdgesA.add(edge);
                    }
                }
                else
                {
                    swapEdgesA.add(edge);
                }
            }
        }

        LinkedList<E> swapEdgesB = new LinkedList<E>();

        for(E edge : getIncidentEdges.execute(tree, internalEdgeNodeB))
        {
            if(!edge.equals(internalEdge))
            {
                if(treeRooted)
                {
                    if(!isDestinationNode.execute(tree, internalEdgeNodeB, edge))
                    {
                       swapEdgesB.add(edge);
                    }
                }
                else
                {
                    swapEdgesB.add(edge);
                }
            }

            if(swapEdgesA.size() == 2 && swapEdgesB.size() == 1)
            {
                break;
            }
        }

        swapHelp(tree, getTreeScore, isBetterScore, givenTreeScore, internalEdge, swapEdgesA, internalEdgeNodeA, swapEdgesB, internalEdgeNodeB, /*ref*/ bestSeenScore,
                          /*ref*/ bestSeenInternalEdgeNodeA, /*ref*/ bestSeenInternalEdgeNodeB, /*ref*/ bestSeenSwapNodeA, /*ref*/ bestSeenSwapEdgeA, /*ref*/ bestSeenSwapNodeB, /*ref*/ bestSeenSwapEdgeB,
                          /*ref*/ newBestSeen);
    }

    private void swapGeneralCase(T tree, Func1<T, S> getTreeScore, Func2<S, S, Boolean> isBetterScore, S givenTreeScore, E internalEdge,
                          N internalEdgeNodeA, N internalEdgeNodeB, Ref<S> bestSeenScore,
                          Ref<N> bestSeenInternalEdgeNodeA, Ref<N> bestSeenInternalEdgeNodeB, Ref<N> bestSeenSwapNodeA, Ref<E> bestSeenSwapEdgeA, Ref<N> bestSeenSwapNodeB, Ref<E> bestSeenSwapEdgeB,
                          Ref<Boolean> newBestSeen)
    {
         LinkedList<E> swapEdgesA = new LinkedList<E>();
         boolean treeRooted = isRooted.execute(tree);

        for(E edge : getIncidentEdges.execute(tree, internalEdgeNodeA))
        {
             if(treeRooted)
                {
                    if(!isDestinationNode.execute(tree, internalEdgeNodeA, edge))
                    {
                        swapEdgesA.add(edge);
                    }
                }
                else
                {
                    swapEdgesA.add(edge);
                }
        }

        LinkedList<E> swapEdgesB = new LinkedList<E>();

        for(E edge : getIncidentEdges.execute(tree, internalEdgeNodeB))
        {
            if(!edge.equals(internalEdge))
            {
                if(treeRooted)
                {
                    if(!isDestinationNode.execute(tree, internalEdgeNodeB, edge))
                    {
                       swapEdgesB.add(edge);
                    }
                }
                else
                {
                    swapEdgesB.add(edge);
                }
            }
        }

        swapHelp(tree, getTreeScore, isBetterScore, givenTreeScore, internalEdge, swapEdgesA, internalEdgeNodeA, swapEdgesB, internalEdgeNodeB, /*ref*/ bestSeenScore, 
                          /*ref*/ bestSeenInternalEdgeNodeA, /*ref*/ bestSeenInternalEdgeNodeB, /*ref*/ bestSeenSwapNodeA, /*ref*/ bestSeenSwapEdgeA, /*ref*/ bestSeenSwapNodeB, /*ref*/ bestSeenSwapEdgeB,
                          /*ref*/ newBestSeen);
    }

    private void swapHelp(T tree, Func1<T, S> getTreeScore, Func2<S, S, Boolean> isBetterScore, S givenTreeScore, E internalEdge, Iterable<E> swapEdgesA,
                          N internalEdgeNodeA, Iterable<E> swapEdgesB, N internalEdgeNodeB, Ref<S> bestSeenScore,
                          Ref<N> bestSeenInternalEdgeNodeA, Ref<N> bestSeenInternalEdgeNodeB, Ref<N> bestSeenSwapNodeA, Ref<E> bestSeenSwapEdgeA, Ref<N> bestSeenSwapNodeB, Ref<E> bestSeenSwapEdgeB,
                          Ref<Boolean> newBestSeen)
    {

        for (E swapA : swapEdgesA)
        {
            Tuple<N,N> nodesOfSwapA = getNodesOfEdge.execute(tree, swapA);
            N transplantInA = nodesOfSwapA.Item1.equals(internalEdgeNodeA) ? nodesOfSwapA.Item2 : nodesOfSwapA.Item1;
            _removeEdge.execute(tree, swapA);
            E newEdge1 = _makeEdge.execute(tree, internalEdgeNodeB, transplantInA);
            _addEdge.execute(tree, newEdge1);
            for (E swapB : swapEdgesB)
            {
                Tuple<N, N> nodesOfSwapB = getNodesOfEdge.execute(tree, swapB);
                N transplantInB = nodesOfSwapB.Item1.equals(internalEdgeNodeB) ? nodesOfSwapB.Item2 : nodesOfSwapB.Item1;

                _removeEdge.execute(tree, swapB);
                E newEdge2 = _makeEdge.execute(tree, internalEdgeNodeA, transplantInB);
                _addEdge.execute(tree, newEdge2);


                S newScore = getTreeScore.execute(tree);

                if (isBetterScore.execute(newScore, givenTreeScore))
                {
                    bestSeenScore.set(newScore);
                    bestSeenSwapEdgeA.set(swapA);
                    bestSeenSwapEdgeB.set(swapB);
                    bestSeenSwapNodeA.set(transplantInA);
                    bestSeenSwapNodeB.set(transplantInB);
                    bestSeenInternalEdgeNodeA.set(internalEdgeNodeA);
                    bestSeenInternalEdgeNodeB.set(internalEdgeNodeB);
                    newBestSeen.set(true);
                }

                _removeEdge.execute(tree, newEdge2);
                _addEdge.execute(tree, swapB);
            }
            _removeEdge.execute(tree, newEdge1);
            _addEdge.execute(tree, swapA);
        }
    }
}

