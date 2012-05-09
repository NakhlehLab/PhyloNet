/*
 * Copyright (c) 2012 Rice University.
 *
 * This file is part of PhyloNet.
 *
 * PhyloNet is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PhyloNet is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PhyloNet.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.rice.cs.bioinfo.programs.phylonet.algos.tbr;

import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.EdgeContractor;
import edu.rice.cs.bioinfo.programs.phylonet.algos.NodeInjector;

import java.util.LinkedList;
import java.util.List;


public class TBRSearchNeighborsExpansiveInPlace<T, N, E, S> extends TreeBisectionAndReconnectionSearcherBase<T, N, E, S>
{
    private Proc2<T, E> _removeEdge;

    private Proc2<T, N> _removeNode;

    private Proc2<T, E> _addEdge;

    private Proc2<T, N> _addNode;

    private Func3<T, N, N, E> _makeEdge;

    private Func1<T, N> _makeNewNodeInGraph;

    private Func2<T, N, Iterable<E>> _getEdgesReachableFromNode;

    private Proc1<T> BeginningToSearchNeighborhoodDelegate;

    public TBRSearchNeighborsExpansiveInPlace(Predicate1<T> isRooted, Func3<T, N, E, Boolean> isDestinationNode, Func1<T, Iterable<N>> getNodes, Func1<T, Iterable<E>> getEdges, Func2<T, N, Iterable<E>> getIncidentEdges,
                                              Func2<T, E, Tuple<N, N>> getNodesOfEdge, Proc2<T, E> removeEdge, Proc2<T, N> addNode, Proc2<T, N> removeNode, Proc2<T, E> addEdge, Func3<T, N, N, E> makeEdge, Func1<T, N> makeNewNodeInGraph,
                                              Func2<T,N,Iterable<E>> getEdgesReachableFromNode)
    {
        super(isRooted, isDestinationNode, getNodes, getEdges, getIncidentEdges, getNodesOfEdge);

        _removeEdge = removeEdge;
        _getEdgesReachableFromNode = getEdgesReachableFromNode;
        _removeNode = removeNode;
        _addEdge = addEdge;
        _makeEdge = makeEdge;
        _makeNewNodeInGraph = makeNewNodeInGraph;
        _addNode = addNode;
    }

    protected /*override*/ T searchWithoutTreeValidation(final T tree, Func1<T, S> getTreeScore, Func2<S, S, Boolean> isBetterScore, long maxIterations, S givenTreeScore)
    {
        S bestSeenScore = givenTreeScore;
        E bestSeenInternalEdge = null /* default(E)*/ ;
        N bestSeenInternalEdgeNode1 = null /* default(N)*/ , bestSeenInternalEdgeNode2 = null /* default(N)*/ ;
        E bestSeenSubtreeEdge1 = null /* default(E)*/ , bestSeenSubtreeEdge2 = null /* default(E)*/ ;
        final boolean isTreeRooted = this.isRooted.execute(tree);

        for (int i = 0; i < maxIterations; i++)
        {


            boolean newBestSeen = false;

            //  bool bestSeenSubtreeEdge1OutOfDate = false, bestSeenSubtreeEdge2OutOfDate = false;

            LinkedList<E> internalEdges = new LinkedList<E>();

            for (E edge : getEdges.execute(tree))
            {
                if (this.isInternalEdge.execute(tree, edge))
                {
                    internalEdges.add(edge);
                }
            }


            N glueNode1 = _makeNewNodeInGraph.execute(tree);
            N glueNode2 = _makeNewNodeInGraph.execute(tree);
            for (E internalEdge : (Iterable<E>) new LinkedList<E>(internalEdges)) // Linked list to avoid potential concurrent modification exception
            {

                Tuple<N, N> internalEdgeNodes = this.getNodesOfEdge.execute(tree, internalEdge);
                N internalEdgeNode1 = internalEdgeNodes.Item1;
                N internalEdgeNode2 = internalEdgeNodes.Item2;
                boolean rootInSubTree1 = isTreeRooted && isDestinationNode.execute(tree, internalEdgeNode2, internalEdge);
                _removeEdge.execute(tree, internalEdge);



                N internalEdgeNode1SmoothedNode1 = null, internalEdgeNode1SmoothedNode2 = null;
                E internalEdgeNode1SmoothedEdge1 = null, internalEdgeNode1SmoothedEdge2 = null;
                E forget = null;
                Ref<N> referenceToInternalEdgeNode1SmoothedNode1 = new Ref<N>(internalEdgeNode1SmoothedNode1);
                Ref<N> referenceToInternalEdgeNode1SmoothedNode2 = new Ref<N>(internalEdgeNode1SmoothedNode2);
                Ref<E> referenceToForget = new Ref<E>(forget);
                Ref<E> referenceToInternalEdgeNode1SmoothedEdge1 = new Ref<E>(internalEdgeNode1SmoothedEdge1);
                Ref<E> referenceToInternalEdgeNode1SmoothedEdge2 = new Ref<E>(internalEdgeNode1SmoothedEdge2);
                boolean internalEdgeNode1Smoothed = EdgeContractor.smoothInlineNodeIfNeeded(tree, isTreeRooted, this.isDestinationNode, internalEdgeNode1, getIncidentEdges, getNodesOfEdge, _removeEdge, _addEdge, _makeEdge,
                        /*out*/ referenceToInternalEdgeNode1SmoothedNode1, /*out*/ referenceToInternalEdgeNode1SmoothedNode2,
                        /*out*/ referenceToForget, /*out*/ referenceToInternalEdgeNode1SmoothedEdge1, /*out*/ referenceToInternalEdgeNode1SmoothedEdge2);
                internalEdgeNode1SmoothedNode1 = referenceToInternalEdgeNode1SmoothedNode1.get();
                internalEdgeNode1SmoothedNode2 = referenceToInternalEdgeNode1SmoothedNode2.get();
                forget = referenceToForget.get();
                internalEdgeNode1SmoothedEdge1 = referenceToInternalEdgeNode1SmoothedEdge1.get();
                internalEdgeNode1SmoothedEdge2 = referenceToInternalEdgeNode1SmoothedEdge2.get();
                if (internalEdgeNode1Smoothed)
                    _removeNode.execute(tree, internalEdgeNode1);

                N internalEdgeNode2SmoothedNode1 = null, internalEdgeNode2SmoothedNode2 = null;
                E internalEdgeNode2SmoothedEdge1 = null, internalEdgeNode2SmoothedEdge2 = null;
                Ref<N> referenceToInternalEdgeNode2SmoothedNode1 = new Ref<N>(internalEdgeNode2SmoothedNode1);
                Ref<N> referenceToInternalEdgeNode2SmoothedNode2 = new Ref<N>(internalEdgeNode2SmoothedNode2);
                referenceToForget.set(forget);
                Ref<E> referenceToInternalEdgeNode2SmoothedEdge1 = new Ref<E>(internalEdgeNode2SmoothedEdge1);
                Ref<E> referenceToInternalEdgeNode2SmoothedEdge2 = new Ref<E>(internalEdgeNode2SmoothedEdge2);
                boolean internalEdgeNode2Smoothed = EdgeContractor.smoothInlineNodeIfNeeded(tree, isTreeRooted, this.isDestinationNode, internalEdgeNode2, getIncidentEdges, getNodesOfEdge, _removeEdge, _addEdge, _makeEdge,
                        /*out*/ referenceToInternalEdgeNode2SmoothedNode1, /*out*/ referenceToInternalEdgeNode2SmoothedNode2,
                        /*out*/ referenceToForget, /*out*/ referenceToInternalEdgeNode2SmoothedEdge1, /*out*/ referenceToInternalEdgeNode2SmoothedEdge2);
                internalEdgeNode2SmoothedNode1 = referenceToInternalEdgeNode2SmoothedNode1.get();
                internalEdgeNode2SmoothedNode2 = referenceToInternalEdgeNode2SmoothedNode2.get();
                forget = referenceToForget.get();
                internalEdgeNode2SmoothedEdge1 = referenceToInternalEdgeNode2SmoothedEdge1.get();
                internalEdgeNode2SmoothedEdge2 = referenceToInternalEdgeNode2SmoothedEdge2.get();
                if (internalEdgeNode2Smoothed)
                    _removeNode.execute(tree, internalEdgeNode2);

                N subtree1Node = internalEdgeNode1Smoothed ? internalEdgeNode1SmoothedNode1 : internalEdgeNode1;
                N subtree2Node = internalEdgeNode2Smoothed ? internalEdgeNode2SmoothedNode1 : internalEdgeNode2;

                LinkedList<E> subtree1Edges = new LinkedList<E>();

                for(E edge : _getEdgesReachableFromNode.execute(tree, subtree1Node))
                {
                    subtree1Edges.add(edge);
                }

                LinkedList<E> subtree2Edges = new LinkedList<E>();

                for(E edge : _getEdgesReachableFromNode.execute(tree, subtree2Node))
                {
                    subtree2Edges.add(edge);
                }


                for (E subtree1Edge : subtree1Edges)
                {

                    NodeInjector.NodeInjectorUndoAction<T, N, E> subbtree1EdgeInjectionUndo =
                            NodeInjector.injectNodeIntoEdge(tree, isTreeRooted, isDestinationNode, subtree1Edge, glueNode1, getNodesOfEdge, _removeEdge, _addEdge, _makeEdge, true);

                    List<Tuple<E,E>> reverseUndoActionsST1 = null;
                    if(isTreeRooted && !rootInSubTree1)
                    {
                        Tuple<N,N> nodesOfEdge = getNodesOfEdge.execute(tree, subtree1Edge);
                        N sourceNode = isDestinationNode.execute(tree, nodesOfEdge.Item1, subtree1Edge) ? nodesOfEdge.Item2 : nodesOfEdge.Item1;

                        reverseUndoActionsST1 = reverseEdgePath(tree, sourceNode, glueNode1, getIncidentEdges, isDestinationNode, _removeEdge, _makeEdge);
                    }

                    for (E subtree2Edge : (Iterable<E>) subtree2Edges)
                    {


                        NodeInjector.NodeInjectorUndoAction<T, N, E> subtree2EdgeInjectionUndo =
                                NodeInjector.injectNodeIntoEdge(tree, isTreeRooted, isDestinationNode, subtree2Edge, glueNode2, getNodesOfEdge, _removeEdge, _addEdge, _makeEdge, true);

                        List<Tuple<E,E>> reverseUndoActionsST2 = null;
                        if(isTreeRooted && rootInSubTree1)
                        {
                            Tuple<N,N> nodesOfEdge = getNodesOfEdge.execute(tree, subtree2Edge);
                            N sourceNode = isDestinationNode.execute(tree, nodesOfEdge.Item1, subtree2Edge) ? nodesOfEdge.Item2 : nodesOfEdge.Item1;

                            reverseUndoActionsST2 = reverseEdgePath(tree, sourceNode, glueNode2, getIncidentEdges, isDestinationNode, _removeEdge, _makeEdge);
                        }


                        E reconnectionEdge = rootInSubTree1 ?  _makeEdge.execute(tree, glueNode1, glueNode2) : _makeEdge.execute(tree, glueNode2, glueNode1);
                        _addEdge.execute(tree, reconnectionEdge);




                        S newScore = getTreeScore.execute(tree);

                        if (isBetterScore.execute(newScore, givenTreeScore))
                        {
                            bestSeenScore = newScore;
                            bestSeenInternalEdge = internalEdge;
                            bestSeenInternalEdgeNode1 = internalEdgeNode1;
                            bestSeenInternalEdgeNode2 = internalEdgeNode2;
                            bestSeenSubtreeEdge1 = subtree1Edge;
                            bestSeenSubtreeEdge2 = subtree2Edge;
                            newBestSeen = true;

                        }

                        _removeEdge.execute(tree, reconnectionEdge);

                        if(reverseUndoActionsST2 != null)
                        {
                            for(Tuple<E,E> undoAction : reverseUndoActionsST2)
                            {
                                _addEdge.execute(tree, undoAction.Item1);
                                _removeEdge.execute(tree, undoAction.Item2);
                            }
                        }

                        subtree2EdgeInjectionUndo.undoInjection();

                    }

                    if(reverseUndoActionsST1 != null)
                    {
                        for(Tuple<E,E> undoAction : reverseUndoActionsST1)
                        {
                            _addEdge.execute(tree, undoAction.Item1);
                            _removeEdge.execute(tree, undoAction.Item2);
                        }
                    }

                    subbtree1EdgeInjectionUndo.undoInjection();


                }


                if (internalEdgeNode1Smoothed)
                {
                    _addNode.execute(tree, internalEdgeNode1);
                    EdgeContractor.undoSmoothing(tree, internalEdgeNode1, getIncidentEdges, getNodesOfEdge, _removeEdge, _addEdge, _makeEdge, internalEdgeNode1SmoothedNode1, internalEdgeNode1SmoothedNode2, internalEdgeNode1SmoothedEdge1, internalEdgeNode1SmoothedEdge2);
                }

                if (internalEdgeNode2Smoothed)
                {
                    _addNode.execute(tree, internalEdgeNode2);
                    EdgeContractor.undoSmoothing(tree, internalEdgeNode2, getIncidentEdges, getNodesOfEdge, _removeEdge, _addEdge, _makeEdge, internalEdgeNode2SmoothedNode1, internalEdgeNode2SmoothedNode2, internalEdgeNode2SmoothedEdge1, internalEdgeNode2SmoothedEdge2);
                }

                this._addEdge.execute(tree, internalEdge);
            }
            _removeNode.execute(tree, glueNode1);
            _removeNode.execute(tree, glueNode2);

            if (!newBestSeen)
            {
                return tree;
            }
            else
            {
                _removeEdge.execute(tree, bestSeenInternalEdge);
                _addNode.execute(tree, glueNode1);
                _addNode.execute(tree, glueNode2);

                Func3<N, Tuple<N, N>, E, E> smoothAndRemoveInlineNodeIfNeeded = new Func3<N, Tuple<N, N>, E, E>() {
                    public E execute(N bestSeenInternalEdgeNode, Tuple<N,N> bestSeenSubtreeNodesFromInternalEdgeNodeSubTree, E bestSeenSubtreeEdge) {

                        N smoothNode1 = null, smoothNode2 = null;
                        E smoothEdge = null, forget1 = null, forget2 = null;
                        Ref<N> referenceToSmoothNode1 = new Ref<N>(smoothNode1);
                        Ref<N> referenceToSmoothNode2 = new Ref<N>(smoothNode2);
                        Ref<E> referenceToSmoothEdge = new Ref<E>(smoothEdge);
                        Ref<E> referenceToForget1 = new Ref<E>(forget1);
                        Ref<E> referenceToForget2 = new Ref<E>(forget2);
                        boolean smoothed = EdgeContractor.smoothInlineNodeIfNeeded(tree, isTreeRooted, isDestinationNode, bestSeenInternalEdgeNode, getIncidentEdges, getNodesOfEdge, _removeEdge, _addEdge, _makeEdge,
                                /*out*/ referenceToSmoothNode1, /*out*/ referenceToSmoothNode2, /*out*/ referenceToSmoothEdge, /*out*/ referenceToForget1, /*out*/ referenceToForget2);
                        smoothNode1 = referenceToSmoothNode1.get();
                        smoothNode2 = referenceToSmoothNode2.get();
                        smoothEdge = referenceToSmoothEdge.get();
                        forget1 = referenceToForget1.get();
                        forget2 = referenceToForget2.get();
                        if (smoothed)
                        {
                            _removeNode.execute(tree, bestSeenInternalEdgeNode);

                            if ((smoothNode1.equals(bestSeenSubtreeNodesFromInternalEdgeNodeSubTree.Item1) && smoothNode2.equals(bestSeenSubtreeNodesFromInternalEdgeNodeSubTree.Item2) ||
                                    (smoothNode2.equals(bestSeenSubtreeNodesFromInternalEdgeNodeSubTree.Item1) && smoothNode1.equals(bestSeenSubtreeNodesFromInternalEdgeNodeSubTree.Item2))))
                                return smoothEdge;
                        }

                        return bestSeenSubtreeEdge;
                    }
                };



                Tuple<N, N> bestSeenSubtreeEdge1Nodes = getNodesOfEdge.execute(tree, bestSeenSubtreeEdge1);
                Tuple<N, N> bestSeenSubtreeEdge2Nodes = getNodesOfEdge.execute(tree, bestSeenSubtreeEdge2);



                bestSeenSubtreeEdge1 = smoothAndRemoveInlineNodeIfNeeded.execute(bestSeenInternalEdgeNode1, bestSeenSubtreeEdge1Nodes, bestSeenSubtreeEdge1);
                bestSeenSubtreeEdge2 = smoothAndRemoveInlineNodeIfNeeded.execute(bestSeenInternalEdgeNode2, bestSeenSubtreeEdge2Nodes, bestSeenSubtreeEdge2);



                NodeInjector.injectNodeIntoEdge(tree, isTreeRooted, isDestinationNode, bestSeenSubtreeEdge1, glueNode1, getNodesOfEdge, _removeEdge, _addEdge, _makeEdge, false);
                NodeInjector.injectNodeIntoEdge(tree, isTreeRooted, isDestinationNode, bestSeenSubtreeEdge2, glueNode2, getNodesOfEdge, _removeEdge, _addEdge, _makeEdge, false);
                _addEdge.execute(tree, _makeEdge.execute(tree, glueNode1, glueNode2));
                // return SearchWithoutTreeValidation(tree, getTreeScore, isBetterScore, maxIterations - 1, bestSeenScore);
            }
        }

        return tree;
    }

    private List<Tuple<E,E>> reverseEdgePath(T tree, N sourceNode, N reversePathParent, Func2<T, N, Iterable<E>> getIncidentEdges, Func3<T, N, E, Boolean> destinationNode, Proc2<T, E> removeEdge, Func3<T, N, N, E> makeEdge) {

        LinkedList<Tuple<E,E>> accum = new LinkedList<Tuple<E, E>>();

        for(E incidentEdge : getIncidentEdges.execute(tree, sourceNode))
        {
            Tuple<N,N> nodesOfEdge = getNodesOfEdge.execute(tree, incidentEdge);
            N notSourceNode = nodesOfEdge.Item1 == sourceNode ? nodesOfEdge.Item2 : nodesOfEdge.Item1;

            if(reversePathParent != notSourceNode)
            {
                if(!isDestinationNode.execute(tree, notSourceNode, incidentEdge))
                {
                    removeEdge.execute(tree, incidentEdge);
                    E edge = makeEdge.execute(tree, sourceNode, notSourceNode);
                    _addEdge.execute(tree, edge);
                    accum.add(new Tuple<E, E>(incidentEdge, edge));
                }
            }
        }

        return accum;
    }


}
