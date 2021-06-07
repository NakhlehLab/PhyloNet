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

package edu.rice.cs.bioinfo.programs.phylonet.algos.spr;

import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.EdgeContractor;
import edu.rice.cs.bioinfo.programs.phylonet.algos.NodeInjector;

import java.util.LinkedList;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 3/23/12
 * Time: 3:26 PM
 * To change this template use File | Settings | File Templates.
 */
public class SPRSearcherNeighborsExpansiveInPlace<T, N, E, S> extends SubtreePruningAndRegraftingSearcherBase<T, N, E, S>
{
    private /*final*/ Proc2<T, E> _removeEdge;

    private /*final*/ Proc2<T, E> _addEdge;

    private /*final*/ Func3<T, N, N, E> _makeEdge;

    private /*final*/ Proc2<T, N> _removeNode;

    private /*final*/ Proc2<T, N> _addNode;

    private /*final*/ Func2<T, N, Iterable<E>> _getEdgesReachableFromNode;

    private /*final*/ Func1<T, N> _makeNewNodeInGraph;

    private Proc1<T> BeginningToSearchNeighborhoodDelegate;

    public SPRSearcherNeighborsExpansiveInPlace(Predicate1<T> isRooted, Func3<T, N, E, Boolean> isDestinationNode, Func1<T, Iterable<N>> getNodes, Func1<T, Iterable<E>> getEdges, Func2<T, N, Iterable<E>> getIncidentEdges,
                                       Func2<T, E, Tuple<N, N>> getNodesOfEdge, Proc2<T, E> removeEdge, Proc2<T, N> addNode, Proc2<T, N> removeNode, Proc2<T, E> addEdge, Func3<T, N, N, E> makeEdge, Func1<T, N> makeNewNodeInGraph,
                                       Func2<T,N,Iterable<E>> getEdgesReachableFromNode)
    {
    	super(isRooted, isDestinationNode, getNodes, getEdges, getIncidentEdges, getNodesOfEdge);

        _removeEdge = removeEdge;
        _getEdgesReachableFromNode = getEdgesReachableFromNode;
        _removeNode = removeNode;
        _addEdge = addEdge;
        _makeEdge = makeEdge;
        _addNode = addNode;
        _makeNewNodeInGraph = makeNewNodeInGraph;
    }

    protected /*override*/ T searchWithoutTreeValidation(T tree, Func1<T, S> getTreeScore, Func2<S, S, Boolean> isBetterScore, long maxIterations, S givenTreeScore)
    {
        S bestSeenScore = givenTreeScore;
        E bestSeenRemoveEdge = null /* default(E)*/ ;
        N bestSeenPlantNode = null /* default(N)*/ ;
        N bestSeenMigrantNode = null /* default(N)*/ ;
        E bestSeenPlantEdge = null /* default(E)*/ ;
        boolean isTreeRooted = isRooted.execute(tree);

        for(int i = 0; i<maxIterations; i++)
        {
            boolean newBestSeen = false;

            LinkedList<E> edges = new LinkedList<E>();   // Linked list to avoid potential concurrent modification exception

            for(E edge : this.getEdges.execute(tree))
            {
                edges.add(edge);
            }

            for (E removeEdge : edges)
            {
                Tuple<N, N> removeEdgeNodes = this.getNodesOfEdge.execute(tree, removeEdge);
                N removeEdgeNode1 = removeEdgeNodes.Item1;
                N removeEdgeNode2 = removeEdgeNodes.Item2;
                _removeEdge.execute(tree, removeEdge);

                Ref<S> referenceToBestSeenScore = new Ref<S>(bestSeenScore);
                Ref<Boolean> referenceToNewBestSeen = new Ref<Boolean>(newBestSeen);
                Ref<E> referenceToBestSeenRemoveEdge = new Ref<E>(bestSeenRemoveEdge);
                Ref<E> referenceToBestSeenPlantEdge = new Ref<E>( bestSeenPlantEdge);
                Ref<N> referenceToBestSeenPlantNode = new Ref<N>( bestSeenPlantNode);
                Ref<N> referenceToBestSeenMigrantNode = new Ref<N>( bestSeenMigrantNode);

                if(isRooted.execute(tree))
                {
                   if(isDestinationNode.execute(tree, removeEdgeNode1, removeEdge))
                   {
                        examineAllPlants(tree, true, getTreeScore, isBetterScore, givenTreeScore, removeEdgeNode1, removeEdgeNode2, removeEdge, /*ref*/ referenceToBestSeenScore, /*ref*/ referenceToNewBestSeen, /*ref*/ referenceToBestSeenRemoveEdge,
                                 /*ref*/ referenceToBestSeenPlantEdge, /*ref*/ referenceToBestSeenPlantNode, /*ref*/ referenceToBestSeenMigrantNode);
                   }
                   else
                   {
                        examineAllPlants(tree, true, getTreeScore, isBetterScore, givenTreeScore, removeEdgeNode2, removeEdgeNode1, removeEdge, /*ref*/ referenceToBestSeenScore, /*ref*/ referenceToNewBestSeen, /*ref*/ referenceToBestSeenRemoveEdge,
                                 /*ref*/ referenceToBestSeenPlantEdge, /*ref*/ referenceToBestSeenPlantNode, /*ref*/ referenceToBestSeenMigrantNode);
                   }
                }
                else
                {

                     examineAllPlants(tree, false, getTreeScore, isBetterScore, givenTreeScore, removeEdgeNode1, removeEdgeNode2, removeEdge, /*ref*/ referenceToBestSeenScore, /*ref*/ referenceToNewBestSeen, /*ref*/ referenceToBestSeenRemoveEdge,
                                 /*ref*/ referenceToBestSeenPlantEdge, /*ref*/ referenceToBestSeenPlantNode, /*ref*/ referenceToBestSeenMigrantNode);


                     examineAllPlants(tree, false, getTreeScore, isBetterScore, givenTreeScore, removeEdgeNode2, removeEdgeNode1, removeEdge, /*ref*/ referenceToBestSeenScore, /*ref*/ referenceToNewBestSeen, /*ref*/ referenceToBestSeenRemoveEdge,
                                 /*ref*/ referenceToBestSeenPlantEdge, /*ref*/ referenceToBestSeenPlantNode, /*ref*/ referenceToBestSeenMigrantNode);
                }
                bestSeenScore = referenceToBestSeenScore.get();
                newBestSeen = referenceToNewBestSeen.get();
                bestSeenRemoveEdge = referenceToBestSeenRemoveEdge.get();
                bestSeenPlantEdge = referenceToBestSeenPlantEdge.get();
                bestSeenPlantNode = referenceToBestSeenPlantNode.get();
                bestSeenMigrantNode = referenceToBestSeenMigrantNode.get();

                _addEdge.execute(tree, removeEdge);

            }



            if (!newBestSeen)
            {
                return tree;
            }
            else
            {
                _removeEdge.execute(tree, bestSeenRemoveEdge);
                N forget1 = null;
                E forget2 = null;
                Ref<N> referenceToForget1 = new Ref<N>(forget1);
                Ref<E> referenceToForget2 = new Ref<E>(forget2);
                boolean smoothPerformed = EdgeContractor.smoothInlineNodeIfNeeded(tree, isTreeRooted, this.isDestinationNode, bestSeenPlantNode, this.getIncidentEdges, this.getNodesOfEdge, _removeEdge, _addEdge, _makeEdge,
                        /*out*/ referenceToForget1, /*out*/ referenceToForget1, /*out*/ referenceToForget2, /*out*/ referenceToForget2, /*out*/ referenceToForget2);

                if (smoothPerformed)
                    _removeNode.execute(tree, bestSeenPlantNode);

                N glueNode = _makeNewNodeInGraph.execute(tree);
                NodeInjector.NodeInjectorUndoAction<T, N, E> forget = null;
                Ref<NodeInjector.NodeInjectorUndoAction<T, N, E>> referenceToForget = new Ref<NodeInjector.NodeInjectorUndoAction<T, N, E>>(forget);
                plantAtEdge(tree, isTreeRooted, bestSeenMigrantNode, bestSeenPlantEdge, glueNode, /*out*/ referenceToForget);
            }
        }
        return tree;
    }

    private void examineAllPlants(T tree, boolean isTreeRooted, Func1<T, S> getTreeScore, Func2<S, S, Boolean> isBetterScore, S givenTreeScore, N migrantTreeNode, N plantTreeNode, E removedEdge,
                                  /*ref*/ Ref<S> bestSeenScore, /*ref*/ Ref<Boolean> newBestSeen, /*ref*/ Ref<E> bestSeenRemoveEdge, /*ref*/ Ref<E> bestSeenPlantEdge, /*ref*/ Ref<N> bestSeenPlantNode, /*ref*/ Ref<N> bestSeenMigrantNode)
    {
        N smoothedNode1 = null, smoothedNode2 = null;
        E illegalReplantEdge = null, smoothedEdge1 = null, smoothedEdge2 = null;
        Ref<N> referenceToSmoothedNode1 = new Ref<N>(smoothedNode1);
        Ref<N> referenceToSmoothedNode2 = new Ref<N>(smoothedNode2);
        Ref<E> referenceToIllegalReplantEdge = new Ref<E>(illegalReplantEdge);
        Ref<E> referenceToSmoothedEdge1 = new Ref<E>(smoothedEdge1);
        Ref<E> referenceToSmoothedEdge2 = new Ref<E>(smoothedEdge2);
        boolean smoothPerformed = EdgeContractor.smoothInlineNodeIfNeeded(tree, isTreeRooted, this.isDestinationNode, plantTreeNode, this.getIncidentEdges, this.getNodesOfEdge, _removeEdge, _addEdge, _makeEdge,
                                                                       /*out*/ referenceToSmoothedNode1, /*out*/ referenceToSmoothedNode2, /*out*/ referenceToIllegalReplantEdge, /*out*/ referenceToSmoothedEdge1, /*out*/ referenceToSmoothedEdge2);
        smoothedNode1 = referenceToSmoothedNode1.get();
        smoothedNode2 = referenceToSmoothedNode2.get();
        illegalReplantEdge = referenceToIllegalReplantEdge.get();
        smoothedEdge1 = referenceToSmoothedEdge1.get();
        smoothedEdge2 = referenceToSmoothedEdge2.get();

        if (smoothPerformed)
            _removeNode.execute(tree, plantTreeNode);

        LinkedList<E> plantEdgeCandidates = new LinkedList<E>();


        if(smoothPerformed)
        {
            for(E edge : _getEdgesReachableFromNode.execute(tree, smoothedNode1))
            {
                if(!illegalReplantEdge.equals(edge))
                {
                    plantEdgeCandidates.add(edge);
                }

            }
        }
        else
        {
            for(E edge : _getEdgesReachableFromNode.execute(tree, plantTreeNode))
            {
                    plantEdgeCandidates.add(edge);
            }
        }



        N glueNode = _makeNewNodeInGraph.execute(tree);
        for (E plantEdge : (Iterable<E>) plantEdgeCandidates)
        {
            NodeInjector.NodeInjectorUndoAction<T, N, E> injectionUndoAction = null;
            Ref<NodeInjector.NodeInjectorUndoAction<T, N, E>> referenceToInjectionUndoAction = new Ref<NodeInjector.NodeInjectorUndoAction<T, N, E>>(injectionUndoAction);
            E regraftEdge = plantAtEdge(tree, isTreeRooted, migrantTreeNode, plantEdge, glueNode, /*out*/ referenceToInjectionUndoAction);
            injectionUndoAction = referenceToInjectionUndoAction.get();

           S newScore = getTreeScore.execute(tree);

            if (isBetterScore.execute(newScore, givenTreeScore))
            {
                bestSeenScore.set(newScore);
                bestSeenRemoveEdge.set(removedEdge);
                bestSeenPlantNode.set(plantTreeNode);
                bestSeenPlantEdge.set(plantEdge);
                bestSeenMigrantNode.set(migrantTreeNode);
                newBestSeen.set(true);
            }

            undoPlantAtEdge(tree, glueNode, regraftEdge, injectionUndoAction);

        }
        _removeNode.execute(tree, glueNode);


        if (smoothPerformed)
        {
            _addNode.execute(tree, plantTreeNode);
            EdgeContractor.undoSmoothing(tree, plantTreeNode, this.getIncidentEdges, getNodesOfEdge, _removeEdge, _addEdge, _makeEdge, smoothedNode1, smoothedNode2, smoothedEdge1, smoothedEdge2);
        }


    }

    private E plantAtEdge(T tree, boolean isRooted, N migrantTreeNode, E plantEdge, N glueNode, /*out*/ Ref<NodeInjector.NodeInjectorUndoAction<T,N,E>> injectionUndoAction)
    {
        injectionUndoAction.set(NodeInjector.injectNodeIntoEdge(tree, isRooted, isDestinationNode, plantEdge, glueNode, getNodesOfEdge, _removeEdge, _addEdge, _makeEdge, false));
        E regraftEdge = _makeEdge.execute(tree, glueNode, migrantTreeNode);
        _addEdge.execute(tree, regraftEdge);

        return regraftEdge;
    }

    private void undoPlantAtEdge(T tree, N glueNode, E regraftEdge, NodeInjector.NodeInjectorUndoAction<T, N, E> injectionUndoAction)
    {
        _removeEdge.execute(tree, regraftEdge);
        injectionUndoAction.undoInjection();
       /* N forget1, forget2;
            E forget3, forget4, forget5;
            EdgeContractor.SmoothInlineNodeIfNeeded(tree, glueNode, this.getIncidentEdges, this.getNodesOfEdge, _removeEdge, _addEdge, _makeEdge, out forget1, out forget2, out forget3, out forget4, out forget5); */

    }



}
