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

package edu.rice.cs.bioinfo.programs.phylonet.algos;

import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.Func3;
import edu.rice.cs.bioinfo.library.programming.Proc2;
import edu.rice.cs.bioinfo.library.programming.Tuple;

public class NodeInjector
    {
        public static class NodeInjectorUndoAction<T, N, E>
        {
             E __removedEdge;private E get_removedEdge(){ return __removedEdge; }; private void set_removedEdge(E value){__removedEdge = value; };

             E __addedEdge1;private E get_addedEdge1(){ return __addedEdge1; }; private void set_addedEdge1(E value){__addedEdge1 = value; };

             E __addedEdge2;private E get_addedEdge2(){ return __addedEdge2; }; private void set_addedEdge2(E value){__addedEdge2 = value; };

             T __tree;private T get_tree(){ return __tree; }; private void set_tree(T value){__tree = value; };

            private Proc2<T, E> _removeEdge;

            private Proc2<T, E> _addEdge;

            NodeInjectorUndoAction(E removedEdge, E addedEdge1, E addedEdge2, T tree, Proc2<T, E> removeEdge, Proc2<T, E> addEdge)
            {
                set_removedEdge(removedEdge);
                set_addedEdge1(addedEdge1);
                set_addedEdge2(addedEdge2);
                set_tree(tree);
                _addEdge = addEdge;
                _removeEdge = removeEdge;
            }

            public void undoInjection()
            {
                _removeEdge.execute(get_tree(), get_addedEdge1());
                _removeEdge.execute(get_tree(), get_addedEdge2());
                _addEdge.execute(get_tree(), get_removedEdge());
            }


        }

        public static <T, N, E> NodeInjectorUndoAction<T,N,E> injectNodeIntoEdge(T tree, boolean isRooted, Func3<T, N, E, Boolean> isDestinationNode, E edge, N node,
                                                                                 Func2<T, E, Tuple<N, N>> getNodesOfEdge, Proc2<T, E> removeEdge, Proc2<T, E> addEdge,
                                                                                 Func3<T, N, N, E> makeEdge, boolean makeNodeOnlySource)
        {

            Tuple<N, N> edgeNodes = getNodesOfEdge.execute(tree, edge);
            boolean item1isSource = isRooted && isDestinationNode.execute(tree, edgeNodes.Item2, edge);
            removeEdge.execute(tree, edge);
            E edgeToAdd1 = null;
            E edgeToAdd2 = null;

            if(isRooted && makeNodeOnlySource)
            {
                edgeToAdd1 = makeEdge.execute(tree, node, edgeNodes.Item1);
                edgeToAdd2 = makeEdge.execute(tree, node, edgeNodes.Item2);
            }
            else if(item1isSource)
            {
                edgeToAdd1 = makeEdge.execute(tree, edgeNodes.Item1, node);
                edgeToAdd2 = makeEdge.execute(tree, node, edgeNodes.Item2);
            }
            else
            {
                edgeToAdd1 = makeEdge.execute(tree, edgeNodes.Item2, node);
                edgeToAdd2 = makeEdge.execute(tree, node, edgeNodes.Item1);
            }

            addEdge.execute(tree, edgeToAdd1);
            addEdge.execute(tree, edgeToAdd2);

            return new NodeInjectorUndoAction<T, N, E>(edge, edgeToAdd1, edgeToAdd2, tree, removeEdge, addEdge);
        }
    }
