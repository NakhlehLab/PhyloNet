package edu.rice.cs.bioinfo.library.phylogenetics;


import edu.rice.cs.bioinfo.library.programming.Func;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.Func5;

import java.util.LinkedList;
import java.util.Queue;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/25/12
 * Time: 1:11 PM
 * To change this template use File | Settings | File Templates.
 */
public class MULTreeFactory<N,NE,TE>
{
    public static class MULTreeNode<N>
    {
        public final N Content;

        public MULTreeNode(N content)
        {
            Content = content;
        }
    }


    private final Func2<GraphReadOnly<N,NE>,N,Integer> _getInDegreeStrategy = new GetInDegree<N, NE>();

    private final Func2<GraphReadOnly<N,NE>,N,Integer> _getOutDegreeStrategy = new GetOutDegree<N, NE>();

    private final Func1<GraphReadOnly<N,?>, N> _getRootStrategy = new FindRoot<N>();

    private final Func2<GraphReadOnly<N,NE>,N,Iterable<N>> _getDirectSuccessorsStrategy = new GetDirectSuccessors<N, NE>();

    Graph<MULTreeNode<N>,TE> makeMULTree(GraphReadOnly<N,NE> network, Func<Graph<MULTreeNode<N>,TE>> makeEmptyGraph,
                                        Func5<GraphReadOnly<N,NE>, NE, Graph<MULTreeNode<N>,TE>, MULTreeNode<N>, MULTreeNode<N>,TE> makeEdge)
    {

         for(N node : network.getNodes())
         {
             int inDegree = _getInDegreeStrategy.execute(network, node);
             int outDegree = _getOutDegreeStrategy.execute(network, node);

             if(inDegree == 1 && outDegree == 1)
             {
                 throw new IllegalArgumentException("Given network has a node with indegree and outderee = 1 (" + node.toString() + ").");
             }
         }

        N netRoot = _getRootStrategy.execute(network);
        MULTreeNode<N> mulTreeRoot = new MULTreeNode<N>(netRoot);

        Graph<MULTreeNode<N>,TE> mulTree = makeEmptyGraph.execute();
        mulTree.addNode(mulTreeRoot);

   //     ((STINode<Double>)(_mulTree.getRoot())).setData(1.0);
        Queue<N> source = new LinkedList<N>();
        Queue<MULTreeNode<N>> dest = new LinkedList<MULTreeNode<N>>();
        source.offer(netRoot);
        dest.offer(mulTreeRoot);
  //      long nameid = System.currentTimeMillis();
        //long nameid = 0;
        while(!source.isEmpty()){
            N parent = source.poll();
            MULTreeNode<N> peer = dest.poll();
            for (N child : _getDirectSuccessorsStrategy.execute(network, parent)) {

                NE parentToChildEdge = network.getEdge(parent, child);
                /*
                if (child.getName().equals(NetNode.NO_NAME)) {
                    child.setName("hnode" + (nameid++));
                }
                String name = child.getName();
                if(child.isNetworkNode()){
                    name = child.getName()+"TO"+parent.getName();
                }
                Integer amount = _nname2tamount.get(name);
                if(amount==null){
                    amount = 0;
                }
                _nname2tamount.put(name, ++amount);
                String newname = name + "_" + amount;  */
                MULTreeNode<N> copy = new MULTreeNode<N>(child);
                mulTree.addNode(copy);
                TE treeEdge = makeEdge.execute(network, parentToChildEdge, mulTree, peer, copy);
                mulTree.addEdge(treeEdge);
               /* if(child.isLeaf()) {
                    _stTaxa.add(newname);
                }

                // Update the distance and data for this child.
                double distance = child.getParentDistance(parent);
                if (distance == NetNode.NO_DISTANCE) {
                    //copy.setParentDistance(TNode.NO_DISTANCE);
                    copy.setParentDistance(0);
                }
                else {
                    copy.setParentDistance(distance);
                }

                double gamma = child.getParentProbability(parent);
                ((STINode<Double>)copy).setData(gamma);     */

                // Continue to iterate over the children of nn and tn.
                source.offer(child);
                dest.offer(copy);
                //index ++;
            }
        }

        return mulTree;

    }

}
