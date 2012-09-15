///*
// * Copyright (c) 2012 Rice University.
// *
// * This file is part of PhyloNet.
// *
// * PhyloNet is free software: you can redistribute it and/or modify
// * it under the terms of the GNU General Public License as published by
// * the Free Software Foundation, either version 3 of the License, or
// * (at your option) any later version.
// *
// * PhyloNet is distributed in the hope that it will be useful,
// * but WITHOUT ANY WARRANTY; without even the implied warranty of
// * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// * GNU General Public License for more details.
// *
// * You should have received a copy of the GNU General Public License
// * along with PhyloNet.  If not, see <http://www.gnu.org/licenses/>.
// */
//
//package edu.rice.cs.bioinfo.programs.phylonet.structs.network;
//
//import com.sun.org.apache.xalan.internal.xsltc.runtime.Node;
//import edu.rice.cs.bioinfo.library.phylogenetics.*;
//import edu.rice.cs.bioinfo.library.programming.Tuple;
//import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
//
//import java.util.HashSet;
//import java.util.LinkedList;
//import java.util.List;
//
///**
// * Created by IntelliJ IDEA.
// * User: Matt
// * Date: 5/15/12
// * Time: 1:14 PM
// * To change this template use File | Settings | File Templates.
// */
//class NetworkToGraphAdapter<G extends Network<T>,T> implements Graph<NetNode<T>, Tuple<NetNode<T>, NetNode<T>>>
//{
//    private LinkedList<NetNode<T>> _disconnectedNodes = new LinkedList<NetNode<T>>();
//
//    public final G Network;
//
//    public final double DefaultBranchLength;
//
//    public NetworkToGraphAdapter(G network, double defaultBranchLength)
//    {
//        Network = network;
//        DefaultBranchLength = defaultBranchLength;
//    }
//
//    public void addNode(NetNode<T> node) {
//
//        for(NetNode<T> netNode : Network.dfs())
//        {
//            if(node.equals(netNode))
//            {
//                return;
//            }
//        }
//
//        _disconnectedNodes.add(node);
//
//    }
//
//    public void removeNode(NetNode<T> node)
//    {
//        for(NetNode<T> parent : node.getParents())
//        {
//            parent.removeChild(node);
//        }
//        _disconnectedNodes.remove(node);
//    }
//
//    public void addEdge(Tuple<NetNode<T>, NetNode<T>> edge)
//    {
//        for(NetNode netNode : Network.dfs())
//        {
//            if(edge.Item1.equals(netNode))
//            {
//                netNode.adoptChild(edge.Item2, DefaultBranchLength);
//                return;
//            }
//        }
//
//        throw new IllegalArgumentException("Given source node is not a member of the inner network.  Node: " + edge.Item1);
//    }
//
//    public void removeEdge(Tuple<NetNode<T>, NetNode<T>> edge)
//    {
//        boolean isChild = false;
//
//        for(NetNode<T> child : edge.Item1.getChildren())
//        {
//            if(edge.Item2.equals(child))
//            {
//                isChild = true;
//                break;
//            }
//        }
//
//        if(isChild)
//        {
//            edge.Item1.removeChild(edge.Item2);
//            _disconnectedNodes.add(edge.Item2);
//        }
//        else
//        {
//            throw new IllegalArgumentException("No such edge in network.");
//        }
//    }
//
//    public Iterable<NetNode<T>> getNodes() {
//
//        List<NetNode<T>> nodes = IterableHelp.toList(Network.dfs());
//        nodes.addAll(_disconnectedNodes);
//        return nodes;
//    }
//
//    public Iterable<Tuple<NetNode<T>, NetNode<T>>> getEdges()
//    {
//        LinkedList<Tuple<NetNode<T>, NetNode<T>>> edges = new LinkedList<Tuple<NetNode<T>, NetNode<T>>>();
//
//        HashSet<NetNode<T>> seenNodes = new HashSet<NetNode<T>>();
//        for(NetNode<T> node : Network.dfs())
//        {
//            if(!seenNodes.contains(node))
//            {
//                seenNodes.add(node);
//
//                for(NetNode<T> directDesc : node.getChildren())
//                {
//                    edges.add(new Tuple<NetNode<T>, NetNode<T>>(node, directDesc));
//                }
//            }
//        }
//
//        return edges;
//    }
//
//    public Tuple<NetNode<T>, NetNode<T>> getNodesOfEdge(Tuple<NetNode<T>, NetNode<T>> edge)
//    {
//        return edge;  //To change body of implemented methods use File | Settings | File Templates.
//    }
//
//    public Iterable<Tuple<NetNode<T>, NetNode<T>>> getIncidentEdges(NetNode<T> node)
//    {
//        LinkedList<Tuple<NetNode<T>, NetNode<T>>> edges = new LinkedList<Tuple<NetNode<T>, NetNode<T>>>();
//
//        for(NetNode<T> directPred : node.getParents())
//        {
//            edges.add(new Tuple<NetNode<T>, NetNode<T>>(directPred, node));
//        }
//
//        for(NetNode<T> directSuc : node.getChildren())
//        {
//            edges.add(new Tuple<NetNode<T>, NetNode<T>>(node, directSuc));
//        }
//
//        return edges;
//    }
//
//    public Tuple<NetNode<T>, NetNode<T>> getEdge(NetNode<T> source, NetNode<T> destination) {
//
//        if(IterableHelp.toList(source.getChildren()).contains(destination))
//        {
//            return new Tuple<NetNode<T>, NetNode<T>>(source, destination);
//        }
//        else
//        {
//            throw new IllegalArgumentException("No such edge.");
//        }
//    }
//
//    public boolean isRooted() {
//        return true;  //To change body of implemented methods use File | Settings | File Templates.
//    }
//}
