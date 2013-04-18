//package edu.rice.cs.bioinfo.programs.soranus.models.analysis;
//
//import java.util.Set;
//import java.util.Stack;
//
///**
// * Created with IntelliJ IDEA.
// * User: Matt
// * Date: 4/4/13
// * Time: 3:15 PM
// * To change this template use File | Settings | File Templates.
// */
//public abstract class MinSpanTreeSnpInferrerTemplate<S,E,D,G> implements MinSpanTreeSnpInferrer<S,E>
//{
//
//    public Set<Set<E>> inferMinTrees(Set<S> sequencings)
//    {
//        G completeGraph = makeEmptyGraph();
//        for(S sequencing : sequencings)
//        {
//            addSequencingNode(completeGraph, sequencing);
//        }
//
//        Stack<S> toLink = new Stack<S>();
//        toLink.addAll(sequencings);
//
//        while(!toLink.isEmpty())
//        {
//            S node1 = toLink.pop();
//            for(S node2 : toLink)
//            {
//                D snpCount = getSnpCount(node1, node2);
//                E edge = makeEdge(node1, node2, snpCount);
//                addEdge(edge, completeGraph);
//            }
//        }
//
//        return inferMinTrees(completeGraph);
//    }
//
//    protected abstract Set<Set<E>> inferMinTrees(G completeGraph);
//
//    protected abstract E makeEdge(S node1, S node2, D snpCount);
//
//    protected abstract void addEdge(E edge, G graph);
//
//    protected abstract D getSnpCount(S sequencing1, S sequencing2);
//
//    protected abstract D getSnpCount(E edge);
//
//    protected abstract void addSequencingNode(G graph, S sequencing);
//
//    protected abstract G makeEmptyGraph();
//
//
//}
