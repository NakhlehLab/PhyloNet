//package edu.rice.cs.bioinfo.programs.soranus.models.analysis;
//
//import edu.rice.cs.bioinfo.library.programming.Tuple;
//import edu.uci.ics.jung.algorithms.shortestpath.PrimMinimumSpanningTree;
//import edu.uci.ics.jung.graph.Graph;
//import edu.uci.ics.jung.graph.SparseGraph;
//import org.apache.commons.collections15.Factory;
//import org.apache.commons.collections15.Transformer;
//
//import java.util.HashSet;
//import java.util.Set;
//
///**
// * Created with IntelliJ IDEA.
// * User: Matt
// * Date: 4/4/13
// * Time: 5:24 PM
// * To change this template use File | Settings | File Templates.
// */
//public abstract class MinSpanTreeSnpInferrerJungBase<S,D>
//        extends MinSpanTreeSnpInferrerTemplate<S,MinSpanTreeSnpInferrerJungBase.DistanceEdge<S,D>,D,
//        Graph<S,MinSpanTreeSnpInferrerJungBase.DistanceEdge<S,D>>>
//{
//    public static class DistanceEdge<S,D> extends Tuple<S,S>
//    {
//        private D _distance;
//
//        public D getDistance()
//        {
//            return _distance;
//        }
//
//        public DistanceEdge(S node1, S node2, D distance)
//        {
//            super(node1, node2);
//            _distance = distance;
//        }
//    }
//
//    @Override
//    protected void addEdge(DistanceEdge<S,D> edge, Graph<S,DistanceEdge<S,D>> graph)
//    {
//        graph.addEdge(edge, edge.Item1, edge.Item2);
//    }
//
//    @Override
//    protected void addSequencingNode(Graph<S,DistanceEdge<S,D>> graph, S sequencing)
//    {
//        graph.addVertex(sequencing);
//    }
//
//    @Override
//    protected Graph<S,DistanceEdge<S,D>> makeEmptyGraph()
//    {
//        return new SparseGraph<S,DistanceEdge<S,D>>(); // we aren't sparse, but JUNG doesn't provide any other concrete impl.
//    }
//
//    @Override
//    protected DistanceEdge<S,D> makeEdge(S node1, S node2, D snpCount)
//    {
//        return new DistanceEdge<S,D>(node1, node2, snpCount);
//    }
//
//    @Override
//    protected Set<Set<DistanceEdge<S,D>>> inferMinTrees(Graph<S,DistanceEdge<S,D>> completeGraph)
//    {
//        PrimMinimumSpanningTree<S,DistanceEdge<S,D>> mst =
//                new PrimMinimumSpanningTree<S, DistanceEdge<S,D>>(new Factory<Graph<S,DistanceEdge<S,D>>>()
//                {
//                    public Graph<S, DistanceEdge<S,D>> create()
//                    {
//                        return makeEmptyGraph();
//                    }
//                },
//                new Transformer<DistanceEdge<S,D>,Double>()
//                {
//                    public Double transform(DistanceEdge edge)
//                    {
//                        D snpCount = getSnpCount(edge);
//                        return DtoDouble(snpCount);
//                    }
//                });
//
//        HashSet<DistanceEdge<S,D>> edgeOfOneMst = new HashSet<DistanceEdge<S,D>>(mst.transform(completeGraph).getEdges());
//        HashSet<Set<DistanceEdge<S,D>>> allMsts = new HashSet<Set<DistanceEdge<S,D>>>();
//        allMsts.add(edgeOfOneMst);
//
//        return allMsts;
//    }
//
//    @Override
//    protected D getSnpCount(DistanceEdge<S,D> edge)
//    {
//        return edge.getDistance();
//    }
//
//    protected abstract double DtoDouble(D snpCount);
//
//}
