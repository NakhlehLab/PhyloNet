//package edu.rice.cs.bioinfo.programs.soranus.models.analysis;
//
//import edu.rice.cs.bioinfo.library.graph.algorithms.generation.simple.CompleteGraphFactory;
//import edu.rice.cs.bioinfo.library.graph.algorithms.spanningTree.FindAllMinimumSpanningTreesCombinationsAsc;
//import edu.rice.cs.bioinfo.library.graph.algorithms.spanningTree.GraphDisconnectedException;
//import edu.rice.cs.bioinfo.library.programming.Tuple;
//
//import java.util.HashSet;
//import java.util.Set;
//
///**
// * Created with IntelliJ IDEA.
// * User: Matt
// * Date: 4/12/13
// * Time: 5:56 PM
// * To change this template use File | Settings | File Templates.
// */
//public abstract class MinSpanTreesSnpInferrerBacktrack<S> implements MinSpanTreesSnpInferrer<S,SequencingEdge<S>,GraphDisconnectedException>
//{
//
//
//    public Set<Set<SequencingEdge<S>>> inferMinTrees(final Set<S> sequencings) throws GraphDisconnectedException
//    {
//        final Set<SequencingEdge<S>> edgeOfCompleteGraph = new CompleteGraphFactory<S,SequencingEdge<S>>()
//        {
//            @Override
//            public SequencingEdge makeEdge(S node1, S node2)
//            {
//                return new SequencingEdge(node1, node2, getSnpDistance(node1, node2));
//            }
//
//        }.makeCompleteGraph(sequencings);
//
//        final Set<Set<SequencingEdge<S>>> minSpanningTrees = new HashSet<Set<SequencingEdge<S>>>();
//        FindAllMinimumSpanningTreesCombinationsAsc<Void,S,SequencingEdge<S>,Long> backtrack =
//                new FindAllMinimumSpanningTreesCombinationsAsc<Void,S,SequencingEdge<S>,Long>()
//        {
//            @Override
//            protected Set<S> getNodesOfGraph(Void graph)
//            {
//                return sequencings;
//            }
//
//            @Override
//            protected Set<SequencingEdge<S>> getEdgesOfGraph(Void graph)
//            {
//                return edgeOfCompleteGraph;
//            }
//
//            @Override
//            protected Iterable<? extends SequencingEdge<S>> getIncidentEdges(S node, Void graph)
//            {
//                HashSet<SequencingEdge<S>> incidentEdges = new HashSet<SequencingEdge<S>>();
//                for(SequencingEdge edge : edgeOfCompleteGraph)
//                {
//                    if(edge.Sequencing1.equals(node) || edge.Sequencing2.equals(node))
//                        incidentEdges.add(edge);
//                }
//
//                return incidentEdges;
//            }
//
//            @Override
//            protected Tuple<? extends S, ? extends S> getNodesOfEdge(SequencingEdge<S> edge, Void graph)
//            {
//                return new Tuple<S, S>(edge.Sequencing1, edge.Sequencing2);
//            }
//
//            @Override
//            protected Long add(Long term1, Long term2)
//            {
//                return term1 + term2;
//            }
//
//            @Override
//            protected Long makeZero()
//            {
//                return 0L;
//            }
//
//            @Override
//            protected Long getWeight(SequencingEdge edge)
//            {
//                return edge.SnpDistance;  //To change body of implemented methods use File | Settings | File Templates.
//            }
//        };
//      /*  backtrack.addMinSpanTreeFoundListener(new Proc1<Set<SequencingEdge<S>>>()
//        {
//            public void execute(Set<SequencingEdge<S>> minSpanTree)
//            {
//                minSpanningTrees.add(minSpanTree);
//            }
//        }); */
//        backtrack.execute(null);
//
//        return minSpanningTrees;
//    }
//
//    protected abstract Long getSnpDistance(S sequencing1, S sequencing2);
//}
