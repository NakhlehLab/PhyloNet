package edu.rice.cs.bioinfo.programs.mdcnetworksearch;
import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseErrorsException;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.RichNewickReadResult;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.graphbuilding.jung.GraphBuilderDirectedSparse;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.parsers.antlr.ast.*;
import edu.rice.cs.bioinfo.library.phylogenetics.Graph;
import edu.rice.cs.bioinfo.library.phylogenetics.GraphReadOnly;
import edu.rice.cs.bioinfo.library.phylogenetics.PhyloEdge;
import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.rea.ReticulateEdgeAddition;
import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.rea.ReticulateEdgeAdditionInPlace;
import edu.rice.cs.bioinfo.library.phylogenetics.scoring.network.MDCOnNetworkYF;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.HillClimber;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.HillClimberObservable;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.network.rea.*;
import edu.rice.cs.bioinfo.library.phylogenetics.graphadapters.jung.*;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.*;
import edu.rice.cs.bioinfo.library.phylogenetics.search.pseudomcmc.network.srna.SrnaPseudoMetropolisHastings;
import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.uci.ics.jung.graph.DirectedSparseGraph;
import org.apache.log4j.*;
import org.w3c.dom.traversal.NodeIterator;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.math.BigDecimal;
import java.util.*;


/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 6/1/12
 * Time: 4:08 PM
 * To change this template use File | Settings | File Templates.
 */
public class Program
{
    private static final  Func1<String, String> _makeNode = new Func1<String, String>()
    {
        public String execute(String node) {
            return node;  //To change body of implemented methods use File | Settings | File Templates.
        }
    };

    private  static final Func5<String, String, BigDecimal, BigDecimal, BigDecimal, PhyloEdge<String>> _makeEdgeFromBD =
            new Func5<String, String, BigDecimal, BigDecimal, BigDecimal, PhyloEdge<String>>()
            {
                public PhyloEdge<String> execute(String source, String dest, BigDecimal branchLength, BigDecimal support, BigDecimal probability) {
                    return new PhyloEdge<String>(source, dest);
                }
            };

    private static final Func4<String,String,Double,Double, PhyloEdge<String>> _makeEdge = new Func4<String, String, Double, Double, PhyloEdge<String>>() {
        public PhyloEdge<String> execute(String source, String dest, Double bl, Double prob) {
            PhyloEdge<String> tbr = new PhyloEdge<String>(source, dest);
            tbr.setBranchLength(bl);
            tbr.setProbability(prob);
            return tbr;
        }
    };

    private static final Func1<String,String> _getNodeLabel = new Func1<String, String>() {
        public String execute(String input) {
            return input;
        }
    };

    private static final Func2<GraphReadOnly<String,PhyloEdge<String>>, PhyloEdge<String>, Double> _getDistance = new Func2<GraphReadOnly<String, PhyloEdge<String>>, PhyloEdge<String>, Double>() {
        public Double execute(GraphReadOnly<String, PhyloEdge<String>> input1, PhyloEdge<String> input2) {
            return  input2.getBranchLength();
        }
    };

    private static final Func2<GraphReadOnly<String,PhyloEdge<String>>, PhyloEdge<String>, Double> _getProbability = new Func2<GraphReadOnly<String, PhyloEdge<String>>, PhyloEdge<String>, Double>() {
        public Double execute(GraphReadOnly<String, PhyloEdge<String>> input1, PhyloEdge<String> input2) {
            return input2.getProbabilty();
        }
    };

    private static final Func1<PhyloEdge<String>, Tuple<String,String>> _edgeToTuple = new Func1<PhyloEdge<String>, Tuple<String, String>>() {
        public Tuple<String, String> execute(PhyloEdge<String> input) {
            return new Tuple<String, String>(input.Source, input.Destination);
        }
    };

    static class GraphAdapter extends DirectedGraphToGraphAdapter<String,PhyloEdge<String>> implements DeepCopyable<GraphAdapter>
    {

        public GraphAdapter(edu.uci.ics.jung.graph.Graph<String, PhyloEdge<String>> stringPhyloEdgeGraph, Func1<PhyloEdge<String>, Tuple<String, String>> edgeToTuple) {
            super(stringPhyloEdgeGraph, edgeToTuple);
        }
        @Override
        public Iterable<PhyloEdge<String>> getEdges()
        {
            List<PhyloEdge<String>> edges = new ArrayList<PhyloEdge<String>>(Graph.getEdges());
            Collections.shuffle(edges);
            return edges;
        }

        @Override
        public GraphAdapter DeepCopy() {

            edu.uci.ics.jung.graph.Graph<String, PhyloEdge<String>> cloneGraph = new DirectedSparseGraph<String, PhyloEdge<String>>();

            for(String node : this.Graph.getVertices())
            {
                cloneGraph.addVertex(node);
            }

            for(PhyloEdge<String> edge : this.Graph.getEdges())
            {
                PhyloEdge<String> cloneEdge = new PhyloEdge<String>(edge.Source, edge.Destination);
                cloneEdge.setBranchLength(edge.getBranchLength());
                cloneEdge.setProbability(edge.getProbabilty());
                cloneEdge.setSupport(edge.getSupport());

                cloneGraph.addEdge(cloneEdge, cloneEdge.Source, cloneEdge.Destination);
            }

            return new GraphAdapter(cloneGraph, this.edgeToTuple);

        }
    }

    public static void main(String[] args) throws Exception
    {
         long maxExaminations = Long.parseLong(args[0]);

        Appender fbAppender = new FileAppender(new PatternLayout("%m\n") ,"expout-firstBetter.txt", false);
        Appender saAppender = new FileAppender(new PatternLayout("%m\n") ,"expout-steepest.txt", false);
        Appender mhAppender = new FileAppender(new PatternLayout("%m\n") ,"expout-pmh.txt", false);
        final Logger logger = Logger.getLogger(Program.class);


        Random rand = new Random(42);


        for(int algoNum = 0; algoNum<3; algoNum++)
        {
            logger.removeAllAppenders();
            if(algoNum == 0)
            {
                logger.addAppender(saAppender);
            }
            else if(algoNum == 1)
            {
                logger.addAppender(fbAppender);
            }
            else if(algoNum == 2)
            {
                logger.addAppender(mhAppender);
            }
            else
            {
                break;
            }

            for(int trial = 0; trial < 1000; trial++)
            {

                GraphAdapter geneTree1 = new GraphAdapter(readSingleNetwork("(4,(1,(2,3)));"), _edgeToTuple);
                GraphAdapter geneTree2 = new GraphAdapter(readSingleNetwork("(3,(4,(1,2)));"), _edgeToTuple);
                GraphAdapter geneTree3 = new GraphAdapter(readSingleNetwork("(4,(3,(1,2)));"), _edgeToTuple);
                GraphAdapter geneTree4 = new GraphAdapter(readSingleNetwork("(3,(4,(1,2)));"), _edgeToTuple);
                GraphAdapter geneTree5 = new GraphAdapter(readSingleNetwork("(4,(3,(1,2)));"), _edgeToTuple);
                GraphAdapter geneTree6 = new GraphAdapter(readSingleNetwork("((1,4),(2,3));"), _edgeToTuple);
                GraphAdapter geneTree7 = new GraphAdapter(readSingleNetwork("(4,(3,(1,2)));"), _edgeToTuple);
                GraphAdapter geneTree8 = new GraphAdapter(readSingleNetwork("(4,(3,(1,2)));"), _edgeToTuple);
                GraphAdapter geneTree9 = new GraphAdapter(readSingleNetwork("(4,(3,(1,2)));"), _edgeToTuple);
                GraphAdapter geneTree10 = new GraphAdapter(readSingleNetwork("(3,(4,(1,2)));"), _edgeToTuple);

                final List<Graph<String,PhyloEdge<String>>> geneTrees = new ArrayList<Graph<String, PhyloEdge<String>>>();
                geneTrees.add(geneTree1);
                geneTrees.add(geneTree2);
                geneTrees.add(geneTree3);
                geneTrees.add(geneTree4);
                geneTrees.add(geneTree5);
                geneTrees.add(geneTree6);
                geneTrees.add(geneTree7);
                geneTrees.add(geneTree8);
                geneTrees.add(geneTree9);
                geneTrees.add(geneTree10);



                GraphAdapter startTree = new GraphAdapter(readSingleNetwork("(4,(3,(1,2)));"), _edgeToTuple);

                ReticulateEdgeAddition<GraphAdapter, String, PhyloEdge<String>> reaStrategy =
                        new ReticulateEdgeAdditionInPlace(new Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String>()
                        {
                            private int _count = 0;
                            public String execute(DirectedGraphToGraphAdapter<String, PhyloEdge<String>> input) {
                                return "_" + (_count++);
                            }
                        },
                                new Func3 <DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, String, PhyloEdge<String>>()
                                {

                                    public PhyloEdge<String> execute(DirectedGraphToGraphAdapter<String, PhyloEdge<String>> arg1, String source, String dest)
                                    {
                                        return  new PhyloEdge<String>(source, dest);
                                    }
                                });


                Func1<GraphAdapter,Integer> getScore = new Func1<GraphAdapter, Integer>() {


                    public Integer execute(GraphAdapter input)
                    {
                        List<Integer> xls = new MDCOnNetworkYF().countExtraCoal(input, geneTrees, null, _getNodeLabel, _getNodeLabel, _getDistance, _getProbability, _getDistance, _getProbability, _makeEdge, _makeEdge);

                        int sum = 0;
                        for(Integer i : xls)
                        {
                            sum+=i;
                        }
                        return sum;
                    }
                };

                Comparator<Integer> isBetterScore = new  Comparator<Integer>() {


                    public int compare(Integer input1, Integer input2) {
                        return input2.compareTo(input1);
                    }
                };

                final int trialCapture = trial;



                Proc4<GraphAdapter,Integer,Long,Long> betterFound = new Proc4<GraphAdapter,Integer,Long,Long>() {
                    @Override
                    public void execute(GraphAdapter input1, Integer score, Long examinationsCount, Long depth) {
                        logger.info(trialCapture + ", " + depth + ", " + examinationsCount + ", " + score);
                    }
                };

                Proc2<GraphAdapter, Integer> initialFound = new Proc2<GraphAdapter, Integer>() {
                    @Override
                    public void execute(GraphAdapter input1, Integer score) {
                        logger.info(trialCapture + ", 0, 1, " + score);
                    }
                };

                Func2<Integer,Integer,Double> divideScore = new Func2<Integer, Integer, Double>()
                {
                    @Override
                    public Double execute(Integer input1, Integer input2) {
                        return ((double)input1) / ((double)input2);
                    }
                };

                 switch (algoNum)
                 {
                     case 0:
                         searchSteepestAscent(startTree, getScore, isBetterScore, maxExaminations, initialFound, betterFound, reaStrategy);
                         break;
                     case 1:
                         searchFirstBetterAscent(startTree, getScore, isBetterScore, maxExaminations, initialFound, betterFound, reaStrategy);
                         break;
                     case 2:
                         searchPMH(startTree, getScore, divideScore, rand, maxExaminations, initialFound, betterFound, reaStrategy);
                         break;

                 }
            }
        }


    }

    private static void searchPMH(GraphAdapter startTree, Func1<GraphAdapter,Integer> getScore, Func2<Integer,Integer,Double> divideScore,
                                  Random rand, long maxExaminations, Proc2<GraphAdapter,Integer> initialFound,
                                  Proc4<GraphAdapter,Integer,Long,Long> betterFound,
                                  ReticulateEdgeAddition<GraphAdapter,String,PhyloEdge<String>> reaStrategy)
    {
        SrnaPseudoMetropolisHastings<GraphAdapter, String, PhyloEdge<String>,Integer>
                searcher = new SrnaPseudoMetropolisHastings<GraphAdapter, String, PhyloEdge<String>,Integer>(reaStrategy, true);

        searcher.addInitialSolutionScoreComputedListener(initialFound);
        searcher.addBetterSolutionFoundListener(betterFound);

        searcher.search(startTree, getScore, divideScore, false, rand, maxExaminations);
    }

    private static void searchSteepestAscent(GraphAdapter startTree, Func1<GraphAdapter, Integer> getScore, Comparator<Integer> betterScore,
                                             long maxExaminations, Proc2<GraphAdapter, Integer> initialFound,
                                             Proc4<GraphAdapter, Integer, Long, Long> betterFound,
                                             ReticulateEdgeAddition<GraphAdapter, String, PhyloEdge<String>> reaStrategy) {

        HillClimberObservable<GraphAdapter, Integer> climber =
                new ReaHillClimberSteepestAscent<GraphAdapter, String, PhyloEdge<String>,Integer>(reaStrategy, true);

        climber.addInitialSolutionScoreComputedListener(initialFound);
        climber.addBetterSolutionFoundListener(betterFound);

        climber.search(startTree, getScore, betterScore, maxExaminations);
    }

    private static void searchFirstBetterAscent(GraphAdapter startTree, Func1<GraphAdapter, Integer> getScore, Comparator<Integer> betterScore,
                                                long maxExaminations, Proc2<GraphAdapter, Integer> initialFound,
                                                Proc4<GraphAdapter, Integer, Long, Long> betterFound,
                                                ReticulateEdgeAddition<GraphAdapter, String, PhyloEdge<String>> reaStrategy) {

        HillClimberObservable<GraphAdapter, Integer> climber =
                new ReaHillClimberFirstBetter<GraphAdapter, String, PhyloEdge<String>,Integer>(reaStrategy, true);

        climber.addInitialSolutionScoreComputedListener(initialFound);
        climber.addBetterSolutionFoundListener(betterFound);

        climber.search(startTree, getScore, betterScore, maxExaminations);
    }

    private static DirectedSparseGraph<String,PhyloEdge<String>> readSingleNetwork(String richNewickString) throws IOException, CoordinateParseErrorsException
    {
        RichNewickReaderAST_ANTLR rnReader = new RichNewickReaderAST_ANTLR();

        GraphBuilderDirectedSparse<String, PhyloEdge<String>> graphBuilder =
                new GraphBuilderDirectedSparse<String, PhyloEdge<String>>(_makeNode,_makeEdgeFromBD);

        RichNewickReadResult<Networks> result = rnReader.read(new ByteArrayInputStream(richNewickString.getBytes()), graphBuilder);

        if(IterableHelp.count(result.getContextErrors()) > 0)
        {
            throw new RuntimeException(result.getContextErrors().iterator().next().Message);
        }

        if(IterableHelp.count(result.getNetworks().Networks) != 1)
        {
            throw new RuntimeException("More than one network generated.");
        }

        return graphBuilder.Graph;
    }


}
