package edu.rice.cs.bioinfo.programs.mdcnetworksearch;
import edu.rice.cs.bioinfo.library.language.parsing.CoordinateParseErrorsException;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.printing.jung.JungRichNewickPrinterCompact;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.RichNewickReadResult;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.graphbuilding.jung.GraphBuilderDirectedSparse;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.parsers.antlr.ast.RichNewickReaderAST_ANTLR;
import edu.rice.cs.bioinfo.library.phylogenetics.Graph;
import edu.rice.cs.bioinfo.library.phylogenetics.GraphReadOnly;
import edu.rice.cs.bioinfo.library.phylogenetics.PhyloEdge;
import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.rea.ReticulateEdgeAddition;
import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.rea.ReticulateEdgeAdditionInPlace;
import edu.rice.cs.bioinfo.library.phylogenetics.scoring.network.MDCOnNetworkYF;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.HillClimbResult;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.HillClimberObservable;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.KSteepestAscentBase;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.network.rea.*;
import edu.rice.cs.bioinfo.library.phylogenetics.graphadapters.jung.*;
import edu.rice.cs.bioinfo.library.phylogenetics.search.pseudomcmc.PMHGenerationLimitRestartSearcher;
import edu.rice.cs.bioinfo.library.phylogenetics.search.pseudomcmc.PseudoMetropolisHastingsResult;
import edu.rice.cs.bioinfo.library.phylogenetics.search.pseudomcmc.network.srna.SrnaPseudoMetropolisHastings;
import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.uci.ics.jung.graph.DirectedGraph;
import edu.uci.ics.jung.graph.DirectedSparseGraph;
import org.apache.log4j.*;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.StringWriter;
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
        private int _count = 0;

        public String execute(String node) {

            if(node == null)
            {
                return "i" + (_count++);
            }
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

        public GraphAdapter(DirectedGraph<String, PhyloEdge<String>> stringPhyloEdgeGraph, Func1<PhyloEdge<String>, Tuple<String, String>> edgeToTuple) {
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

            DirectedGraph<String, PhyloEdge<String>> cloneGraph = new DirectedSparseGraph<String, PhyloEdge<String>>();

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
         long maxExaminations = 2000;

        Appender fbAppender = new FileAppender(new PatternLayout("%m\n") ,"expout-firstBetter.txt", false);
        Appender saAppender = new FileAppender(new PatternLayout("%m\n") ,"expout-steepest.txt", false);
        Appender mhAppender = new FileAppender(new PatternLayout("%m\n") ,"expout-pmh.txt", false);
        Appender kSteepest = new FileAppender(new PatternLayout("%m\n") ,"expout-ksteepest.txt", false);
        final Logger logger = Logger.getLogger(Program.class);


        Random rand = new Random(42);


        for(int algoNum = 0; algoNum<4; algoNum++)
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
            else if(algoNum == 3)
            {
                logger.addAppender(kSteepest);
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


                final Ref<Integer> getScoreCallsCount = new Ref<Integer>(0);
                Func1<GraphAdapter,Integer> getScore = new Func1<GraphAdapter, Integer>() {


                    public Integer execute(GraphAdapter input)
                    {
                        getScoreCallsCount.set(getScoreCallsCount.get() + 1);
                        Func1<String,String> getStringString = new Func1<String, String>() {
                            @Override
                            public String execute(String input) {
                                return input;
                            }
                        };
                        StringWriter w = new StringWriter();
                        new JungRichNewickPrinterCompact().print( (DirectedGraph) input.Graph, getStringString, w);
                        String rn = w.toString();
                        List<Integer> xls = new MDCOnNetworkYF().countExtraCoal(input, geneTrees, null, _getNodeLabel, _getNodeLabel, _getDistance, _getProbability, _getDistance, _getProbability, _makeEdge, _makeEdge);

                        int sum = 0;
                        for(Integer i : xls)
                        {
                            sum+=i;
                        }
                        if(sum == 0)
                        {
                            int j = 0;
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
                        logger.info(trialCapture + ", " + depth + ", " + getScoreCallsCount.get() + ", " + score);
                    }
                };

                Proc2<GraphAdapter, Integer> initialFound = new Proc2<GraphAdapter, Integer>() {
                    @Override
                    public void execute(GraphAdapter input1, Integer score) {
                        logger.info(trialCapture + ", 0, " + getScoreCallsCount.get() + ", " + score);
                    }
                };

                Func2<Integer,Integer,Double> divideScore = new Func2<Integer, Integer, Double>()
                {
                    @Override
                    public Double execute(Integer input1, Integer input2) {
                        return ((double)input1) / ((double)input2);
                    }
                };

                Func2<Integer,Integer,Integer> getScoreDelta = new Func2<Integer, Integer, Integer>() {
                    @Override
                    public Integer execute(Integer input1, Integer input2) {
                        return input2 - input1;
                    }
                };

                long numExams = -1;
                long generation = -1;
                int minScore = -1;
                 switch (algoNum)
                 {

                     case 0:
                         HillClimbResult<GraphAdapter,Integer> result1 = searchSteepestAscent(startTree, getScore, isBetterScore, maxExaminations, initialFound, betterFound, reaStrategy);
                         numExams = result1.ExaminationsCount;
                         generation = result1.GenerationCount;
                         minScore = result1.BestExaminedScore;
                         break;
                     case 1:
                         HillClimbResult<GraphAdapter,Integer> result2 = searchFirstBetterAscent(startTree, getScore, isBetterScore, maxExaminations, initialFound, betterFound, reaStrategy);
                         numExams = result2.ExaminationsCount;
                         generation = result2.GenerationCount;
                         minScore = result2.BestExaminedScore;
                         break;
                     case 2:
                         PseudoMetropolisHastingsResult<GraphAdapter, Integer> result3 = searchPMH(startTree, getScore, divideScore, rand, maxExaminations, initialFound, betterFound, reaStrategy);
                         numExams = result3.ExaminationsCount;
                         generation = result3.GenerationCount;
                         minScore = result3.BestExaminedScore;
                         break;
                     case 3:
                        HillClimbResult<GraphAdapter,Integer> result4 = searchKSteepest(startTree, getScore, isBetterScore, getScoreDelta, maxExaminations, initialFound, betterFound, reaStrategy);
                        numExams = result4.ExaminationsCount;
                         generation = result4.GenerationCount;
                         minScore = result4.BestExaminedScore;
                         break;

                 }

                logger.info("*" + trialCapture + ", " + generation + ", " + numExams + ", " + minScore);
            }
        }

        fbAppender.close();
        saAppender.close();
        mhAppender.close();


    }

    private static HillClimbResult<GraphAdapter,Integer> searchKSteepest(GraphAdapter startTree, Func1<GraphAdapter, Integer> getScore, Comparator<Integer> betterScore, Func2<Integer, Integer, Integer> getScoreDelta, long maxExaminations, Proc2<GraphAdapter, Integer> initialFound, Proc4<GraphAdapter, Integer, Long, Long> betterFound, ReticulateEdgeAddition<GraphAdapter, String, PhyloEdge<String>> reaStrategy)
    {
        KSteepestAscentBase<GraphAdapter,Integer> searcher = new SrnaKSteepestAscent(50, getScoreDelta, reaStrategy, true);
        searcher.addInitialSolutionScoreComputedListener(initialFound);
        searcher.addBetterSolutionFoundListener(betterFound);
        return searcher.search(startTree, getScore, betterScore, maxExaminations);
    }

    private static PseudoMetropolisHastingsResult<GraphAdapter, Integer> searchPMH(GraphAdapter startTree, Func1<GraphAdapter,Integer> getScore, Func2<Integer,Integer,Double> divideScore,
                                  Random rand, long maxExaminations, final Proc2<GraphAdapter,Integer> initialFound,
                                  final Proc4<GraphAdapter,Integer,Long,Long> betterFound,
                                  final ReticulateEdgeAddition<GraphAdapter,String,PhyloEdge<String>> reaStrategy)
    {

        int generationLimit = 10;

        Func<SrnaPseudoMetropolisHastings<GraphAdapter, String, PhyloEdge<String>,Integer>> makeSearcher = new Func<SrnaPseudoMetropolisHastings<GraphAdapter, String, PhyloEdge<String>, Integer>>() {
            @Override
            public SrnaPseudoMetropolisHastings<GraphAdapter, String, PhyloEdge<String>, Integer> execute() {

                SrnaPseudoMetropolisHastings<GraphAdapter, String, PhyloEdge<String>,Integer>
                searcher = new SrnaPseudoMetropolisHastings<GraphAdapter, String, PhyloEdge<String>,Integer>(reaStrategy, true);

                searcher.addInitialSolutionScoreComputedListener(initialFound);
                searcher.addBetterSolutionFoundListener(betterFound);

                return searcher;  //To change body of implemented methods use File | Settings | File Templates.
            }
        };



        return new PMHGenerationLimitRestartSearcher(makeSearcher).search(startTree, getScore, divideScore, false, rand, maxExaminations, generationLimit);
    }

    private static HillClimbResult<GraphAdapter,Integer> searchSteepestAscent(GraphAdapter startTree, Func1<GraphAdapter, Integer> getScore, Comparator<Integer> betterScore,
                                             long maxExaminations, Proc2<GraphAdapter, Integer> initialFound,
                                             Proc4<GraphAdapter, Integer, Long, Long> betterFound,
                                             ReticulateEdgeAddition<GraphAdapter, String, PhyloEdge<String>> reaStrategy) {

        HillClimberObservable<GraphAdapter, Integer> climber =
                new ReaHillClimberSteepestAscent<GraphAdapter, String, PhyloEdge<String>,Integer>(reaStrategy, true);

        climber.addInitialSolutionScoreComputedListener(initialFound);
        climber.addBetterSolutionFoundListener(betterFound);

        return climber.search(startTree, getScore, betterScore, maxExaminations);
    }

    private static HillClimbResult<GraphAdapter,Integer> searchFirstBetterAscent(GraphAdapter startTree, Func1<GraphAdapter, Integer> getScore, Comparator<Integer> betterScore,
                                                long maxExaminations, Proc2<GraphAdapter, Integer> initialFound,
                                                Proc4<GraphAdapter, Integer, Long, Long> betterFound,
                                                ReticulateEdgeAddition<GraphAdapter, String, PhyloEdge<String>> reaStrategy) {

        HillClimberObservable<GraphAdapter, Integer> climber =
                new ReaHillClimberFirstBetter<GraphAdapter, String, PhyloEdge<String>,Integer>(reaStrategy, true);

        climber.addInitialSolutionScoreComputedListener(initialFound);
        climber.addBetterSolutionFoundListener(betterFound);

        return climber.search(startTree, getScore, betterScore, maxExaminations);
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
