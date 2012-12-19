package edu.rice.cs.bioinfo.programs.addreticulationedges;


import com.sun.org.apache.bcel.internal.generic.NEW;
import com.sun.xml.internal.ws.server.StatefulInstanceResolver;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.printing.jung.JungRichNewickPrinterCompact;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.RichNewickReaderAST;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.graphbuilding.jung.GraphBuilderDirectedOrderedSparse;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.parsers.antlr.ast.ANTLRRichNewickParser;
import edu.rice.cs.bioinfo.library.phylogenetics.*;
import edu.rice.cs.bioinfo.library.phylogenetics.graphadapters.jung.DirectedGraphToGraphAdapter;
import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.uci.ics.jung.graph.DirectedGraph;

import javax.print.attribute.standard.MediaSize;
import java.io.ByteArrayInputStream;
import java.io.StringWriter;
import java.math.BigDecimal;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 6/26/12
 * Time: 12:15 PM
 * To change this template use File | Settings | File Templates.
 */
public class Program {

    private  Func1<Graph<String,PhyloEdge2<String,BigDecimal>>,String> _makeNodeFromGraph = new Func1<Graph<String,PhyloEdge2<String,BigDecimal>>,String>()
    {
        private int _nextNodeNumber = 1;

        private String makeNodeName()
        {
                return "i" + (_nextNodeNumber++);
        }

        public String execute(Graph<String, PhyloEdge2<String,BigDecimal>> network) {
            List<String> existingNodeNames = IterableHelp.toList(network.getNodes());
            String nodeName = makeNodeName();

            while(existingNodeNames.contains(nodeName))
            {
                nodeName = makeNodeName();
            }
            return  nodeName;
        }
    };

    private  Func3<Graph<String,PhyloEdge2<String,BigDecimal>>, String, String, PhyloEdge2<String,BigDecimal>> _makeEdgeFromGraph = new Func3<Graph<String, PhyloEdge2<String,BigDecimal>>, String, String, PhyloEdge2<String,BigDecimal>>() {
        public PhyloEdge2<String,BigDecimal> execute(Graph<String, PhyloEdge2<String,BigDecimal>> arg1, String arg2, String arg3) {
            return new PhyloEdge2<String,BigDecimal>(arg2, arg3);
        }
    };

    private Func5<String, String, BigDecimal, BigDecimal, BigDecimal, PhyloEdge2<String,BigDecimal>> _makeEdge =
            new  Func5<String, String, BigDecimal, BigDecimal, BigDecimal,PhyloEdge2<String,BigDecimal>>()
            {
                @Override
                public PhyloEdge2<String,BigDecimal> execute(String source, String dest, BigDecimal branchLength, BigDecimal support, BigDecimal prob) {
                    PhyloEdge2<String,BigDecimal> edge = new  PhyloEdge2<String,BigDecimal>(source, dest);
                    edge.setBranchLength(branchLength);
                    edge.setSupport(support);
                    edge.setProbability(prob);
                    return edge;
                }
            };

    private Func2<GraphReadOnly<String,PhyloEdge2<String,BigDecimal>>, PhyloEdge2<String,BigDecimal>, BigDecimal> _getBranchLength = new Func2<GraphReadOnly<String, PhyloEdge2<String,BigDecimal>>, PhyloEdge2<String,BigDecimal>, BigDecimal>() {
        public BigDecimal execute(GraphReadOnly<String, PhyloEdge2<String,BigDecimal>> input1, PhyloEdge2<String,BigDecimal> edge) {
            return edge.getBranchLength();
        }
    };

    private Proc3<GraphReadOnly<String,PhyloEdge2<String,BigDecimal>>, PhyloEdge2<String,BigDecimal>, BigDecimal> _setBranchLength = new Proc3<GraphReadOnly<String, PhyloEdge2<String,BigDecimal>>, PhyloEdge2<String,BigDecimal>, BigDecimal>() {
        public void execute(GraphReadOnly<String, PhyloEdge2<String,BigDecimal>> input1, PhyloEdge2<String,BigDecimal> edge, BigDecimal branchLength) {
            edge.setBranchLength(branchLength);
        }
    };

    private Func<BigDecimal> _makeZero = new Func<BigDecimal>() {
        public BigDecimal execute() {
            return BigDecimal.ZERO;
        }
    };

    private Func2<BigDecimal,BigDecimal,BigDecimal> _add = new Func2<BigDecimal, BigDecimal, BigDecimal>() {
        public BigDecimal execute(BigDecimal input1, BigDecimal input2) {
            return input1.add(input2);
        }
    };

    private Func2<BigDecimal,BigDecimal,BigDecimal> _subtract = new Func2<BigDecimal, BigDecimal, BigDecimal>() {
        public BigDecimal execute(BigDecimal input1, BigDecimal input2) {
            return input1.subtract(input2);
        }
    };

    private Func1<BigDecimal,BigDecimal> _half = new Func1<BigDecimal, BigDecimal>() {
        public BigDecimal execute(BigDecimal input) {
            return input.divide(new BigDecimal(2));
        }
    };

    private static Proc3<GraphReadOnly<String,PhyloEdge2<String,BigDecimal>>, PhyloEdge2<String,BigDecimal>, BigDecimal> _setProb = new Proc3<GraphReadOnly<String,PhyloEdge2<String,BigDecimal>>, PhyloEdge2<String,BigDecimal>, BigDecimal>() {
        public void execute(GraphReadOnly<String,PhyloEdge2<String,BigDecimal>> input1, PhyloEdge2<String,BigDecimal> edge, BigDecimal prob) {
            edge.setProbability(prob);
        }
    };

    private Func2<GraphReadOnly<String,PhyloEdge2<String,BigDecimal>>, PhyloEdge2<String,BigDecimal>, BigDecimal> _getProb = new Func2<GraphReadOnly<String, PhyloEdge2<String,BigDecimal>>, PhyloEdge2<String,BigDecimal>, BigDecimal>() {
           public BigDecimal execute(GraphReadOnly<String, PhyloEdge2<String,BigDecimal>> input1, PhyloEdge2<String,BigDecimal> edge) {
               return edge.getProbability();
           }
       };

    private static Func2<GraphReadOnly, PhyloEdge2<String,BigDecimal>, Boolean> _isEdgeProbUnset = new Func2<GraphReadOnly, PhyloEdge2<String,BigDecimal>, Boolean>() {
        public Boolean execute(GraphReadOnly input1, PhyloEdge2<String,BigDecimal> edge) {
            return edge.getProbability() == null;
        }
    };


    private  AddReticulationEdge<String,PhyloEdge2<String,BigDecimal>, BigDecimal> _addReticulationEdge =
            new AddReticulationEdge<String,PhyloEdge2<String,BigDecimal>,BigDecimal>(_makeNodeFromGraph, _makeEdgeFromGraph, _getBranchLength, _setBranchLength,
                                                                                     _getProb, _setProb, _makeZero, _add, _subtract, _half);


    private DirectedGraph<String,PhyloEdge2<String,BigDecimal>> parseNewick(final String networkNewick)
    {
        Func1<String,String> makeNode = new Func1<String, String>()
            {
                private int _nextNodeNumber = 1;

            @Override
                public String execute(String label) {
                    if(label == null)
                    {
                        while(networkNewick.contains("i" + _nextNodeNumber))
                        {
                            _nextNodeNumber++;
                        }
                        label = "i" + _nextNodeNumber;
                        _nextNodeNumber++;
                    }
                    return label;
                }
            };

        GraphBuilderDirectedOrderedSparse<String,PhyloEdge2<String,BigDecimal>> graphBuilder = new GraphBuilderDirectedOrderedSparse<String, PhyloEdge2<String,BigDecimal>>(makeNode, _makeEdge);

        RichNewickReaderAST rnReader =  new RichNewickReaderAST(ANTLRRichNewickParser.MAKE_DEFAULT_PARSER);
        rnReader.readAnyErrorToRuntimeException(networkNewick, graphBuilder);

        return graphBuilder.Graph;
    }

    String addEdgesAndEnforceHeight(String networkNewick, int numEdgesToAdd, double desiredHeight, BigDecimal ultrametricThreshold, Random rand)
    {
        DirectedGraph<String,PhyloEdge2<String,BigDecimal>> digraph = parseNewick(networkNewick);
        final Graph<String, PhyloEdge2<String,BigDecimal>> graph = new DirectedGraphToGraphAdapter(digraph,
                new Func1<PhyloEdge2<String,BigDecimal>, Tuple<String,String>>()
                {

                    public Tuple<String, String> execute(PhyloEdge2<String, BigDecimal> edge) {
                        return edge.NodesOfEdge;
                    }
                });
        IsUltrametric.Result ultrametricResult = new IsUltrametric(PhyloEdge2.GetBranchLength()).execute(graph, ultrametricThreshold);
        if(ultrametricResult.IsUltrametricWithinThreshold)
        {
            double scaleFactor = desiredHeight / ultrametricResult.Height.doubleValue();
            return addEdges(digraph, numEdgesToAdd, scaleFactor, rand);
        }
        else
        {
            throw new IllegalArgumentException("Network does not meet ultrametric threshold.");
        }

    }

    String addEdgesAndScale(String networkNewick, int numEdgesToAdd, double scaleFactor, Random rand)
    {
        DirectedGraph<String,PhyloEdge2<String,BigDecimal>> graph = parseNewick(networkNewick);
        return  addEdges(graph, numEdgesToAdd, scaleFactor, rand);
    }

    private String addEdges(final DirectedGraph<String,PhyloEdge2<String,BigDecimal>> digraph, int numEdgesToAdd, double scaleFactor, Random rand)
    {
        final Graph<String, PhyloEdge2<String,BigDecimal>> graph = new DirectedGraphToGraphAdapter(digraph,
                new Func1<PhyloEdge2<String,BigDecimal>, Tuple<String,String>>()
                {

                    public Tuple<String, String> execute(PhyloEdge2<String, BigDecimal> edge) {
                        return edge.NodesOfEdge;
                    }
                });

        int numEdgesAdded = 0;


        while(numEdgesAdded < numEdgesToAdd)
        {
            List<PhyloEdge2<String,BigDecimal>> edges = IterableHelp.toList(digraph.getEdges());
            Collections.shuffle(edges, rand);

            final PhyloEdge2<String,BigDecimal> e1 = edges.remove(0);
            final PhyloEdge2<String,BigDecimal> e2 = edges.remove(0);



            try
            {
                _addReticulationEdge.execute(graph, e1, e2, false);
                numEdgesAdded++;
            }
            catch(IllegalArgumentException e)
            {
            }

        }


        AssignProbToUnnotatedHybridEdges.setRandom(rand);
        new AssignProbToUnnotatedHybridEdges().execute(graph , AssignProbToUnnotatedHybridEdges.UNIFORM_RANDOM, _setProb, _isEdgeProbUnset);

        Func1<String,String> getLabel = new Func1<String, String>() {
            public String execute(String input) {
                return input;
            }
        };

        // scale edges if requested
        if(scaleFactor != 1.0)
        {
            BigDecimal scaleFactorBD = new BigDecimal(scaleFactor);
            for(PhyloEdge2<String,BigDecimal> edge : graph.getEdges())
            {
                edge.setBranchLength(edge.getBranchLength().multiply(scaleFactorBD));
            }
        }


        StringWriter writer = new StringWriter();
        JungRichNewickPrinterCompact printer = new JungRichNewickPrinterCompact<String>();
        printer.setGetBranchLength(new Func2<String,String,String>()
        {

            public String execute(String source, String dest) {
                for(PhyloEdge2<String,BigDecimal> edge : digraph.getOutEdges(source))
                {
                    if(edge.Destination.equals(dest))
                    {
                        return edge.getBranchLength() + "";
                    }
                }
                throw new RuntimeException("Unknown edge " + source + " " + dest);
            }
        });
        printer.setGetProbability(new Func2<String,String,String>()
        {
            public String execute(String source, String dest) {

                if(new GetInDegree().execute(graph, dest) < 2)
                {
                    return null;
                }

                for(PhyloEdge2<String,BigDecimal> edge : digraph.getOutEdges(source))
                {
                    if(edge.Destination.equals(dest))
                    {
                        return edge.getProbability() + "";
                    }
                }
                throw new RuntimeException("Unknown edge " + source + " " + dest);
            }
        });



        printer.print(digraph, getLabel, writer);

        return writer.toString();
    }

    public static void main(String[] args)
    {


        String networkNewick = args[0];
        int numReticulationEdges = Integer.parseInt(args[1]);


        String resultNetwork;
        if(args.length == 1)
        {
            double scaleFactor = 1.0;
            resultNetwork = new Program().addEdgesAndScale(networkNewick, numReticulationEdges, scaleFactor, new Random());
        }
        else if(args.length == 3)
        {
            String switchArg = args[2];
            if(switchArg.equalsIgnoreCase("-s"))
            {
                double scaleFactor = Double.parseDouble(args[3]);
                resultNetwork = new Program().addEdgesAndScale(networkNewick, numReticulationEdges, scaleFactor, new Random());
            }
            else
            {
                throw new IllegalArgumentException("Unexpected switch " + switchArg);
            }
        }
        else if(args.length == 5)
        {
            String switchArg = args[2];
            if(switchArg.equalsIgnoreCase("-h"))
            {
                double height = Double.parseDouble(args[3]);
                BigDecimal ultrametricThreshold = new BigDecimal(args[4]);
                resultNetwork = new Program().addEdgesAndEnforceHeight(networkNewick, numReticulationEdges, height, ultrametricThreshold, new Random());
            }
            else
            {
                throw new IllegalArgumentException("Unexpected switch " + switchArg);
            }
        }
        else
        {
            throw new IllegalArgumentException("Expected number of program arguments to be 1, 3, or 4.");
        }


        System.out.println(resultNetwork);

    }
}
