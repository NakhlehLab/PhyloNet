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

    private int _nextNodeNumber = 1;

    private String makeNodeName()
    {
        return "i" + (_nextNodeNumber++);
    }

     private Func1<String,String> _makeNode = new Func1<String, String>()
    {
        @Override
        public String execute(String label) {
            if(label == null)
                return makeNodeName();
            return label;
        }
    };

    private  Func1<Graph<String,PhyloEdge2<String,BigDecimal>>,String> _makeNodeFromGraph = new Func1<Graph<String,PhyloEdge2<String,BigDecimal>>,String>()
                        {
                            public String execute(Graph<String, PhyloEdge2<String,BigDecimal>> input) {
                                return makeNodeName();
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
                public PhyloEdge2<String,BigDecimal> execute(String source, String dest, BigDecimal branchLength, BigDecimal arg4, BigDecimal prob) {
                    PhyloEdge2<String,BigDecimal> edge = new  PhyloEdge2<String,BigDecimal>(source, dest);
                    edge.setBranchLength(branchLength);
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

    private static Proc3<GraphReadOnly, PhyloEdge2<String,BigDecimal>, BigDecimal> _setProb = new Proc3<GraphReadOnly, PhyloEdge2<String,BigDecimal>, BigDecimal>() {
          public void execute(GraphReadOnly input1, PhyloEdge2<String,BigDecimal> edge, BigDecimal prob) {
              edge.setProbability(prob);
          }
      };

     private static Func2<GraphReadOnly, PhyloEdge2<String,BigDecimal>, Boolean> _isEdgeProbUnset = new Func2<GraphReadOnly, PhyloEdge2<String,BigDecimal>, Boolean>() {
        public Boolean execute(GraphReadOnly input1, PhyloEdge2<String,BigDecimal> edge) {
            return edge.getProbabilty() == null;
        }
    };


    private  AddReticulationEdge<String,PhyloEdge2<String,BigDecimal>, BigDecimal> _addReticulationEdge =
                new AddReticulationEdge<String,PhyloEdge2<String,BigDecimal>,BigDecimal>(_makeNodeFromGraph, _makeEdgeFromGraph, _getBranchLength, _setBranchLength, _makeZero, _add, _subtract, _half);



    String addEdges(String networkNewick, int numEdgesToAdd, double scaleFactor, Random rand)
    {
        GraphBuilderDirectedOrderedSparse<String,PhyloEdge2<String,BigDecimal>> graphBuilder = new GraphBuilderDirectedOrderedSparse<String, PhyloEdge2<String,BigDecimal>>(_makeNode, _makeEdge);

       RichNewickReaderAST rnReader =  new RichNewickReaderAST(ANTLRRichNewickParser.MAKE_DEFAULT_PARSER);
        rnReader.readAnyErrorToRuntimeException(new ByteArrayInputStream(networkNewick.getBytes()), graphBuilder);

        final DirectedGraphToGraphAdapter<String,PhyloEdge2<String,BigDecimal>> graph = new DirectedGraphToGraphAdapter(graphBuilder.Graph, new Func1<PhyloEdge2<String,BigDecimal>, Tuple<String,String>>()
        {

            public Tuple<String, String> execute(PhyloEdge2<String,BigDecimal> input) {
                return new Tuple<String, String>(input.Source, input.Destination);
            }
        });



        int numEdgesAdded = 0;


        while(numEdgesAdded < numEdgesToAdd)
        {
            List<PhyloEdge2<String,BigDecimal>> edges = IterableHelp.toList(graph.getEdges());
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
                for(PhyloEdge2<String,BigDecimal> edge : graph.Graph.getOutEdges(source))
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

                for(PhyloEdge2<String,BigDecimal> edge : graph.Graph.getOutEdges(source))
                {
                    if(edge.Destination.equals(dest))
                    {
                        return edge.getProbabilty() + "";
                    }
                }
                throw new RuntimeException("Unknown edge " + source + " " + dest);
            }
        });
        printer.print(graph.Graph, getLabel, writer);

        return writer.toString();
    }

    public static void main(String[] args)
    {


        String networkNewick = args[0];
        int numReticulationEdges = Integer.parseInt(args[1]);

        double scaleFactor = 1.0;
        if(args.length > 2)
        {
            scaleFactor = Double.parseDouble(args[2]);
        }

        String resultNetwork = new Program().addEdges(networkNewick, numReticulationEdges, scaleFactor, new Random());

        System.out.println(resultNetwork);

    }
}
