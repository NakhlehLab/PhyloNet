package edu.rice.cs.bioinfo.programs.addreticulationedges;


import com.sun.org.apache.bcel.internal.generic.NEW;
import com.sun.xml.internal.ws.server.StatefulInstanceResolver;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.printing.jung.JungRichNewickPrinterCompact;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.graphbuilding.jung.GraphBuilderDirectedOrderedSparse;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.parsers.antlr.ast.RichNewickReaderAST_ANTLR;
import edu.rice.cs.bioinfo.library.phylogenetics.AddReticulationEdge;
import edu.rice.cs.bioinfo.library.phylogenetics.Graph;
import edu.rice.cs.bioinfo.library.phylogenetics.GraphReadOnly;
import edu.rice.cs.bioinfo.library.phylogenetics.PhyloEdge;
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
        return "_" + (_nextNodeNumber++);
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

    private  Func1<Graph<String,PhyloEdge<String>>,String> _makeNodeFromGraph = new Func1<Graph<String,PhyloEdge<String>>,String>()
                        {
                            public String execute(Graph<String, PhyloEdge<String>> input) {
                                return makeNodeName();
                            }
                        };

    private  Func3<Graph<String,PhyloEdge<String>>, String, String, PhyloEdge<String>> _makeEdgeFromGraph = new Func3<Graph<String, PhyloEdge<String>>, String, String, PhyloEdge<String>>() {
        public PhyloEdge<String> execute(Graph<String, PhyloEdge<String>> arg1, String arg2, String arg3) {
            return new PhyloEdge<String>(arg2, arg3);
        }
    };

    private Func5<String, String, BigDecimal, BigDecimal, BigDecimal, PhyloEdge<String>> _makeEdge =
            new  Func5<String, String, BigDecimal, BigDecimal, BigDecimal,PhyloEdge<String>>()
            {
                @Override
                public PhyloEdge<String> execute(String source, String dest, BigDecimal branchLength, BigDecimal arg4, BigDecimal prob) {
                    PhyloEdge<String> edge = new  PhyloEdge<String>(source, dest);
                    edge.setBranchLength(branchLength.doubleValue());
                    return edge;
                }
            };

    private Func2<GraphReadOnly<String,PhyloEdge<String>>, PhyloEdge<String>, Double> _getBranchLength = new Func2<GraphReadOnly<String, PhyloEdge<String>>, PhyloEdge<String>, Double>() {
        public Double execute(GraphReadOnly<String, PhyloEdge<String>> input1, PhyloEdge<String> edge) {
            return edge.getBranchLength();
        }
    };

    private Proc3<GraphReadOnly<String,PhyloEdge<String>>, PhyloEdge<String>, Double> _setBranchLength = new Proc3<GraphReadOnly<String, PhyloEdge<String>>, PhyloEdge<String>, Double>() {
        public void execute(GraphReadOnly<String, PhyloEdge<String>> input1, PhyloEdge<String> edge, Double branchLength) {
            edge.setBranchLength(branchLength);
        }
    };

    private Func<Double> _makeZero = new Func<Double>() {
        public Double execute() {
            return new Double(0);
        }
    };

    private Func2<Double,Double,Double> _add = new Func2<Double, Double, Double>() {
        public Double execute(Double input1, Double input2) {
            return new Double(input1.doubleValue() + input2.doubleValue());
        }
    };

    private Func2<Double,Double,Double> _subtract = new Func2<Double, Double, Double>() {
        public Double execute(Double input1, Double input2) {
            return new Double(input1.doubleValue() - input2.doubleValue());
        }
    };

    private Func1<Double,Double> _half = new Func1<Double, Double>() {
        public Double execute(Double input) {
            return input.doubleValue() / 2.0;
        }
    };

    private  AddReticulationEdge<String,PhyloEdge<String>, Double> _addReticulationEdge =
                new AddReticulationEdge<String,PhyloEdge<String>,Double>(_makeNodeFromGraph, _makeEdgeFromGraph, _getBranchLength, _setBranchLength, _makeZero, _add, _subtract, _half);



    String addEdges(String networkNewick, int numReticulationEdges, double scaleFactor, Random rand)
    {
        GraphBuilderDirectedOrderedSparse<String,PhyloEdge<String>> graphBuilder = new GraphBuilderDirectedOrderedSparse<String, PhyloEdge<String>>(_makeNode, _makeEdge);

        RichNewickReaderAST_ANTLR rnReader = new RichNewickReaderAST_ANTLR();
        rnReader.readAnyErrorToRuntimeException(new ByteArrayInputStream(networkNewick.getBytes()), graphBuilder);

        final DirectedGraphToGraphAdapter<String,PhyloEdge<String>> graph = new DirectedGraphToGraphAdapter(graphBuilder.Graph, new Func1<PhyloEdge<String>, Tuple<String,String>>()
        {

            public Tuple<String, String> execute(PhyloEdge<String> input) {
                return new Tuple<String, String>(input.Source, input.Destination);
            }
        });

        int numReticulationEdgesInNetwork = 0;


        while(numReticulationEdgesInNetwork < numReticulationEdges)
        {
            List<PhyloEdge<String>> edges = IterableHelp.toList(graph.getEdges());
            Collections.shuffle(edges, rand);

            final PhyloEdge<String> e1 = edges.remove(0);
            final PhyloEdge<String> e2 = edges.remove(0);


            try
            {
                _addReticulationEdge.execute(graph, e1, e2, true);
                numReticulationEdgesInNetwork++;
            }
            catch(IllegalArgumentException e)
            {
            }

        }

        Func1<String,String> getLabel = new Func1<String, String>() {
            public String execute(String input) {
                return input;
            }
        };

        // scale edges if requested
        if(scaleFactor != 1.0)
        {
            for(PhyloEdge<String> edge : graph.getEdges())
            {
                edge.setBranchLength(edge.getBranchLength() * scaleFactor);
            }
        }


        StringWriter writer = new StringWriter();
        JungRichNewickPrinterCompact printer = new JungRichNewickPrinterCompact<String>();
        printer.setGetBranchLength(new Func2<String,String,String>()
        {

            public String execute(String source, String dest) {
                for(PhyloEdge<String> edge : graph.Graph.getOutEdges(source))
                {
                    if(edge.Destination.equals(dest))
                    {
                        return edge.getBranchLength() + "";
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
