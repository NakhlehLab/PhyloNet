package edu.rice.cs.bioinfo.programs.networksearchgen;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.printing.jung.JungRichNewickPrinterCompact;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.RichNewickReaderAST;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.graphbuilding.jung.GraphBuilderDirectedOrderedSparse;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.parsers.antlr.ast.ANTLRRichNewickParser;
import edu.rice.cs.bioinfo.library.phylogenetics.GetInDegree;
import edu.rice.cs.bioinfo.library.phylogenetics.PhyloEdge2;
import edu.rice.cs.bioinfo.library.phylogenetics.graphadapters.jung.DirectedGraphToGraphAdapter;
import edu.rice.cs.bioinfo.library.programming.*;

import javax.media.j3d.LineStripArray;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.StringWriter;
import java.math.BigDecimal;
import java.util.LinkedList;
import java.util.Map;
import java.util.Scanner;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 8/21/12
 * Time: 2:34 PM
 * To change this template use File | Settings | File Templates.
 */
public class Program
{
    public static void main(String[] args)
    {

        run(args, new Proc1<String>()
        {

            public void execute(String input) {
                System.out.print(input);
            }
        });
    }

    private static Func5<String, String, BigDecimal, BigDecimal, BigDecimal, PhyloEdge2<String,BigDecimal>> _makeEdge =
            new  Func5<String, String, BigDecimal, BigDecimal, BigDecimal,PhyloEdge2<String,BigDecimal>>()
            {
                @Override
                public PhyloEdge2<String,BigDecimal> execute(String source, String dest, BigDecimal branchLength, BigDecimal arg4, BigDecimal prob) {
                    PhyloEdge2<String,BigDecimal> edge = new  PhyloEdge2<String,BigDecimal>(source, dest);
                    edge.setBranchLength(branchLength);
                    return edge;
                }
            };

    private static Func1<String,String> _makeNode = new Func1<String, String>()
    {
        private int _nodeNum = 0;
        @Override
        public String execute(String label) {
            if(label == null)
                return "_" + _nodeNum++;
            return label;
        }
    };

    public static void run(String[] args, Proc1<String> out)
    {
        if(args[0].equals("ms"))
        {
            int num_gt = Integer.parseInt(args[1]);
            String networkNewick = args[2];
            BigDecimal ultrametricThreshold = new BigDecimal(args[3]);
            Tuple<String,Map<String,Integer>> toMSScriptResult = edu.rice.cs.bioinfo.programs.rn2ms.Program.toMSScript(num_gt, networkNewick, ultrametricThreshold);
            out.execute(toMSScriptResult.Item2.toString() + "\n");
            out.execute(toMSScriptResult.Item1 + "\n");
        }
        else if(args[0].equals("st"))
        {
            String gtText;
            if(new File(args[1]).exists())
            {
                try
                {
                    gtText = new Scanner( new File(args[1]) ).useDelimiter("\\A").next();
                }
                catch(FileNotFoundException e)
                {
                    throw new RuntimeException(e);
                }
            }
            else
            {
                gtText = args[1];
            }

            String[] lines = gtText.split("//");
            for(String line : lines)
            {
                line = line.trim();
                if(line.startsWith("("))
                {
                    final GraphBuilderDirectedOrderedSparse<String,PhyloEdge2<String,BigDecimal>> graphBuilder = new GraphBuilderDirectedOrderedSparse<String, PhyloEdge2<String,BigDecimal>>(_makeNode, _makeEdge);

                    final DirectedGraphToGraphAdapter<String,PhyloEdge2<String,BigDecimal>> graph =
                            new DirectedGraphToGraphAdapter(graphBuilder.Graph, new Func1<PhyloEdge2<String,BigDecimal>, Tuple<String,String>>()
                            {

                                public Tuple<String, String> execute(PhyloEdge2<String, BigDecimal> edge) {
                                    return edge.NodesOfEdge;
                                }
                            });

                    RichNewickReaderAST rnReader =  new RichNewickReaderAST(ANTLRRichNewickParser.MAKE_DEFAULT_PARSER);
                    rnReader.readAnyErrorToRuntimeException(new ByteArrayInputStream(line.getBytes()), graphBuilder);

                    for(PhyloEdge2<String,BigDecimal> edge : graphBuilder.Graph.getEdges())
                    {
                        edge.setBranchLength(edge.getBranchLength().multiply(new BigDecimal("2")));
                    }

                    StringWriter writer = new StringWriter();
                    JungRichNewickPrinterCompact printer = new JungRichNewickPrinterCompact<String>();
                    printer.setGetBranchLength(new Func2<String,String,String>()
                    {

                        public String execute(String source, String dest) {
                            for(PhyloEdge2<String,BigDecimal> edge : graphBuilder.Graph.getOutEdges(source))
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

                            for(PhyloEdge2<String,BigDecimal> edge : graphBuilder.Graph.getOutEdges(source))
                            {
                                if(edge.Destination.equals(dest))
                                {
                                    return edge.getProbability() + "";
                                }
                            }
                            throw new RuntimeException("Unknown edge " + source + " " + dest);
                        }
                    });

                    Func1<String,String> getLabel = new Func1<String, String>() {
                        public String execute(String input) {
                            return input;
                        }
                    };
                    printer.print(graphBuilder.Graph, getLabel, writer);

                    out.execute(writer.toString() + "\n");
                }
            }
        }
        else
        {
            throw new IllegalArgumentException(args[0]);
        }
    }
}
