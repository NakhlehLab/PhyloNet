package edu.rice.cs.bioinfo.programs.richnewick2hybridsimnewick;

import edu.rice.cs.bioinfo.library.language.hybridsimnewick._2012_6_22.printing.HybridSimNewickPrinterCompact;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.RichNewickReaderAST;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.graphbuilding.phylolib.GraphBuilderPhyloGraph2;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.parsers.antlr.ast.ANTLRRichNewickParser;
import edu.rice.cs.bioinfo.library.phylogenetics.*;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func1Identity;
import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.Tuple;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.StringWriter;
import java.math.BigDecimal;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/24/12
 * Time: 11:25 AM
 * To change this template use File | Settings | File Templates.
 */
public class Program
{
    public static void main(String[] args) throws IOException
    {
        if(args.length != 2)
        {
            System.err.println("usage: java -jar [jar file] [rich newick string] [root branch length]");
            System.err.println("found " + args.length + " args.");
            return;
        }
        String richNewickInput = args[0];
        String rootBranchLength = args[1];
        String hybridSimNewickOutput = convertString(richNewickInput, rootBranchLength);


        System.out.println(hybridSimNewickOutput);
    }

    static String convertString(String richNewickInput, String rootBranchLength)
    {
        GraphBuilderPhyloGraph2<String,BigDecimal> graphBuilder = new GraphBuilderPhyloGraph2<String, BigDecimal>(new Func1Identity<String>(), new Func1Identity<BigDecimal>());
        RichNewickReaderAST reader =  new RichNewickReaderAST(ANTLRRichNewickParser.MAKE_DEFAULT_PARSER);
        reader.readAnyErrorToRuntimeException(new ByteArrayInputStream(richNewickInput.getBytes()), graphBuilder);
        final PhyloGraph2<String,BigDecimal> graph = graphBuilder.Graph;

        final Map<String,String> nodeToHybridIndex = new HashMap<String, String>();

        GetInDegree<String, PhyloEdge2<String, BigDecimal>> getInDegree = new GetInDegree<String, PhyloEdge2<String, BigDecimal>>();
        int nextFreeHybidIndex = 1;
        for(String node : graph.getNodes())
        {
            if(getInDegree.execute(graph, node) > 1)
            {
                nodeToHybridIndex.put(node, nextFreeHybidIndex + "");
                nextFreeHybidIndex++;
            }
            else
            {
                nodeToHybridIndex.put(node, null);
            }
        }

        Func2<String,String,String> getBranchLength = new Func2<String, String, String>() {
            public String execute(String tail, String tip) {
                return graph.getEdge(tail, tip).getBranchLength().toString();
            }
        };

        Func2<String,String,String> getProbability = new Func2<String, String, String>() {
            public String execute(String tail, String tip) {
                return graph.getEdge(tail, tip).getProbability().toString();
            }
        };

         Func1<String, String> getHybridIndex = new Func1<String, String>() {
             public String execute(String node) {
                 return nodeToHybridIndex.get(node);
             }
         };

        Func1<String, Iterable<String>> getDestinationNodes = new Func1<String, Iterable<String>>() {
            public Iterable<String> execute(String node) {
                return new GetDirectSuccessors<String, PhyloEdge2<String, BigDecimal>>().execute(graph, node);
            }
        };

        Func1<String, Tuple<String,String>> getHybridParents = new Func1<String, Tuple<String, String>>() {
            public Tuple<String, String> execute(String node) {

                Iterator<String> parents = new GetDirectPredecessors<String,PhyloEdge2<String,BigDecimal>>().execute(graph, node).iterator();
                return new Tuple<String, String>(parents.next(), parents.next());

            }
        };

        HybridSimNewickPrinterCompact<String> printer = new HybridSimNewickPrinterCompact<String>(new Func1Identity<String>(), getBranchLength, getProbability, getHybridIndex);

        String rootNode = new FindRoot<String>().execute(graph);

        StringWriter stringAccum = new StringWriter();
        printer.print(rootNode, getDestinationNodes, getHybridParents, rootBranchLength, stringAccum);

        return stringAccum.toString();
    }
}
