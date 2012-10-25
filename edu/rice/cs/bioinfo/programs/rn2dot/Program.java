package edu.rice.cs.bioinfo.programs.rn2dot;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.RichNewickReaderAST;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.graphbuilding.jung.GraphBuilderDirectedOrderedSparse;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.parsers.antlr.ast.ANTLRRichNewickParser;
import edu.rice.cs.bioinfo.library.phylogenetics.PhyloEdge2;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func5;

import java.io.ByteArrayInputStream;
import java.math.BigDecimal;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 8/28/12
 * Time: 2:46 PM
 * To change this template use File | Settings | File Templates.
 */
public class Program
{
    private static int _nextNodeName = 0;

    private static Func1<String,String> _makeNode = new Func1<String, String>()
    {
        @Override
        public String execute(String label) {
            if(label == null)
                return "_" + _nextNodeName++;
            return label;
        }
    };

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

    public static void main(String[] args)
    {
        String networkNewick = args[0];

        GraphBuilderDirectedOrderedSparse<String,PhyloEdge2<String,BigDecimal>> graphBuilder = new GraphBuilderDirectedOrderedSparse<String, PhyloEdge2<String,BigDecimal>>(_makeNode, _makeEdge);

        RichNewickReaderAST rnReader =  new RichNewickReaderAST(ANTLRRichNewickParser.MAKE_DEFAULT_PARSER);
        rnReader.readAnyErrorToRuntimeException(new ByteArrayInputStream(networkNewick.getBytes()), graphBuilder);

        Map<String,Integer> nodeLabelToDOTName = new HashMap<String,Integer>();

         System.out.println("digraph network {");

        int nextAvailDotName = 0;
        for(String networkNode : graphBuilder.Graph.getVertices())
        {
           int dotName = nextAvailDotName++;
           nodeLabelToDOTName.put(networkNode, dotName);
           System.out.println(dotName + "[label=\"" + networkNode + "\"];");
        }

         for(PhyloEdge2<String,BigDecimal> edge : graphBuilder.Graph.getEdges())
         {
              String tailNodeDotName =  nodeLabelToDOTName.get(edge.Source).toString();
              String tipNodeDotName = nodeLabelToDOTName.get(edge.Destination).toString();
              System.out.println( tailNodeDotName + "->" + tipNodeDotName + "[label=\"" + edge.getBranchLength() + "\"];");
         }

        System.out.println("}");
    }
}
