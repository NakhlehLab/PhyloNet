package edu.rice.cs.bioinfo.programs.rn2ms;

import com.sun.org.apache.bcel.internal.generic.RETURN;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.graphbuilding.jung.GraphBuilderDirectedSparse;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.parsers.antlr.ast.RichNewickReaderAST_ANTLR;
import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.uci.ics.jung.algorithms.shortestpath.DijkstraDistance;
import edu.uci.ics.jung.graph.DirectedGraph;
import edu.uci.ics.jung.graph.util.Pair;
import org.apache.commons.collections15.Transformer;

import javax.swing.text.html.HTMLDocument;
import java.io.ByteArrayInputStream;
import java.io.Console;
import java.io.Reader;
import java.math.BigDecimal;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 6/19/12
 * Time: 1:38 PM
 * To change this template use File | Settings | File Templates.
 */
public class Program
{
    static class NetworkNode implements Comparable<NetworkNode>
    {
        public BigDecimal BackTime;

        public String RNLabel = "";

        @Override
        public int compareTo(NetworkNode o) {
            return this.BackTime.compareTo(o.BackTime);
        }
    }

    static class NetworkEdge
    {
        public BigDecimal BranchLength;

        public Integer PopulationNumber;
    }

    private static  int _nextPopNumber = 1;

    private static Func1<String,NetworkNode> _makeNode = new Func1<String, NetworkNode>()
    {
        @Override
        public NetworkNode execute(String label) {
            NetworkNode tbr = new NetworkNode();
            if(label != null)
                tbr.RNLabel = label;
            return tbr;
        }
    };

    private static  Func5<NetworkNode, NetworkNode, BigDecimal, BigDecimal, BigDecimal, NetworkEdge> _makeEdge =
            new  Func5<NetworkNode, NetworkNode, BigDecimal, BigDecimal, BigDecimal,NetworkEdge>()
            {
                @Override
                public NetworkEdge execute(NetworkNode arg1, NetworkNode arg2, BigDecimal branchLength, BigDecimal arg4, BigDecimal arg5) {
                    NetworkEdge edge = new NetworkEdge();
                    edge.BranchLength = branchLength;
                    return edge;
                }
            };

    public static void main(String[] args)
    {
        if(args.length != 1)
        {
            System.err.println("Usage: java -jar rn2ms.jar rich_newick_string");
        }

        String networkNewick = args[0];

        GraphBuilderDirectedSparse<NetworkNode,NetworkEdge> graphBuilder = new GraphBuilderDirectedSparse<NetworkNode, NetworkEdge>(_makeNode, _makeEdge);

        RichNewickReaderAST_ANTLR rnReader = new RichNewickReaderAST_ANTLR();
        rnReader.readAnyErrorToRuntimeException(new ByteArrayInputStream(networkNewick.getBytes()), graphBuilder);

        DirectedGraph<NetworkNode,NetworkEdge> graph = graphBuilder.Graph;


        Map<String,Integer> taxonLabelToPopNumber = new HashMap<String,Integer>();
        NetworkNode root = null;
        for(NetworkNode n : graph.getVertices())
        {
            if(graph.getOutEdges(n).size() == 0) // leaf node
            {
                taxonLabelToPopNumber.put(n.RNLabel, _nextPopNumber);
                _nextPopNumber++;
            }

            if(graph.getInEdges(n).size()  == 0) // root node
            {
                root = n;
            }
        }

        assignBackTimes(graph, root);

        List<NetworkNode> netNodesByBackTimeAsc = new LinkedList<NetworkNode>(graph.getVertices());
        Collections.sort(netNodesByBackTimeAsc);


        assignPopulationNumbersToEdges(netNodesByBackTimeAsc, graph, taxonLabelToPopNumber);
        List<PopCommand> commands = generateSplitAndJoinCommands(netNodesByBackTimeAsc, graph, taxonLabelToPopNumber.keySet());

        System.out.println("Node to population:");
        System.out.println(taxonLabelToPopNumber);
        System.out.println("MS command:");

        for(PopCommand command : commands)
        {
            command.execute(new PopCommandAlgo<Object, Object, RuntimeException>() {
                public Object forSplit(Split split, Object input) throws RuntimeException {
                    System.out.print("-es " + split.BackTime + " " + split.OldPopulationNumber + " ");
                    return null;  //To change body of implemented methods use File | Settings | File Templates.
                }

                public Object forJoin(Join join, Object input) throws RuntimeException {
                    System.out.print("-ej " + join.BackTime + " " + join.SecondaryPopulationNumber + " " + join.PrimaryPopulationNumber + " ");
                    return null;  //To change body of implemented methods use File | Settings | File Templates.
                }
            }, null);
        }



    }

    private static List<PopCommand> generateSplitAndJoinCommands(List<NetworkNode> netNodesByBackTimeAsc, DirectedGraph<NetworkNode,NetworkEdge> graph, Set<String> taxonLabels) {

        LinkedList<PopCommand> commands = new LinkedList<PopCommand>();
        for(NetworkNode node : netNodesByBackTimeAsc)
        {
            List<NetworkEdge> inEdges = new ArrayList<NetworkEdge>(graph.getInEdges(node));
            List<NetworkEdge> outEdges = new ArrayList<NetworkEdge>(graph.getOutEdges(node));
            String command;
            if(inEdges.size() == 0)
            {
                commands.add(new Join(node.BackTime, outEdges.get(1).PopulationNumber, outEdges.get(0).PopulationNumber));
            }
            else if(inEdges.size() == 1 && !taxonLabels.contains(node.RNLabel))
            {
                int primaryPopNumber = inEdges.iterator().next().PopulationNumber;
                int secondaryPopNumber = (primaryPopNumber == outEdges.get(0).PopulationNumber ?
                                outEdges.get(1).PopulationNumber :
                                outEdges.get(0).PopulationNumber);
                commands.add(new Join(node.BackTime, primaryPopNumber, secondaryPopNumber));

            }
            else if(inEdges.size() == 2 )
            {
                int oldPopNumber = outEdges.get(0).PopulationNumber;
                int newPopNumber = (inEdges.get(0).PopulationNumber == oldPopNumber ?
                                inEdges.get(1).PopulationNumber :
                                inEdges.get(0).PopulationNumber);
                commands.add(new Split(node.BackTime, oldPopNumber, newPopNumber));
            }

        }

        return commands;
    }

    private static void assignPopulationNumbersToEdges(List<NetworkNode> netNodesByBackTimeAsc, DirectedGraph<NetworkNode,NetworkEdge> graph, Map<String, Integer> taxonLabelToPopNumber) {

         netNodesByBackTimeAsc = new LinkedList<NetworkNode>(netNodesByBackTimeAsc);

        while(netNodesByBackTimeAsc.size() > 1)
        {
            BigDecimal headBackTime = netNodesByBackTimeAsc.get(0).BackTime;
            for(int i = 0; i<netNodesByBackTimeAsc.size() && netNodesByBackTimeAsc.get(i).BackTime.equals(headBackTime); i++)
            {
                NetworkNode iNode = netNodesByBackTimeAsc.get(i);

                for(NetworkEdge outEdge : graph.getOutEdges(iNode))
                {
                    if(outEdge.PopulationNumber == null)
                    {
                        continue;
                    }
                }

                Collection<NetworkEdge> inEdges = graph.getInEdges(iNode);
                if(inEdges.size() == 1)
                {
                    if(taxonLabelToPopNumber.containsKey(iNode.RNLabel))
                    {
                        inEdges.iterator().next().PopulationNumber = taxonLabelToPopNumber.get(iNode.RNLabel);
                    }
                    else
                    {
                        inEdges.iterator().next().PopulationNumber = graph.getOutEdges(iNode).iterator().next().PopulationNumber;
                    }
                    netNodesByBackTimeAsc.remove(iNode);
                }
                else if(inEdges.size() > 1)
                {
                    Iterator<NetworkEdge> inEdgesIterator = inEdges.iterator();
                    NetworkEdge inEdge1 = inEdgesIterator.next();

                    inEdge1.PopulationNumber = graph.getOutEdges(iNode).iterator().next().PopulationNumber;

                    while(inEdgesIterator.hasNext())
                    {
                        inEdgesIterator.next().PopulationNumber = _nextPopNumber;
                        _nextPopNumber++;
                    }
                    netNodesByBackTimeAsc.remove(iNode);
                }
                else
                {
                    throw new RuntimeException("Non root nodes must have 1 or 2 in edges.  Violation: " + iNode.RNLabel);
                }
            }

        }
    }

    private static void assignBackTimes(DirectedGraph<NetworkNode, NetworkEdge> graph, NetworkNode root)
    {
        Transformer<NetworkEdge, BigDecimal> edgeToBranchLength = new Transformer<NetworkEdge, BigDecimal>() {
            @Override
            public BigDecimal transform(NetworkEdge edge) {
                return edge.BranchLength;
            }
        };
        DijkstraDistance<NetworkNode, NetworkEdge> dd = new DijkstraDistance<NetworkNode, NetworkEdge>(graph, edgeToBranchLength);

        Map<NetworkNode,BigDecimal> nodeToForwardTimeValue = new HashMap<NetworkNode, BigDecimal>();

        BigDecimal rootToLeafForwardTime = null;
        for(NetworkNode node : graph.getVertices())
        {
            BigDecimal forwardTime = new BigDecimal(dd.getDistance(root, node).doubleValue());
            nodeToForwardTimeValue.put(node, forwardTime);

            if(graph.getOutEdges(node).size() == 0)
            {
                rootToLeafForwardTime = forwardTime;
            }

        }

        for(NetworkNode node : graph.getVertices())
        {
            node.BackTime =  rootToLeafForwardTime.subtract(nodeToForwardTimeValue.get(node));
        }

    }
}
