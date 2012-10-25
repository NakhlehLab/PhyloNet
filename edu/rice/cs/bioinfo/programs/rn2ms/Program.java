package edu.rice.cs.bioinfo.programs.rn2ms;

import com.sun.org.apache.bcel.internal.generic.RETURN;
import com.sun.org.apache.xpath.internal.NodeSet;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.Network;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.RichNewickReaderAST;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.graphbuilding.jung.GraphBuilderDirectedOrderedSparse;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.parsers.antlr.ast.ANTLRRichNewickParser;
import edu.rice.cs.bioinfo.library.phylogenetics.*;
import edu.rice.cs.bioinfo.library.phylogenetics.graphadapters.jung.*;
import edu.rice.cs.bioinfo.library.programming.*;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.uci.ics.jung.algorithms.shortestpath.DijkstraDistance;
import edu.uci.ics.jung.graph.DirectedGraph;
import edu.uci.ics.jung.graph.util.Pair;
import org.apache.commons.collections15.Transformer;

import javax.media.j3d.IndexedGeometryArray;
import javax.swing.text.html.HTMLDocument;
import java.io.ByteArrayInputStream;
import java.io.Console;
import java.io.PrintStream;
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

        @Override
        public String toString()
        {
            return RNLabel;
        }
    }

    static class NetworkEdge
    {
        private NetworkNode _source, _dest;

        NetworkEdge(NetworkNode source, NetworkNode dest)
        {
            _source = source;
            _dest = dest;
        }

        public BigDecimal BranchLength;


        private BigDecimal _probability;

        public void setProbability(BigDecimal prob)
        {
            _probability = prob;
        }

        public BigDecimal getProbability()
        {
            return _probability;
        }

        public Integer PopulationNumber;

        public Tuple<NetworkNode,NetworkNode> getNodes()
        {
            return new Tuple<NetworkNode, NetworkNode>(_source, _dest);
        }
    }

    private static Func1<NetworkEdge, Tuple<NetworkNode,NetworkNode>> _edgeToTuple = new Func1<NetworkEdge, Tuple<NetworkNode, NetworkNode>>() {
        public Tuple<NetworkNode, NetworkNode> execute(NetworkEdge input) {
            return input.getNodes();
        }
    };


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
                public NetworkEdge execute(NetworkNode source, NetworkNode dest, BigDecimal branchLength, BigDecimal arg4, BigDecimal prob) {
                    NetworkEdge edge = new NetworkEdge(source, dest);
                    edge.BranchLength = branchLength;
                    edge.setProbability(prob);
                    return edge;
                }
            };

    public static void main(String[] args)
    {
        Proc1<String> out = new Proc1<String>() {
            public void execute(String input) {
                   System.out.print(input);
            }
        };

        Proc1<String> error = new Proc1<String>() {
            public void execute(String input) {
                System.err.print(input);
            }
        };

        run(args, out, error);
    }

    public static Tuple<String,Map<String,Integer>> toMSScript(int numGT, String networkNewick, BigDecimal ultrametricThreshold)
    {
         GraphBuilderDirectedOrderedSparse<NetworkNode,NetworkEdge> graphBuilder = new GraphBuilderDirectedOrderedSparse<NetworkNode, NetworkEdge>(_makeNode, _makeEdge);


            RichNewickReaderAST rnReader =  new RichNewickReaderAST(ANTLRRichNewickParser.MAKE_DEFAULT_PARSER);
            rnReader.readAnyErrorToRuntimeException(new ByteArrayInputStream(networkNewick.getBytes()), graphBuilder);


        DirectedGraph<NetworkNode,NetworkEdge> graph = graphBuilder.Graph;



        Map<String,Integer> taxonLabelToPopNumber = new HashMap<String,Integer>();
        NetworkNode root = null;
        HashSet<NetworkNode> leaves = new HashSet<NetworkNode>();
        for(NetworkNode n : graph.getVertices())
        {
            if(graph.getOutEdges(n).size() == 0) // leaf node
            {
                leaves.add(n);
                taxonLabelToPopNumber.put(n.RNLabel, _nextPopNumber);
                _nextPopNumber++;
            }

            if(graph.getInEdges(n).size()  == 0) // root node
            {
                root = n;
            }
        }

        if(!isUltrametric(graph, leaves, ultrametricThreshold))
        {
            throw new IllegalArgumentException("Given newtwork must be ultrametric within threshold.");
        }


        assignBackTimes(graph, root);

        List<NetworkNode> netNodesByBackTimeAsc = new ArrayList<NetworkNode>(graph.getVertices());
        Collections.sort(netNodesByBackTimeAsc);

        assignPopulationNumbersToEdges(netNodesByBackTimeAsc, graph, taxonLabelToPopNumber);
        List<PopCommand> commands = generateSplitAndJoinCommands(netNodesByBackTimeAsc, graph, taxonLabelToPopNumber.keySet());
        Collections.sort(commands);  // already sorted by backtime, but need to tie break the split commands by new pop number




        String ones = " ";

        for(int i = 0; i<taxonLabelToPopNumber.size(); i++)
        {
            ones+= "1 ";
        }

        final StringBuffer msCommand =  new StringBuffer("ms " + taxonLabelToPopNumber.size() + " " + numGT + " -T -I " + taxonLabelToPopNumber.size() + ones);

        for(PopCommand command : commands)
        {
            command.execute(new PopCommandAlgo<Object, Object, RuntimeException>() {
                public Object forSplit(Split split, Object input) throws RuntimeException {
                    BigDecimal backTimeMSUnits = split.BackTime.divide(new BigDecimal("2"));
                    msCommand.append("-es " + backTimeMSUnits + " " + split.OldPopulationNumber + " " + split.OldPopulationProbability + " ");
                    return null;  //To change body of implemented methods use File | Settings | File Templates.
                }

                public Object forJoin(Join join, Object input) throws RuntimeException {
                    BigDecimal backTimeMSUnits = join.BackTime.divide(new BigDecimal("2"));
                    msCommand.append("-ej " + backTimeMSUnits + " " + join.SecondaryPopulationNumber + " " + join.PrimaryPopulationNumber + " ");
                    return null;  //To change body of implemented methods use File | Settings | File Templates.
                }
            }, null);
        }

        return new Tuple<String, Map<String, Integer>>(msCommand.toString(), taxonLabelToPopNumber);

    }

    static void run(String[] args, final Proc1<String> out, Proc1<String> error)
    {
        if(args.length != 3)
        {
            error.execute("Usage: java -jar rn2ms.jar num_gt rich_newick_string ultrametric_threshold");
        }

        int num_gt = Integer.parseInt(args[0]);
        String networkNewick = args[1];
        BigDecimal ultrametricThreshold = new BigDecimal(args[2]);


        try
        {
            Tuple<String,Map<String,Integer>> toMSScriptResult = toMSScript(num_gt, networkNewick, ultrametricThreshold);
            Map<String,Integer> taxonLabelToPopNumber = toMSScriptResult.Item2;
            String msCommand = toMSScriptResult.Item1;

            out.execute("Node to population:\n");
            out.execute(taxonLabelToPopNumber.toString() + "\n");
            out.execute("MS command:\n");
            out.execute(msCommand);
        }
        catch(IllegalArgumentException e)
        {
            error.execute(e.getMessage());
            return;
        }



    }

    private static boolean isUltrametric(DirectedGraph<NetworkNode, NetworkEdge> graph, Set<NetworkNode> leafs, BigDecimal ultrametricThreshold)
    {
        Map<NetworkNode, BigDecimal> pathLengthToLeafs = new HashMap<NetworkNode, BigDecimal>();

        LinkedList<NetworkNode> workingList = new LinkedList<NetworkNode>();

        for(NetworkNode leaf : leafs)
        {
            pathLengthToLeafs.put(leaf, BigDecimal.ZERO);
            workingList.add(leaf);
        }

        while(!workingList.isEmpty())
        {
            NetworkNode destNode = workingList.remove();
            BigDecimal lengthFromDestNodeToLeaf = pathLengthToLeafs.get(destNode);

            for(NetworkEdge incomingEdge : graph.getInEdges(destNode))
            {
                Pair<NetworkNode> nodesOfEdge = graph.getEndpoints(incomingEdge);
                NetworkNode sourceNode = nodesOfEdge.getFirst();
                BigDecimal foundPathLength = pathLengthToLeafs.get(sourceNode);

                BigDecimal expectedPathLength =  lengthFromDestNodeToLeaf.add(incomingEdge.BranchLength);

                if(foundPathLength == null)
                {
                    pathLengthToLeafs.put(sourceNode, expectedPathLength);
                    workingList.add(sourceNode);
                }
                else if(foundPathLength.subtract(expectedPathLength).abs().compareTo(ultrametricThreshold) == 1)
                {
                    return false;
                }
            }

        }


        return true;

    }

    private static BigDecimal computeArbitraryRootToLeafPathLength(DirectedGraph<NetworkNode, NetworkEdge> graph, NetworkNode root) {

        BigDecimal totalPathLength = BigDecimal.ZERO;

        NetworkNode i = root;
        Collection<NetworkNode> iDirectSucs = graph.getSuccessors(i);

        while(iDirectSucs.size() > 0)
        {
            NetworkEdge arbOutEdge = graph.getOutEdges(i).iterator().next();
            NetworkNode nextI = graph.getEndpoints(arbOutEdge).getSecond();
            totalPathLength.add(arbOutEdge.BranchLength);
            i = nextI;
            iDirectSucs = graph.getSuccessors(i);
        }

        return totalPathLength;

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
                int primaryPopNumber = inEdges.get(0).PopulationNumber;
                int secondaryPopNumber = (primaryPopNumber == outEdges.get(0).PopulationNumber ?
                                outEdges.get(1).PopulationNumber :
                                outEdges.get(0).PopulationNumber);
                commands.add(new Join(node.BackTime, primaryPopNumber, secondaryPopNumber));

            }
            else if(inEdges.size() == 2 )
            {
                int oldPopNumber = outEdges.get(0).PopulationNumber;
                NetworkEdge newPopInEdge = (inEdges.get(0).PopulationNumber == oldPopNumber ?
                                inEdges.get(1) :
                                inEdges.get(0));
                commands.add(new Split(node.BackTime, oldPopNumber, newPopInEdge.PopulationNumber, BigDecimal.ONE.subtract(newPopInEdge.getProbability())));
            }

        }

        return commands;
    }

    private static void assignPopulationNumbersToEdges(List<NetworkNode> netNodesByBackTimeAsc, DirectedGraph<NetworkNode,NetworkEdge> graph, Map<String, Integer> taxonLabelToPopNumber) {

       int edgesAssigned = 0;

        // assign all leaf in edges their taxon pop number
        for(NetworkNode node : netNodesByBackTimeAsc)
        {
            if(taxonLabelToPopNumber.containsKey(node.RNLabel))
            {
                Collection<NetworkEdge> inEdges = graph.getInEdges(node);
                inEdges.iterator().next().PopulationNumber = taxonLabelToPopNumber.get(node.RNLabel);
                edgesAssigned++;
            }
        }

        // assign new population number to each split event in chronological order
        for(NetworkNode node : netNodesByBackTimeAsc)
        {
            Collection<NetworkEdge> inEdges = graph.getInEdges(node);
            if(inEdges.size() == 2)
            {
                Iterator<NetworkEdge> inEdgesIterator = inEdges.iterator();
                NetworkEdge inEdge = inEdgesIterator.next();
                inEdge.PopulationNumber =  _nextPopNumber;
                _nextPopNumber++;
                edgesAssigned++;
            }
        }

        while(edgesAssigned < graph.getEdgeCount())
        {
            nodes: for(int i = 0; i<netNodesByBackTimeAsc.size(); i++)
            {
                NetworkNode iNode = netNodesByBackTimeAsc.get(i);

                for(NetworkEdge outEdge : graph.getOutEdges(iNode))
                {
                    if(outEdge.PopulationNumber == null)
                    {
                        continue nodes;
                    }
                }

                Collection<NetworkEdge> inEdges = graph.getInEdges(iNode);
                Collection outEdges = graph.getOutEdges(iNode);
                if(inEdges.size() == 1 && outEdges.size() == 2) // join event
                {
                   inEdges.iterator().next().PopulationNumber = graph.getOutEdges(iNode).iterator().next().PopulationNumber;
                   edgesAssigned++;
                }
                else if(inEdges.size() == 2 && outEdges.size() == 1)  // split event
                {
                    Iterator<NetworkEdge> inEdgesIterator = inEdges.iterator();
                    NetworkEdge inEdge = inEdgesIterator.next();

                    while(inEdge.PopulationNumber != null)
                    {
                        inEdge = inEdgesIterator.next();
                    }

                    inEdge.PopulationNumber = graph.getOutEdges(iNode).iterator().next().PopulationNumber;
                    edgesAssigned++;
                }
            }

        }
    }

    private static void assignBackTimes(DirectedGraph<NetworkNode, NetworkEdge> graph, NetworkNode root)
    {
        HashSet<NetworkNode> toAssign = new HashSet<NetworkNode>(graph.getVertices());

        // assign all leafs back time of 0
        for(NetworkNode node : graph.getVertices())
        {
            if(graph.getOutEdges(node).size() == 0)
            {
                node.BackTime = BigDecimal.ZERO;
                toAssign.remove(node);
            }
        }

        while(toAssign.size() > 0)
        {
            for(NetworkNode node : new LinkedList<NetworkNode>(toAssign))
            {
                for(NetworkEdge outEdge  : graph.getOutEdges(node))
                {
                    NetworkNode directDesc = outEdge.getNodes().Item2;

                    if(directDesc.BackTime != null)
                    {
                        node.BackTime = directDesc.BackTime.add(outEdge.BranchLength);
                        toAssign.remove(node);
                    }
                }
            }
        }

        /*
        Transformer<NetworkEdge, BigDecimal> edgeToBranchLength = new Transformer<NetworkEdge, BigDecimal>() {
            @Override
            public BigDecimal transform(NetworkEdge edge) {
                return edge.BranchLength;
            }
        };
        DijkstraDistance<NetworkNode, NetworkEdge> dd = new DijkstraDistance<NetworkNode, NetworkEdge>(graph, edgeToBranchLength);

        Map<NetworkNode,BigDecimal> nodeToForwardTimeValue = new HashMap<NetworkNode, BigDecimal>();

        BigDecimal longestRootToLeafForwardTime = null;
        for(NetworkNode node : graph.getVertices())
        {
            BigDecimal forwardTime = new BigDecimal(dd.getDistance(root, node).doubleValue());
            nodeToForwardTimeValue.put(node, forwardTime);

            if(graph.getOutEdges(node).size() == 0)
            {
                if(longestRootToLeafForwardTime == null || longestRootToLeafForwardTime.compareTo(forwardTime) == -1)
                {
                    longestRootToLeafForwardTime = forwardTime;
                }
            }

        }

        for(NetworkNode node : graph.getVertices())
        {
            node.BackTime =  longestRootToLeafForwardTime.subtract(nodeToForwardTimeValue.get(node));
        }     */

    }
}
