package edu.rice.cs.bioinfo.programs.bdnetworkgen;


import edu.rice.cs.bioinfo.library.language.richnewick._1_0.printing.jung.JungRichNewickPrinterCompact;
import edu.rice.cs.bioinfo.library.phylogenetics.AddReticulationEdge;
import edu.rice.cs.bioinfo.library.phylogenetics.Graph;
import edu.rice.cs.bioinfo.library.phylogenetics.GraphValidator;
import edu.rice.cs.bioinfo.library.phylogenetics.graphadapters.jung.DirectedGraphToGraphAdapter;
import edu.rice.cs.bioinfo.library.programming.Func;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func3;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.uci.ics.jung.graph.DirectedGraph;
import edu.uci.ics.jung.graph.DirectedSparseMultigraph;
import edu.uci.ics.jung.visualization.renderers.Renderer;

import javax.print.DocFlavor;
import java.io.StringWriter;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 6/26/12
 * Time: 12:15 PM
 * To change this template use File | Settings | File Templates.
 */
public class Program {

    private static int _nextNodeNumber = 0;

    private static String makeNodeName()
    {
        return (_nextNodeNumber++) + "";
    }


    public static void main(String[] args)
    {
        double birthProbability = Double.parseDouble(args[0]);
        int minimumNodeCount = Integer.parseInt(args[1]);
        int numReticulationEdges = Integer.parseInt(args[2]);

        Func1<Tuple<String,String>, Tuple<String,String>> edgeToTuple = new Func1<Tuple<String, String>, Tuple<String, String>>() {
            public Tuple<String, String> execute(Tuple<String, String> input) {
                return input;
            }
        };

        DirectedGraphToGraphAdapter<String, Tuple<String,String>> graph = new DirectedGraphToGraphAdapter(new DirectedSparseMultigraph<Integer, Tuple<String,String>>(), edgeToTuple);

        String root = makeNodeName();
        graph.addNode(root);

        LinkedList<String> expansionNodes = new LinkedList<String>();
        expansionNodes.add(root);

        Random rand = new Random();
        while(!expansionNodes.isEmpty())
        {
            String parent = expansionNodes.remove();
            if(rand.nextDouble() >= birthProbability)
            {
                String child1 = makeNodeName();
                String child2 = makeNodeName();
                graph.addNode(child1);
                graph.addNode(child2);
                graph.addEdge(new Tuple<String, String>(parent, child1));
                graph.addEdge(new Tuple<String, String>(parent, child2));
                expansionNodes.addLast(child1);
                expansionNodes.addLast(child2);
            }
            else if(graph.Graph.getVertexCount() < minimumNodeCount)
            {
                expansionNodes.addLast(parent);
            }
        }

        int numReticulationEdgesInNetwork = 0;
        AddReticulationEdge<String,Tuple<String,String>> addReticulationEdge = new AddReticulationEdge<String,Tuple<String,String>>();
        while(numReticulationEdgesInNetwork < numReticulationEdges)
        {
            List<Tuple<String, String>> edges = IterableHelp.toList(graph.getEdges());
            Collections.shuffle(edges, rand);

            Tuple<String,String> e1 = edges.remove(0);
            Tuple<String,String> e2 = edges.remove(0);

            try
            {
                addReticulationEdge.execute(graph, e1, e2, new Func1<Graph<String,Tuple<String,String>>,String>()
                        {
                            public String execute(Graph<String, Tuple<String, String>> input) {
                                return makeNodeName();
                            }
                        }, new Func3<Graph<String,Tuple<String,String>>,String, String, Tuple<String,String>>()
                        {


                            public Tuple<String, String> execute(Graph<String, Tuple<String, String>> arg1, String arg2, String arg3) {
                                return new Tuple<String, String>(arg2, arg3);
                            }
                        }, true);
                numReticulationEdgesInNetwork++;
            }
            catch(IllegalArgumentException e)
            {
            }

        }



        GraphValidator.assertValidGraph(graph);

        Func1<String,String> getLabel = new Func1<String, String>() {
            public String execute(String input) {
                return input;
            }
        };

        StringWriter writer = new StringWriter();
        new JungRichNewickPrinterCompact<String>().print(graph.Graph, getLabel, writer);
        System.out.println(writer.toString());

    }
}
