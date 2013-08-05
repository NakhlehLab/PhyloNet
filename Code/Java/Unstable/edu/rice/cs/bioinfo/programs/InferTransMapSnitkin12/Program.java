package edu.rice.cs.bioinfo.programs.infertransmapsnitkin12;

import edu.rice.cs.bioinfo.library.graph.algorithms.spanningTree.FindAMinimumSpanningTreeKruskal;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.SparseMultigraph;
import org.joda.time.Days;
import org.joda.time.LocalDate;
import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import java.io.File;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.util.HashSet;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 7/25/13
 * Time: 1:24 PM
 * To change this template use File | Settings | File Templates.
 */
public class Program
{
    static class GraphNode
    {

        public final String Sequence;

        public final String Name;

        GraphNode(String sequence, String name)
        {
            Sequence = sequence;
            Name = name;
        }
    }

    static class GraphEdge
    {
        public final GraphNode Node1;

        public final GraphNode Node2;

        public final long Weight;

        GraphEdge(GraphNode node1, GraphNode node2)
        {
            Node1 = node1;
            Node2 = node2;
            Weight = getGeneticDistance(node1, node2);
        }
    }

    public static void main(String[] args) throws Exception
    {
        File inFile = new File(args[0]);

        DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
        DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
        Document doc = dBuilder.parse(inFile);



        NodeList elements = doc.getElementsByTagName("Sequencing");

        Graph<GraphNode, GraphEdge> graph = new SparseMultigraph<GraphNode, GraphEdge>();

        for(int i = 0; i<elements.getLength(); i++)
        {
            Node node = elements.item(i);
            String id = node.getAttributes().getNamedItem("sourceId").getNodeValue();
            GraphNode graphNode = new GraphNode(node.getAttributes().getNamedItem("sequence").getNodeValue(),
                                                id.substring("THMS-SA-".length(), id.indexOf(".fastq")));
            graph.addVertex(graphNode);
        }

        for(GraphNode n1 : graph.getVertices())
        {
            for(GraphNode n2 : graph.getVertices())
            {
                if(n1 != n2)
                {
                    graph.addEdge(new GraphEdge(n1, n2), n1, n2);
                }
            }
        }

        Tuple<Set<GraphEdge>,Long> mst = new FindAMinimumSpanningTreeKruskal<Graph<GraphNode, GraphEdge>, GraphEdge, Long>()
        {

            @Override
            protected Long add(Long term1, Long term2)
            {
                return term1 + term2;
            }

            @Override
            protected Long makeZero()
            {
                return 0L;
            }

            @Override
            protected Tuple<?, ?> getNodesOfEdge(GraphEdge edge, Graph<GraphNode, GraphEdge> graph)
            {
                return new Tuple<Object, Object>(edge.Node1, edge.Node2);
            }

            @Override
            protected Long getWeight(GraphEdge edge, Graph<GraphNode, GraphEdge> graph)
            {
                return edge.Weight;
            }

            @Override
            protected Set<? extends GraphEdge> getEdges(Graph<GraphNode, GraphEdge> graph)
            {
                return new HashSet<GraphEdge>(graph.getEdges());
            }

            @Override
            protected Set<?> getNodes(Graph<GraphNode, GraphEdge> graph)
            {
                return new HashSet<GraphNode>(graph.getVertices());
            }
        }.execute(graph);

       StringBuilder graphML = new StringBuilder();
       graphML.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" +
                      "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\"  \n" +
                      "    xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n" +
                      "    xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns\n" +
                      "     http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">\n" +
                      "<key id=\"d1\" for=\"edge\" attr.name=\"weight\" attr.type=\"double\"/>\n" +
                      "<graph id=\"G\" edgedefault=\"undirected\">");

       for(GraphNode node : graph.getVertices())
       {
           graphML.append("\n\t<node id=\"" + node.Name + "\"/>");
       }

        int id = 0;
        for(GraphEdge edge : mst.Item1)
        {
            graphML.append("\n\t<edge id=\"" + (id++) +"\" source=\"" + edge.Node1.Name + "\" target=\"" + edge.Node2.Name +"\">\n" +
                           "\t\t<data key=\"d1\">" + edge.Weight + "</data>\n" +
                           "\t</edge>");
        }

        graphML.append("\n\t</graph>\n" +
                       "</graphml>");

        System.out.println(graphML);


        Class.forName("com.mysql.jdbc.Driver");

        Connection conn = DriverManager.getConnection("jdbc:mysql://localhost:3306/data", "pipeline", "pipeline");
        PreparedStatement select = conn.prepareStatement("SELECT collectionDateUtc FROM isolates WHERE id=?");

        for(GraphEdge edge : mst.Item1)
        {
            String sourceDate = null;
            select.setString(1, "TMHS-SA-" + edge.Node1.Name);
            ResultSet result = select.executeQuery();
            if(result.next())
            {
                sourceDate = result.getString("collectionDateUtc");
            }

            String destDate = null;
            select.setString(1, "TMHS-SA-" + edge.Node2.Name);
            result = select.executeQuery();
            if(result.next())
            {
                destDate = result.getString("collectionDateUtc");
            }

            if(sourceDate != null & destDate != null)
            {
                LocalDate source = LocalDate.parse(sourceDate);
                LocalDate dest = LocalDate.parse(destDate);
                double daysBetween = Math.abs(Days.daysBetween(source, dest).getDays());
                System.out.println(daysBetween + " " + edge.Weight);
            }
        }



    }

    private static Integer getGeneticDistance(GraphNode item1, GraphNode item2)
    {
        int accum = 0;
        if(item1.Sequence.length() != item2.Sequence.length())
            throw new RuntimeException();

        for(int i = 0; i<item1.Sequence.length(); i++)
        {
            if(item1.Sequence.charAt(i) != item2.Sequence.charAt(i))
                accum++;
        }

        return accum;
    }
}
