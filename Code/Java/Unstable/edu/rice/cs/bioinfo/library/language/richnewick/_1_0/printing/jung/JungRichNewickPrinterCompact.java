package edu.rice.cs.bioinfo.library.language.richnewick._1_0.printing.jung;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.printing.HybridNodeType;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.printing.RichNewickPrinterCompact;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.uci.ics.jung.graph.DirectedGraph;
import edu.uci.ics.jung.graph.Graph;

import java.io.StringWriter;
import java.lang.reflect.Array;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 6/12/12
 * Time: 5:34 PM
 * To change this template use File | Settings | File Templates.
 */
public class JungRichNewickPrinterCompact<V> extends RichNewickPrinterCompact<V>
{

    public <E> void print(final DirectedGraph<V,E> digraph, Func1<V, String> getLabel, StringWriter writer)
    {
        final Map<V,String> hybridNodeToHybridIndex = new HashMap<V, String>();
        V root = null;

        for(V vertex : digraph.getVertices())
        {
            if(digraph.getInEdges(vertex).size() > 1)
            {
                hybridNodeToHybridIndex.put(vertex, hybridNodeToHybridIndex.size() + "");
            }
            else if(digraph.getInEdges(vertex).size() == 0)
            {
                root = vertex;
            }
        }

        Func1<V,String> getHybridNodeIndex = new Func1<V, String>() {
            @Override
            public String execute(V input) {
                return hybridNodeToHybridIndex.get(input);
            }
        };

        Func1<V,HybridNodeType> getHybridNodeType = new Func1<V, HybridNodeType>() {
            @Override
            public HybridNodeType execute(V input) {
                return HybridNodeType.Hybridization;
            }
        };

        Func2<V,V,String> noDetail = new Func2<V, V, String>() {
            @Override
            public String execute(V input1, V input2) {
                return null;
            }
        };

        Func1<V, Iterable<V>> getDestinationNodes = new Func1<V, Iterable<V>>() {
            @Override
            public Iterable<V> execute(V input) {

                LinkedList<V> destinationNodes = new LinkedList<V>();
                for(E outEdge : digraph.getOutEdges(input))
                {
                    destinationNodes.add(digraph.getEndpoints(outEdge).getSecond());
                }

                return destinationNodes;
            }
        };

        this.print(true, root, getLabel, getDestinationNodes, noDetail, noDetail, noDetail, getHybridNodeIndex, getHybridNodeType, writer );


    }
}
