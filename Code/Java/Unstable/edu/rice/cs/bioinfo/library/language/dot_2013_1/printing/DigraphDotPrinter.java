package edu.rice.cs.bioinfo.library.language.dot_2013_1.printing;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/29/13
 * Time: 3:10 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class DigraphDotPrinter<N,E>
{
    public String toDot(Set<E> edges)
    {
        HashSet<N> nodesAccum = new HashSet<N>();

        for(E edge : edges)
        {
            nodesAccum.add(getSource(edge));
            nodesAccum.add(getDestination(edge));
        }

        return toDot(nodesAccum, edges);
    }

    public String toDot(Set<N> nodes, Set<E> edges)
    {
        StringBuffer dotAccum = new StringBuffer("digraph g {\n");

        HashMap<N,String> nodeToDotIdent = new HashMap<N, String>();

        for(N node : nodes)
        {
            String nodeIdent = (nodeToDotIdent.size() + 1) + "";
            nodeToDotIdent.put(node, nodeIdent );
            String nodeLabel = getNodeLabel(node);

            dotAccum.append(nodeIdent);

            if(nodeLabel != null)
            {
                dotAccum.append("[label=\"" + nodeLabel +"\"]");
            }

            dotAccum.append(";\n");

        }


        for(E edge : edges)
        {
            N source = getSource(edge);
            String sourceIdent = nodeToDotIdent.get(source);

            N dest = getDestination(edge);
            String destIdent = nodeToDotIdent.get(dest);

            dotAccum.append(sourceIdent + "->" + destIdent);

            String edgeLabel = getEdgeLabel(edge);

            if(edgeLabel != null)
            {
                dotAccum.append("[label=\"" + edgeLabel +"\"]");
            }

            dotAccum.append(";\n");
        }

        dotAccum.append("}");

        return dotAccum.toString();
    }

    protected String getEdgeLabel(E edge) {
        return null;
    }

    protected String getNodeLabel(N source)
    {
        return null;
    }

    protected abstract N getSource(E edge);

    protected abstract N getDestination(E edge);
}
