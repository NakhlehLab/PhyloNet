package edu.rice.cs.bioinfo.library.language.dot_2013_1.printing;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/5/13
 * Time: 2:13 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class DotPrinterBase<N,E>
{
    public String toDot(Set<E> edges, String graphType)
    {
        HashSet<N> nodesAccum = new HashSet<N>();

        for(E edge : edges)
        {
            nodesAccum.add(getEdgeLhs(edge));
            nodesAccum.add(getEdgeRhs(edge));
        }

        return toDot(nodesAccum, edges, graphType);
    }

    protected String toDot(Set<N> nodes, Set<E> edges, String graphType)
    {
        StringBuffer dotAccum = new StringBuffer(graphType + " g {\n");

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
            N lhs = getEdgeLhs(edge);
            N rhs = getEdgeRhs(edge);
            String edgeString = getEdgeString();

            String sourceIdent = nodeToDotIdent.get(lhs);
            String destIdent = nodeToDotIdent.get(rhs);

            dotAccum.append(sourceIdent + edgeString + destIdent);

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

    protected abstract String getEdgeString();

    protected abstract N getEdgeRhs(E edge);

    protected abstract N getEdgeLhs(E edge);

    protected String getEdgeLabel(E edge) {
        return null;
    }

    protected String getNodeLabel(N node)
    {
        return null;
    }
}
