package edu.rice.cs.bioinfo.library.language.tgf_2013_2.writing;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/6/13
 * Time: 1:50 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class TGFWriter<N,E>
{
    public String toTGF(Set<E> edges)
    {
        HashSet<N> nodesAccum = new HashSet<N>();

        for(E edge : edges)
        {
            nodesAccum.add(getSource(edge));
            nodesAccum.add(getDestination(edge));
        }

        return toTGF(nodesAccum, edges);
    }

    public String toTGF(Set<N> nodes, Set<E> edges)
    {
        StringBuffer tfgAccum = new StringBuffer();

        HashMap<N,String> nodeToDotIdent = new HashMap<N, String>();

        for(N node : nodes)
        {
            String nodeIdent = (nodeToDotIdent.size() + 1) + "";
            nodeToDotIdent.put(node, nodeIdent );
            String nodeLabel = getNodeLabel(node);

            tfgAccum.append(nodeIdent);

            if(nodeLabel != null)
            {
                tfgAccum.append(" " + nodeLabel);
            }

            tfgAccum.append("\n");

        }

        tfgAccum.append("#\n");

        for(E edge : edges)
        {
            N source = getSource(edge);
            String sourceIdent = nodeToDotIdent.get(source);

            N dest = getDestination(edge);
            String destIdent = nodeToDotIdent.get(dest);

            tfgAccum.append(sourceIdent + " " + destIdent);

            String edgeLabel = getEdgeLabel(edge);

            if(edgeLabel != null)
            {
                tfgAccum.append(" " + edgeLabel);
            }

            tfgAccum.append("\n");
        }

        return tfgAccum.toString();
    }

    protected String getEdgeLabel(E edge) {
        return null;
    }

    protected String getNodeLabel(N node)
    {
        return null;
    }

    protected abstract N getSource(E edge);

    protected abstract N getDestination(E edge);
}
