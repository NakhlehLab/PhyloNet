package edu.rice.cs.bioinfo.library.language.richnewick._1_0.graphbuilding;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.HybridNodeType;
import netscape.security.PrivilegeTable;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/26/11
 * Time: 1:03 PM
 * To change this template use File | Settings | File Templates.
 */
public class GraphBuilderDOT implements GraphBuilder<Integer>
{
    private int _nextNodeNumber = 1;

    private  StringBuffer _graph = new StringBuffer("digraph network {");

    public Integer createNode(String label)
    {
        try
        {
            String dotLabel = label == null ? "" : label;
            _graph.append("\n\t" + _nextNodeNumber + "[label=\"" + dotLabel + "\"];");
            return new Integer(_nextNodeNumber);
        }
        finally
        {
            _nextNodeNumber++;
        }

    }

    public Integer createHybridNode(String label, HybridNodeType hybridType, BigInteger hybridNodeIndex) {

        try
        {
            String dotLabel = "(" + hybridNodeIndex + (makeHybridTypeLabel(hybridType)) + ")";

            if(label != null)
                dotLabel = label + " " + dotLabel;

            _graph.append("\n\t" + _nextNodeNumber + "[label=\"" + dotLabel + "\" shape=\"box\"];");
            return new Integer(_nextNodeNumber);
        }
        finally
        {
            _nextNodeNumber++;
        }
    }

    private String makeHybridTypeLabel(HybridNodeType hybridType) {

        if(hybridType == null || hybridType == HybridNodeType.Unspecified)
            return "";

        if(hybridType == HybridNodeType.LateralGeneTransfer)
            return "LGT";

        if(hybridType == HybridNodeType.Hybridization)
            return "H";

        if(hybridType == HybridNodeType.Recombination)
            return "R";

         throw new IllegalArgumentException("Unexpected hybrid type.");
    }

    public void createDirectedEdge(Integer tail, Integer tip, BigDecimal branchLength, BigDecimal support, BigDecimal probability) {

        String label = null;

        if(branchLength != null || support != null || probability != null)
        {
            label = "[label=\"";

            ArrayList<String> adj = new ArrayList<String>();

            if(branchLength != null)
                adj.add("bl:" + branchLength.toPlainString());

            if(support != null)
                adj.add("bs:" + support.toPlainString());

            if(probability != null)
                adj.add("p:" + probability);

            for(int i = 0; i<adj.size(); i++)
            {
                label+= adj.get(i).replace('"', '\'');

                if(i != adj.size() -1)
                    label+=" ";
            }


            label+= "\"]";
        }

        _graph.append("\n\t" + tail + "->" + tip + (label == null ? "" : label) + ";");

    }

    public String getDOT()
    {
        _graph.append("\n}");
        return _graph.toString();
    }
}
