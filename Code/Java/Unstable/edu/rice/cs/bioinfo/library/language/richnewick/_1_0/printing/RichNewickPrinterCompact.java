package edu.rice.cs.bioinfo.library.language.richnewick._1_0.printing;

//import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.Probability;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.Ref;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;

import javax.print.attribute.standard.MediaSize;
import javax.xml.soap.Node;
import java.awt.*;
import java.io.BufferedWriter;
import java.io.StringWriter;
import java.io.Writer;
import java.util.HashSet;
import java.util.LinkedList;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/14/12
 * Time: 4:12 PM
 * To change this template use File | Settings | File Templates.
 */
public class RichNewickPrinterCompact<N> implements RichNewickPrinter<N>
{
    public void print(boolean isRooted, N printRoot, Func1<N, String> getLabel, Func1<N, Iterable<N>> getDestinationNodes, Func2<N,N,String> getBranchLength, Func2<N,N,String> getSupport,
                             Func2<N,N,String> getProbability, Func1<N, String> getHybridIndex, Func1<N,HybridNodeType> getHybridNodeType, StringWriter writer)
    {
       StringBuffer buffer = writer.getBuffer();
       buffer.append(";");

        HashSet<String> seenHybridIndexes = new HashSet<String>();
       if(!isRooted && IterableHelp.count(getDestinationNodes.execute(printRoot)) == 1)
       {
           buffer.insert(0, ")");
           N neighbor = getDestinationNodes.execute(printRoot).iterator().next();
           prependNode(getLabel, getDestinationNodes, getBranchLength, getSupport, getProbability, getHybridIndex, getHybridNodeType, buffer, printRoot, neighbor, seenHybridIndexes);
           buffer.insert(0, ",");
           prependNode(getLabel, getDestinationNodes, getBranchLength, getSupport, getProbability, getHybridIndex, getHybridNodeType, buffer, neighbor, printRoot, seenHybridIndexes);
           buffer.insert(0, "(");

       }
       else
       {
            prependNode(getLabel, getDestinationNodes, getBranchLength, getSupport, getProbability, getHybridIndex, getHybridNodeType, buffer, printRoot, null, seenHybridIndexes);
       }

       if(!isRooted)
       {
            buffer.insert(0, "[&U]");
       }
    }

    private void prependNode(Func1<N, String> getLabel, Func1<N, Iterable<N>> getDestinationNodes, Func2<N,N,String> getBranchLength, Func2<N,N,String> getSupport,
                             Func2<N,N,String> getProbability, Func1<N, String> getHybridIndex, Func1<N,HybridNodeType> getHybridNodeType, StringBuffer buffer, N node, N printParent,
                             HashSet<String> seenHybridIndexes)
    {
        Ref<Boolean> seenHybridNode = new Ref<Boolean>(false);
        String nodeInfo = makeNodeInfoString(node, printParent, getLabel, getBranchLength, getSupport, getProbability, getHybridIndex, getHybridNodeType, seenHybridIndexes,
                seenHybridNode);
        buffer.insert(0, nodeInfo);

        if(!seenHybridNode.get())
        {
            LinkedList<N> nonParentDest = new LinkedList<N>();

            for(N dest : getDestinationNodes.execute(node))
            {
                if(printParent == null || !printParent.equals(dest))
                {
                    nonParentDest.add(dest);
                }
            }

            if(nonParentDest.size() > 0)
            {
                buffer.insert(0, ")");

                for(int i = 0; i<nonParentDest.size(); i++)
                {
                    prependNode(getLabel, getDestinationNodes, getBranchLength, getSupport, getProbability, getHybridIndex, getHybridNodeType, buffer, nonParentDest.get(i), node,
                                seenHybridIndexes);

                     if(i<nonParentDest.size() - 1)
                    {
                        buffer.insert(0, ",");
                    }

                }


                buffer.insert(0, "(");
            }
        }


    }

    private String makeNodeInfoString(N node, N printParent, Func1<N, String> getLabel, Func2<N,N,String> getBranchLength, Func2<N,N,String> getSupport, Func2<N,N,String> getProbability,
                                      Func1<N, String> getHybridIndex, Func1<N, HybridNodeType> getHybridNodeType, HashSet<String> seenHybridIndexes, Ref<Boolean> seenHybridNode) {
        String result = getLabel.execute(node);

        String hybridIndex = getHybridIndex.execute(node);

        if(hybridIndex != null)
        {
            result+="#";

            HybridNodeType hybridType = getHybridNodeType.execute(node);

            if(hybridType == HybridNodeType.Hybridization)
            {
                result+="H";
            }
            else if(hybridType == HybridNodeType.Recombination)
            {
               result+="R";
            }
            else if(hybridType == HybridNodeType.LateralGeneTransfer)
            {
                result+="LGT";
            }
            result+= hybridIndex;

            if(seenHybridIndexes.contains(hybridIndex))
            {
                seenHybridNode.set(true);
            }
            else
            {
                seenHybridNode.set(false);
                seenHybridIndexes.add(hybridIndex);
            }
        }

        if(printParent != null)
        {

            String branchLength = getBranchLength.execute(printParent, node);
            String support = getSupport.execute(printParent, node);
            String probability = getProbability.execute(printParent, node);

            if(branchLength != null)
            {
                result += ":" + branchLength;
            }

            if(support != null)
            {
                if(branchLength != null)
                {
                    result += ":" + support;
                }
                else
                {
                    result += "::" + support;
                }
            }

            if(probability != null)
            {
                if(support != null)
                {
                    result += ":" + probability;
                }
                else if(branchLength != null)
                {
                    result += "::" + probability;
                }
                else
                {
                    result += ":::" + probability;
                }
            }
        }

        return result;
    }
}
