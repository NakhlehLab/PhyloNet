package edu.rice.cs.bioinfo.library.language.richnewick._1_0.printing;

//import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.Probability;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;

import java.io.StringWriter;
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
     public final Func2<N,N,String> NO_DETAIL = new Func2<N, N, String>() {

            public String execute(N input1, N input2) {
                return null;
            }
        };

    private Func2<N,N,String> getBranchLength = NO_DETAIL;

    public void setGetBranchLength(Func2<N,N,String> getBranchLength)
    {
        this.getBranchLength = getBranchLength;
    }

    private Func2<N,N,String> getSupport = NO_DETAIL;

    public void setGetSupport(Func2<N,N,String> getSupport)
    {
        this.getSupport = getSupport;
    }

    private Func2<N,N,String> getProbability = NO_DETAIL;

    public void setGetProbability(Func2<N,N,String> getProbability)
    {
        this.getProbability = getProbability;
    }

    public void print(boolean isRooted, N printRoot, Func1<N, String> getLabel, Func1<N, Iterable<N>> getDestinationNodes, Func1<N, String> getHybridIndex, Func1<N,HybridNodeType> getHybridNodeType, StringWriter writer)
    {
       StringBuffer buffer = writer.getBuffer();
       buffer.append(";");

        HashSet<String> seenHybridIndexes = new HashSet<String>();
       if(!isRooted && IterableHelp.countInt(getDestinationNodes.execute(printRoot)) == 1)
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
        String hybridIndex = getHybridIndex.execute(node);
        boolean isFirstVisitToHybrid =  !seenHybridIndexes.contains(hybridIndex);
        if(hybridIndex != null && isFirstVisitToHybrid)
            seenHybridIndexes.add(hybridIndex);


        String nodeInfo = makeNodeInfoString(node, printParent, getLabel, getBranchLength, getSupport, getProbability, hybridIndex, getHybridNodeType, seenHybridIndexes);
        buffer.insert(0, nodeInfo);

        if(hybridIndex != null && isFirstVisitToHybrid)
            return;


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

    private String makeNodeInfoString(N node, N printParent, Func1<N, String> getLabel, Func2<N,N,String> getBranchLength, Func2<N,N,String> getSupport, Func2<N,N,String> getProbability,
                                      String hybridIndex, Func1<N, HybridNodeType> getHybridNodeType, HashSet<String> seenHybridIndexes) {
        String result = getLabel.execute(node);



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
