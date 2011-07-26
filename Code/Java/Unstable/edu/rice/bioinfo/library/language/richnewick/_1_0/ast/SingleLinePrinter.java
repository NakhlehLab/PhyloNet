package edu.rice.bioinfo.library.language.richnewick._1_0.ast;

import edu.rice.bioinfo.library.programming.extensions.java.lang.iterable.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/19/11
 * Time: 2:09 PM
 * To change this template use File | Settings | File Templates.
 */
public class SingleLinePrinter
{
    public static String toString(Network network)
    {
        StringBuffer accum = new StringBuffer();


        appendDescendantsList(network.PrincipleDescendants, accum);
        appendInfo(network.PrincipleInfo, accum);
        accum.append(';');

        return accum.toString();
    }

    private static void appendDescendantsList(DescendantList descendants, StringBuffer accum)
    {
        Object[] subtrees = IterableHelp.toArray(descendants.Subtrees);

        if(subtrees.length == 0)
            return;

        accum.append("(");
        for(int i = 0; i<subtrees.length; i++)
        {
            Subtree subtree = (Subtree) subtrees[i];
            appendDescendantsList(subtree.Descendants, accum);
            appendInfo(subtree.NetworkInfo, accum);

            if(i != subtrees.length -1)
                accum.append(",");
        }
        accum.append(")");



    }

    private static void appendInfo(NetworkInfo info, final StringBuffer accum)
    {
        String labelPart = info.NodeLabel.execute(new NodeLabelAlgo<String, Object, RuntimeException>() {

           public String forNodeLabelNonEmpty(NodeLabelNonEmpty node, Object input) {

               if(node.Label.OriginallyQuoted)
               {
                    return "'" + node.Label.Content.replace("'", "''") + "'";
               }
               else
               {
                   return node.Label.Content.replace(" ", "_");
               }

           }

           public String forNodeLabelEmpty(NodeLabelEmpty node, Object input)  {
               return "";
           }
       }, null);

       String hybridPart = info.HybridNodeQualifier.execute(new HybridNodeQualifierAlgo<String, Object, RuntimeException>() {

           public String forHybridNodeQualifierEmpty(HybridNodeQualifierEmpty qualifier, Object input)  {
              return "";
           }

           public String forHybridNodeQualifierNonEmpty(HybridNodeQualifierNonEmpty qualifier, Object input)  {
               return "#" + qualifier.HybridNodeIndex.Content;
           }

           public String forHybridNodeQualifierWithType(HybridNodeQualifierWithType qualifier, Object input) throws RuntimeException {
               return "#" + qualifier.HybridNodeType.Content + qualifier.HybridNodeIndex.Content;
           }
       }, null);

        final String branchLengthPart = info.BranchLength.execute(new BranchLengthAlgo<String, Object, RuntimeException>() {

            public String forBranchLengthEmpty(BranchLengthEmpty branchLength, Object input) {
                return "";
            }

            public String forBranchLengthNonEmpty(BranchLengthNonEmpty branchLength, Object input) {
                return ":" + branchLength.Length.Content;
            }
        }, null);

        final String bootstrapPart = info.Bootstrap.execute(new BootstrapAlgo<String, Object, RuntimeException>() {

            public String forBootstrapNonEmpty(BootstrapNonEmpty bootstrap, Object input) {
               return (branchLengthPart == "" ? "::" : ":") + bootstrap.BootstrapValue.Content;
            }

            public String forBootstrapEmpty(BootstrapEmpty bootstrap, Object input) {
                return "";
            }
        }, null);

        String probabilityPart = info.Probability.execute(new ProbabilityAlgo<String, Object, RuntimeException>() {

            public String forProbabilityEmpty(ProbabilityEmpty prob, Object input)  {
                return "";
            }

            public String forProbabilityNonEmpty(ProbabilityNonEmpty prob, Object input) {

                if(branchLengthPart == "" && bootstrapPart == "")
                    return ":::" + prob.ProbabilityValue.Content;
                else if (bootstrapPart == "")
                    return "::" + prob.ProbabilityValue.Content;
                else
                    return ":" + prob.ProbabilityValue.Content;

            }
        }, null);

        accum.append(labelPart + hybridPart + branchLengthPart + bootstrapPart + probabilityPart);
    }
}
