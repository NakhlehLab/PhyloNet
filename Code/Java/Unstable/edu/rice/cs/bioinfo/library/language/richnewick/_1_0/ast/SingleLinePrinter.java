package edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast;

import edu.rice.cs.bioinfo.library.programming.Func;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.*;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/19/11
 * Time: 2:09 PM
 * To change this template use File | Settings | File Templates.
 */
public class SingleLinePrinter
{
    public static String toString(NetworkNonEmpty network)
    {
        return toString(network, new Func1<String, String>() {
            public String execute(String s) {
                return s;
            }
        });
    }

    public static String toString(NetworkNonEmpty network, Func1<String,String> supportTransformer)
    {
        StringBuffer accum = new StringBuffer();


        appendDescendantsList(network.PrincipleDescendants, accum, supportTransformer);
        appendInfo(network.PrincipleInfo, accum, supportTransformer);
        accum.append(';');

        return accum.toString();
    }

    private static void appendDescendantsList(DescendantList descendants, StringBuffer accum, Func1<String,String> supportTransformer)
    {
        Object[] subtrees = IterableHelp.toArray(descendants.Subtrees);

        if(subtrees.length == 0)
            return;

        accum.append("(");
        for(int i = 0; i<subtrees.length; i++)
        {
            Subtree subtree = (Subtree) subtrees[i];
            appendDescendantsList(subtree.Descendants, accum, supportTransformer);
            appendInfo(subtree.NetworkInfo, accum, supportTransformer);

            if(i != subtrees.length -1)
                accum.append(",");
        }
        accum.append(")");



    }

    private static void appendInfo(NetworkInfo info, final StringBuffer accum, final Func1<String,String> supportTransformer)
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

        final String supportPart = info.Support.execute(new SupportAlgo<String, Object, RuntimeException>() {

            public String forSupportNonEmpty(SupportNonEmpty support, Object input) {
               return (branchLengthPart == "" ? "::" : ":") + supportTransformer.execute(support.SupportValue.Content);
            }

            public String forSupportEmpty(SupportEmpty support, Object input) {
                return "";
            }
        }, null);

        String probabilityPart = info.Probability.execute(new ProbabilityAlgo<String, Object, RuntimeException>() {

            public String forProbabilityEmpty(ProbabilityEmpty prob, Object input)  {
                return "";
            }

            public String forProbabilityNonEmpty(ProbabilityNonEmpty prob, Object input) {

                if(branchLengthPart == "" && supportPart == "")
                    return ":::" + prob.ProbabilityValue.Content;
                else if (supportPart == "")
                    return "::" + prob.ProbabilityValue.Content;
                else
                    return ":" + prob.ProbabilityValue.Content;

            }
        }, null);

        accum.append(labelPart + hybridPart + branchLengthPart + supportPart + probabilityPart);
    }
}
