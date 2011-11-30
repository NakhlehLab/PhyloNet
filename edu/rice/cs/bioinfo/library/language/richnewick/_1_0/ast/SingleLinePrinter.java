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
    private boolean _excludeProbability = false;

    public void setExcludeProbability(boolean value)
    {
        _excludeProbability = value;
    }

    private Func1<String,String> _supportTransformer = null;

    public void setSupportTransformer(Func1<String,String> value)
    {
        _supportTransformer = value;
    }

    private Func1<String,String> _hybridQualifierTransformer;

    public void setHybridQualifierTransformer(Func1<String,String> value)
    {
        _hybridQualifierTransformer = value;
    }

    private Func1<String,String> _nodeLabelTransformer;

    public void setNodeLabelTransformer(Func1<String,String> value)
    {
        _nodeLabelTransformer = value;
    }


    public String toString(NetworkNonEmpty network)
    {

        StringBuffer accum = new StringBuffer();


        appendDescendantsList(network.PrincipleDescendants, accum);
        appendInfo(network.PrincipleInfo, accum);
        accum.append(';');

        return accum.toString();
    }

    private void appendDescendantsList(DescendantList descendants, StringBuffer accum)
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

    private void appendInfo(NetworkInfo info, final StringBuffer accum)
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
       labelPart = _nodeLabelTransformer != null ? _nodeLabelTransformer.execute(labelPart) : labelPart;

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
       hybridPart = _hybridQualifierTransformer != null ? _hybridQualifierTransformer.execute(hybridPart) : hybridPart;

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
                String prefix =  (branchLengthPart == "" ? "::" : ":");
                String suffix =  _supportTransformer == null ? support.SupportValue.Content : _supportTransformer.execute(support.SupportValue.Content);
               return  prefix + suffix;
            }

            public String forSupportEmpty(SupportEmpty support, Object input) {
                return "";
            }
        }, null);

        String probabilityPart = _excludeProbability ? "" : info.Probability.execute(new ProbabilityAlgo<String, Object, RuntimeException>() {

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
