package edu.rice.bioinfo.library.language.richnewick._1_0.ast;

import org.junit.*;

import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 7/19/11
 * Time: 5:33 PM
 * To change this template use File | Settings | File Templates.
 */
public class SingleLinePrinterTest {

    @Test
    public void testToString() {
        NetworkInfo a = new NetworkInfo(new NodeLabelNonEmpty(new Text("A", 1, 0, false)), HybridNodeQualifierEmpty.Singleton,
                BranchLengthEmpty.Singleton, SupportEmpty.Singleton, ProbabilityEmpty.Singleton);
        NetworkInfo a1 = new NetworkInfo(new NodeLabelNonEmpty(new Text("A1", 1, 0, false)), HybridNodeQualifierEmpty.Singleton,
                BranchLengthEmpty.Singleton, SupportEmpty.Singleton, ProbabilityEmpty.Singleton);
        NetworkInfo a2 = new NetworkInfo(new NodeLabelNonEmpty(new Text("A2", 1, 0, false)), HybridNodeQualifierEmpty.Singleton,
                BranchLengthEmpty.Singleton, SupportEmpty.Singleton, ProbabilityEmpty.Singleton);

        NetworkInfo b = new NetworkInfo(new NodeLabelNonEmpty(new Text("B", 1, 0, false)), HybridNodeQualifierEmpty.Singleton,
                BranchLengthEmpty.Singleton, SupportEmpty.Singleton, ProbabilityEmpty.Singleton);
        NetworkInfo c = new NetworkInfo(new NodeLabelNonEmpty(new Text("C", 1, 0, false)), HybridNodeQualifierEmpty.Singleton,
                BranchLengthEmpty.Singleton, SupportEmpty.Singleton, ProbabilityEmpty.Singleton);
        NetworkInfo r = new NetworkInfo(new NodeLabelNonEmpty(new Text("R", 1, 0, false)), HybridNodeQualifierEmpty.Singleton,
                BranchLengthEmpty.Singleton, SupportEmpty.Singleton, ProbabilityEmpty.Singleton);

        ArrayList<Subtree> aChildren = new ArrayList<Subtree>();
        aChildren.add(new Subtree(DescendantList.EMPTY_DESCENDANT_LIST, a1));
        aChildren.add(new Subtree(DescendantList.EMPTY_DESCENDANT_LIST, a2));

        ArrayList<Subtree> networkDl = new ArrayList<Subtree>();
        networkDl.add(new Subtree(new DescendantList(aChildren), a));
        networkDl.add(new Subtree(DescendantList.EMPTY_DESCENDANT_LIST, b));
        networkDl.add(new Subtree(DescendantList.EMPTY_DESCENDANT_LIST, c));


        Network network = new Network(RootageQualifierEmpty.Singleton, new DescendantList(networkDl), r);

        Assert.assertEquals("((A1,A2)A,B,C)R;", SingleLinePrinter.toString(network));

    }
}
