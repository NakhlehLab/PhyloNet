package edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast;

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


        NetworkNonEmpty network = new NetworkNonEmpty(RootageQualifierEmpty.Singleton, new DescendantList(networkDl), r);

        Assert.assertEquals("((A1,A2)A,B,C)R;", new SingleLinePrinter().toString(network));

         NetworkInfo z = new NetworkInfo(new NodeLabelNonEmpty(new Text("Z", 1, 0, false)), new HybridNodeQualifierNonEmpty(new Text("1", -1, -1, false)),
                BranchLengthEmpty.Singleton, SupportEmpty.Singleton, ProbabilityEmpty.Singleton);

        ArrayList<Subtree> a1Children = new ArrayList<Subtree>();
        a1Children.add(new Subtree(DescendantList.EMPTY_DESCENDANT_LIST, z));

        ArrayList<Subtree> a2Children = new ArrayList<Subtree>();
        a2Children.add(new Subtree(DescendantList.EMPTY_DESCENDANT_LIST, z));

        aChildren = new ArrayList<Subtree>();
        aChildren.add(new Subtree(new DescendantList(a1Children), a1));
        aChildren.add(new Subtree(new DescendantList(a2Children), a2));

        network = new NetworkNonEmpty(RootageQualifierEmpty.Singleton, new DescendantList(aChildren), a);

         Assert.assertEquals("((Z#1)A1,(Z#1)A2)A;", new SingleLinePrinter().toString(network));


    }
}
