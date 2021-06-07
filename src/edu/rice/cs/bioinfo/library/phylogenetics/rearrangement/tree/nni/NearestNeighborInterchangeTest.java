package edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.tree.nni;

import edu.rice.cs.bioinfo.library.phylogenetics.Graph;
import edu.rice.cs.bioinfo.library.phylogenetics.IsLeaf;
import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.Proc4;
import edu.rice.cs.bioinfo.library.programming.Ref;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import org.junit.Assert;
import org.junit.Test;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/6/12
 * Time: 3:32 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class NearestNeighborInterchangeTest<T extends Graph<String, Tuple<String,String>>>
{
    private Func2<String,String,Tuple<String,String>> _makeEdge;

    public NearestNeighborInterchangeTest(Func2<String,String,Tuple<String,String>> makeEdge)
    {
        _makeEdge = makeEdge;
    }

    @Test
    public void testComputeRearrangementsWithoutValidation1()
    {
        T tree = makeUnrootedTree("T", "J", "A", "B", "C", "D");

        tree.addEdge(new Tuple<String, String>("T", "J"));
        tree.addEdge(new Tuple<String, String>("T", "A"));
        tree.addEdge(new Tuple<String, String>("T", "B"));
        tree.addEdge(new Tuple<String, String>("J", "C"));
        tree.addEdge(new Tuple<String, String>("J", "D"));

        final Ref<Boolean> case1Covered = new Ref<Boolean>(false);
        final Ref<Boolean> case2Covered = new Ref<Boolean>(false);
        final Ref<Integer> numCases = new Ref<Integer>(0);
        final IsLeaf isLeaf = new IsLeaf();
        new NearestNeighborInterchangeInPlace(_makeEdge).computeRearrangementsWithoutValidation(tree, new Proc4<T,Tuple<String,String>,Tuple<String,String>,Tuple<String,String>>()
        {
            public void execute(T tree, Tuple<String,String> internalEdge, Tuple<String,String> swapEdgeA, Tuple<String,String> swapEdgeB)
            {
                Assert.assertTrue(isLeaf.execute(tree, "A"));
                Assert.assertTrue(isLeaf.execute(tree, "B"));
                Assert.assertTrue(isLeaf.execute(tree, "C"));
                Assert.assertTrue(isLeaf.execute(tree, "D"));
                Assert.assertTrue(!isLeaf.execute(tree, "T"));
                Assert.assertTrue(!isLeaf.execute(tree, "J"));
                Assert.assertTrue(containsEdge(tree, "T", "J", false));
                Assert.assertTrue(IterableHelp.countInt(tree.getIncidentEdges("T")) == 3);
                Assert.assertTrue(IterableHelp.countInt(tree.getIncidentEdges("J")) == 3);

                if(tree.getIncidentEdges("A").iterator().next().other("A").equals(tree.getIncidentEdges("D").iterator().next().other("D")))
                {
                    Assert.assertTrue(tree.getIncidentEdges("B").iterator().next().other("B") == tree.getIncidentEdges("C").iterator().next().other("C"));
                    case1Covered.set(true);
                }
                else if(tree.getIncidentEdges("A").iterator().next().other("A").equals(tree.getIncidentEdges("C").iterator().next().other("C")))
                {
                    Assert.assertTrue(tree.getIncidentEdges("B").iterator().next().other("B") == tree.getIncidentEdges("D").iterator().next().other("D"));
                    case2Covered.set(true);
                }

                numCases.set(numCases.get() + 1);

            }
        });

        Assert.assertTrue(case1Covered.get());
        Assert.assertTrue(case2Covered.get());
        Assert.assertTrue(numCases.get() == 2);
    }

    @Test
    public void testComputeRearrangementsWithoutValidation2()
    {
        T tree = makeRootedTree("T", "J", "A", "B", "C", "D");

        tree.addEdge(new Tuple<String, String>("T", "J"));
        tree.addEdge(new Tuple<String, String>("T", "A"));
        tree.addEdge(new Tuple<String, String>("T", "B"));
        tree.addEdge(new Tuple<String, String>("J", "C"));
        tree.addEdge(new Tuple<String, String>("J", "D"));

        final Ref<Boolean> case1Covered = new Ref<Boolean>(false);
        final Ref<Boolean> case2Covered = new Ref<Boolean>(false);
        final Ref<Integer> numCases = new Ref<Integer>(0);
        final IsLeaf isLeaf = new IsLeaf();
        new NearestNeighborInterchangeInPlace(_makeEdge).computeRearrangementsWithoutValidation(tree, new Proc4<T,Tuple<String,String>,Tuple<String,String>,Tuple<String,String>>()
        {
            public void execute(T tree, Tuple<String,String> internalEdge, Tuple<String,String> swapEdgeA, Tuple<String,String> swapEdgeB)
            {
                Assert.assertTrue(isLeaf.execute(tree, "A"));
                Assert.assertTrue(isLeaf.execute(tree, "B"));
                Assert.assertTrue(isLeaf.execute(tree, "C"));
                Assert.assertTrue(isLeaf.execute(tree, "D"));
                Assert.assertTrue(!isLeaf.execute(tree, "T"));
                Assert.assertTrue(!isLeaf.execute(tree, "J"));
                Assert.assertTrue(containsEdge(tree, "T", "J", true));
                Assert.assertTrue(IterableHelp.countInt(tree.getIncidentEdges("T")) == 3);
                Assert.assertTrue(IterableHelp.countInt(tree.getIncidentEdges("J")) == 3);

                if(tree.getIncidentEdges("A").iterator().next().other("A").equals(tree.getIncidentEdges("D").iterator().next().other("D")))
                {
                    Assert.assertTrue(tree.getIncidentEdges("B").iterator().next().other("B") == tree.getIncidentEdges("C").iterator().next().other("C"));
                    case1Covered.set(true);
                }
                else if(tree.getIncidentEdges("A").iterator().next().other("A").equals(tree.getIncidentEdges("C").iterator().next().other("C")))
                {
                    Assert.assertTrue(tree.getIncidentEdges("B").iterator().next().other("B") == tree.getIncidentEdges("D").iterator().next().other("D"));
                    case2Covered.set(true);
                }

                numCases.set(numCases.get() + 1);

            }
        });

        Assert.assertTrue(case1Covered.get());
        Assert.assertTrue(case2Covered.get());
        Assert.assertTrue(numCases.get() == 2);
    }

    @Test
    public void testComputeRearrangementsWithoutValidation3()
    {
        T tree = makeRootedTree("T", "J", "A", "B", "C", "D");

        tree.addEdge(new Tuple<String, String>("A", "X"));
        tree.addEdge(new Tuple<String, String>("A", "Z"));
        tree.addEdge(new Tuple<String, String>("A", "T"));
        tree.addEdge(new Tuple<String, String>("T", "J"));
        tree.addEdge(new Tuple<String, String>("T", "B"));
        tree.addEdge(new Tuple<String, String>("J", "C"));
        tree.addEdge(new Tuple<String, String>("J", "D"));

        final Ref<Boolean> case1Covered = new Ref<Boolean>(false);
        final Ref<Boolean> case2Covered = new Ref<Boolean>(false);
        final Ref<Boolean> case3Covered = new Ref<Boolean>(false);
        final Ref<Boolean> case4Covered = new Ref<Boolean>(false);
        final Ref<Integer> numCases = new Ref<Integer>(0);
        final IsLeaf isLeaf = new IsLeaf();
        new NearestNeighborInterchangeInPlace(_makeEdge).computeRearrangementsWithoutValidation(tree, new Proc4<T,Tuple<String,String>,Tuple<String,String>,Tuple<String,String>>()
        {
            public void execute(T tree, Tuple<String,String> internalEdge, Tuple<String,String> swapEdgeA, Tuple<String,String> swapEdgeB)
            {
                Assert.assertTrue(isLeaf.execute(tree, "X"));
                Assert.assertTrue(isLeaf.execute(tree, "Z"));
                Assert.assertTrue(isLeaf.execute(tree, "B"));
                Assert.assertTrue(isLeaf.execute(tree, "C"));
                Assert.assertTrue(isLeaf.execute(tree, "D"));
                Assert.assertTrue(!isLeaf.execute(tree, "A"));
                Assert.assertTrue(!isLeaf.execute(tree, "T"));
                Assert.assertTrue(!isLeaf.execute(tree, "J"));
                Assert.assertTrue(containsEdge(tree, "A", "T", true));
                Assert.assertTrue(IterableHelp.countInt(tree.getIncidentEdges("A")) == 3);
                Assert.assertTrue(IterableHelp.countInt(tree.getIncidentEdges("T")) == 3);
                Assert.assertTrue(IterableHelp.countInt(tree.getIncidentEdges("J")) == 3);

                if(tree.getIncidentEdges("A").iterator().next().other("A").equals(tree.getIncidentEdges("D").iterator().next().other("D")))
                {
                    Assert.assertTrue(tree.getIncidentEdges("B").iterator().next().other("B") == tree.getIncidentEdges("C").iterator().next().other("C"));
                    case1Covered.set(true);
                }
                else if(tree.getIncidentEdges("A").iterator().next().other("A").equals(tree.getIncidentEdges("C").iterator().next().other("C")))
                {
                    Assert.assertTrue(tree.getIncidentEdges("B").iterator().next().other("B") == tree.getIncidentEdges("D").iterator().next().other("D"));
                    case2Covered.set(true);
                }
                else if(tree.getIncidentEdges("Z").iterator().next().other("Z").equals("A") &&
                        tree.getIncidentEdges("B").iterator().next().other("B").equals("A"))
                {
                    Assert.assertTrue(tree.getIncidentEdges("X").iterator().next().other("X").equals("T"));
                    Assert.assertTrue(tree.getIncidentEdges("C").iterator().next().other("C").equals("J"));
                    Assert.assertTrue(tree.getIncidentEdges("D").iterator().next().other("D").equals("J"));
                    case3Covered.set(true);
                }
                else if(tree.getIncidentEdges("X").iterator().next().other("X").equals("A") &&
                        tree.getIncidentEdges("B").iterator().next().other("B").equals("A"))
                {
                    Assert.assertTrue(tree.getIncidentEdges("Z").iterator().next().other("Z").equals("T"));
                    Assert.assertTrue(tree.getIncidentEdges("C").iterator().next().other("C").equals("J"));
                    Assert.assertTrue(tree.getIncidentEdges("D").iterator().next().other("D").equals("J"));
                    case4Covered.set(true);
                }


                numCases.set(numCases.get() + 1);

            }
        });

        Assert.assertTrue(case1Covered.get());
        Assert.assertTrue(case2Covered.get());
        Assert.assertTrue(case3Covered.get());
        Assert.assertTrue(case4Covered.get());
        Assert.assertTrue(numCases.get() == 4);
    }

    protected abstract T makeUnrootedTree(String... nodes);

    protected abstract T makeRootedTree(String... nodes);

    protected  abstract boolean containsEdge(T tree, String source, String destination, boolean directed);
}
