package edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.rea;

import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.Graph;
import edu.rice.cs.bioinfo.library.programming.*;
import org.junit.*;
import org.junit.internal.matchers.Each;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/9/12
 * Time: 1:28 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class ReticulateEdgeAdditionTest<G extends Graph<String, Tuple<String,String>>>
{
    private Func3<G,String,String, Tuple<String,String>> _makeEdge;

    private Func1<G,String> _makeNode;

    public ReticulateEdgeAdditionTest(Func1<G,String> makeNode, Func3<G,String,String, Tuple<String,String>> makeEdge)
    {
       _makeEdge = makeEdge;
       _makeNode = makeNode;
    }

    protected abstract G makeNetwork(String... nodes);

    protected abstract boolean containsEdge(G network, String source, String destination);

    @Test
    public void testComputeRearrangementsWithoutValidationCase1()
    {
         G network = makeNetwork("R", "A", "B");
         network.addEdge(_makeEdge.execute(network, "R", "A"));
         network.addEdge(_makeEdge.execute(network, "R", "B"));

        final Ref<Integer> numCases = new Ref<Integer>(0);
        final Ref<Boolean> case1Covered = new Ref<Boolean>(false);
        final Ref<Boolean> case2Covered = new Ref<Boolean>(false);
        new ReticulateEdgeAdditionInPlace<G, String, Tuple<String,String>>(_makeNode, _makeEdge).computeRearrangementsWithoutValidation(network,
                new Proc4<G, Tuple<String, String>, Tuple<String, String>, Tuple<String, String>>() {
            public void execute(G network, Tuple<String, String> sourceEdge, Tuple<String, String> destinationEdge, Tuple<String, String> reticulateEdge)
            {
                Assert.assertTrue(containsEdge(network, "R", reticulateEdge.Item1));
                Assert.assertTrue(containsEdge(network, "R", reticulateEdge.Item2));
                Assert.assertFalse(containsEdge(network, "R", "A"));
                Assert.assertFalse(containsEdge(network, "R", "B"));
                Assert.assertTrue(containsEdge(network, reticulateEdge.Item1, reticulateEdge.Item2));
                Assert.assertTrue(containsEdge(network, sourceEdge.Item1, reticulateEdge.Item1));
                Assert.assertTrue(containsEdge(network, destinationEdge.Item1, reticulateEdge.Item2));

                if(sourceEdge.Item2.equals("A"))
                {
                    Assert.assertTrue(containsEdge(network, reticulateEdge.Item1, "A"));
                    Assert.assertTrue(containsEdge(network, reticulateEdge.Item2, "B"));
                    case1Covered.set(true);
                }
                else  if(sourceEdge.Item2.equals("B"))
                {
                    Assert.assertTrue(containsEdge(network, reticulateEdge.Item1, "B"));
                    Assert.assertTrue(containsEdge(network, reticulateEdge.Item2, "A"));
                    case2Covered.set(true);
                }
                numCases.set(numCases.get() + 1);
            }
        });

        Assert.assertTrue(case1Covered.get());
        Assert.assertTrue(case2Covered.get());
        Assert.assertTrue(numCases.get() == 2);

    }
}
