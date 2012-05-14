package edu.rice.cs.bioinfo.library.phylogenetics.scoring.network.acceptancetesting;

import com.sun.org.apache.bcel.internal.generic.NEW;
import edu.rice.cs.bioinfo.library.phylogenetics.Graph;
import edu.rice.cs.bioinfo.library.phylogenetics.GraphReadOnly;
import edu.rice.cs.bioinfo.library.phylogenetics.PhyloEdge;
import edu.rice.cs.bioinfo.library.phylogenetics.scoring.network.MDCOnNetworkYF;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.Func3;
import edu.rice.cs.bioinfo.library.programming.Func4;
import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import junit.framework.Assert;
import org.junit.Test;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/14/12
 * Time: 1:57 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class MDCOnNetworkYFFromRichNewick<G extends Graph<String,PhyloEdge<String>>>
{
    protected abstract G makeNetwork(String richNewick);

     protected Func1<G, String> makeNode = new Func1<G, String>()
    {
        private int counter = 0;

        public String execute(G input1) {
            return IterableHelp.count(input1.getNodes()) + "";
        }
    };

    protected Func3<G, String, String, PhyloEdge<String>> makeEdge = new Func3<G, String, String, PhyloEdge<String>>()
    {
        public PhyloEdge<String> execute(G arg1, String source, String destination) {
            return new PhyloEdge<String>(source, destination);
        }
    };

    private final Func1<String,String> _getNetworkNodeLabel  = new Func1<String,String>()
    {
        public String execute(String node) {
            return node;
        }
    };

    private final Func2<GraphReadOnly<String,PhyloEdge<String>>, PhyloEdge<String>, Double> _getNetworkDistance = new Func2<GraphReadOnly<String, PhyloEdge<String>>, PhyloEdge<String>, Double>() {
        public Double execute(GraphReadOnly<String, PhyloEdge<String>> network, PhyloEdge<String> edge) {
            return edge.getBranchLength();
        }
    };

    private final Func2<GraphReadOnly<String,PhyloEdge<String>>, PhyloEdge<String>, Double> _getProbability = new Func2<GraphReadOnly<String, PhyloEdge<String>>, PhyloEdge<String>, Double>() {
        public Double execute(GraphReadOnly<String, PhyloEdge<String>> network, PhyloEdge<String> edge) {
            return edge.getBranchLength();
        }
    };

    private final Func4<String,String,Double,Double,PhyloEdge<String>> _makeNetworkEdge = new Func4<String, String, Double, Double, PhyloEdge<String>>() {
        public PhyloEdge<String> execute(String source, String destination, Double distance, Double probability) {
            PhyloEdge<String> edge = new PhyloEdge<String>(source,destination);
            edge.setBranchLength(distance);
            edge.setProbability(probability);
            return edge;
        }
    };

    @Test
    public void test1()
    {
        Graph<String,PhyloEdge<String>> network = makeNetwork("((A:2,((B:1,C:1)K:1)X#1:0::0.3)J:1,(D:2,X#1:0::0.7)L:1)M;");
        Graph<String,PhyloEdge<String>> gt1 = makeNetwork("(C,((B,D),A));");
        Graph<String,PhyloEdge<String>> gt2 = makeNetwork("(B,(D,(C,A)));");
        Graph<String,PhyloEdge<String>> gt3 = makeNetwork("(D,(B,(C,A)));");
        Graph<String,PhyloEdge<String>> gt4 = makeNetwork("((B,D),(C,A));");
        Graph<String,PhyloEdge<String>> gt5 = makeNetwork("((B,D),(C,A));");
        LinkedList<Graph<String,PhyloEdge<String>>> geneTrees = new LinkedList<Graph<String,PhyloEdge<String>>>();
        geneTrees.add(gt1);
        geneTrees.add(gt2);
        geneTrees.add(gt3);
        geneTrees.add(gt4);
        geneTrees.add(gt5);

        MDCOnNetworkYF scorer = new MDCOnNetworkYF();
        List<Integer> result = scorer.countExtraCoal(network, geneTrees, null, _getNetworkNodeLabel, _getNetworkNodeLabel, _getNetworkDistance, _getProbability, _getNetworkDistance, _getProbability,
                              _makeNetworkEdge, _makeNetworkEdge);
        Assert.assertTrue(result.get(0) == 2);
        Assert.assertTrue(result.get(1) == 2);
        Assert.assertTrue(result.get(2) == 1);
        Assert.assertTrue(result.get(3) == 1);
        Assert.assertTrue(result.get(4) == 1);

        int i = 0;
    }

}
