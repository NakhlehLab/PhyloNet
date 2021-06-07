package edu.rice.cs.bioinfo.library.phylogenetics.scoring.network.genetreeprobability;

import edu.rice.cs.bioinfo.library.phylogenetics.DirectedPhyloGraphDefault;
import edu.rice.cs.bioinfo.library.phylogenetics.PhyloEdge;
import edu.rice.cs.bioinfo.library.phylogenetics.PhyloGraph;
import edu.rice.cs.bioinfo.library.programming.Func2;
import org.junit.Test;

import java.util.Arrays;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/28/12
 * Time: 6:19 PM
 * To change this template use File | Settings | File Templates.
 */
public class GeneTreeProbabilityYFTest
{
    @Test
    public void test()
    {
        PhyloGraph<String> network = new DirectedPhyloGraphDefault<String>();
        network.addNode("R");
        network.addNode("J");
        network.addNode("K");
        network.addNode("A");
        network.addNode("X");
        network.addNode("D");
        network.addNode("L");
        network.addNode("B");
        network.addNode("C");
        network.addEdge(new PhyloEdge<String>("R", "J", 1));
        network.addEdge(new PhyloEdge<String>("R", "K", 1));
        network.addEdge(new PhyloEdge<String>("K", "D", 2));
        network.addEdge(new PhyloEdge<String>("J", "A", 2));
        network.addEdge(new PhyloEdge<String>("J", "X", 0).setProbability(.3));
        network.addEdge(new PhyloEdge<String>("K", "X", 0).setProbability(.7));
        network.addEdge(new PhyloEdge<String>("X", "L", 1));
        network.addEdge(new PhyloEdge<String>("L", "B", 1));
        network.addEdge(new PhyloEdge<String>("L", "C", 1));



        PhyloGraph<String> gt = new DirectedPhyloGraphDefault<String>();
        gt.addNode("R");
        gt.addNode("K");
        gt.addNode("J");
        gt.addNode("C");
        gt.addNode("A");
        gt.addNode("B");
        gt.addNode("D");
        gt.addEdge(new PhyloEdge<String>("R", "K"));
        gt.addEdge(new PhyloEdge<String>("R", "C"));
        gt.addEdge(new PhyloEdge<String>("K", "A"));
        gt.addEdge(new PhyloEdge<String>("K", "J"));
        gt.addEdge(new PhyloEdge<String>("J", "B"));
        gt.addEdge(new PhyloEdge<String>("J", "D"));


        Func2<String,PhyloGraph<String>,String> getNodeName = new Func2<String, PhyloGraph<String>, String>() {
            public String execute(String name, PhyloGraph<String> input2) {
                return name;
            }
        };

        Func2<PhyloEdge<String>, PhyloGraph<String>, Double> getProbability = new Func2<PhyloEdge<String>, PhyloGraph<String>, Double>() {
            public Double execute(PhyloEdge<String> edge, PhyloGraph<String> input2) {
                return edge.getProbability();
            }
        };

         Func2<PhyloEdge<String>, PhyloGraph<String>, Double> getDistance = new Func2<PhyloEdge<String>, PhyloGraph<String>, Double>() {
            public Double execute(PhyloEdge<String> edge, PhyloGraph<String> input2) {
                return edge.getBranchLength();
            }
        };

        List<Double> result =
        new GeneTreeProbabilityYF().calculateGTDistribution(network, Arrays.asList(gt), getNodeName, getNodeName, getDistance, getProbability, null);
    }
}
