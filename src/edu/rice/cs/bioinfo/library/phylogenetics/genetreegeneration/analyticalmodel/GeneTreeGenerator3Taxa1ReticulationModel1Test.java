package edu.rice.cs.bioinfo.library.phylogenetics.genetreegeneration.analyticalmodel;

import edu.rice.cs.bioinfo.library.phylogenetics.DirectedPhyloGraphDefault;
import edu.rice.cs.bioinfo.library.phylogenetics.Graph;
import edu.rice.cs.bioinfo.library.phylogenetics.PhyloEdge;
import org.junit.Test;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 12/4/12
 * Time: 2:37 PM
 * To change this template use File | Settings | File Templates.
 */
public class GeneTreeGenerator3Taxa1ReticulationModel1Test
{
    @Test
    public void testCanGenerateForNetwork1()
    {
        Graph<String, PhyloEdge<String>> network = new DirectedPhyloGraphDefault();
        network.addNode("R");
        network.addNode("D");
        network.addNode("F");

    }
}
