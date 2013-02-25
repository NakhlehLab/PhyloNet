package edu.rice.cs.bioinfo.library.phylogenetics.phylographfactories.richnewick._1_0.ast;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.*;
import edu.rice.cs.bioinfo.library.phylogenetics.*;
import edu.rice.cs.bioinfo.library.phylogenetics.phylograph2factories.richnewick._1_0.ast.RNNode;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/1/12
 * Time: 1:44 PM
 * To change this template use File | Settings | File Templates.
 */
public class NetworkToPhyloGraphFactoryDefault implements GraphFactory<RNNode, PhyloEdge<RNNode>, Network>
{

    public PhyloGraph<RNNode> make(Network input) {

        return input.execute(new NetworkAlgo<PhyloGraph<RNNode>, Object, RuntimeException>() {
            public PhyloGraph<RNNode> forNetworkEmpty(NetworkEmpty network, Object input) throws RuntimeException {
                return new DirectedPhyloGraphDefault<RNNode>();
            }

            public PhyloGraph<RNNode> forNetworkNonEmpty(NetworkNonEmpty network, Object input) throws RuntimeException {

                PhyloGraphBuilder builder;
                if(network.RootageQualifier.execute(new IsRooted(), null))
                {
                    builder = new PhyloGraphBuilder(new DirectedPhyloGraphDefault<RNNode>());
                }
                else
                {
                    builder = new PhyloGraphBuilder(new UndirectedPhyloGraphDefault<RNNode>());
                }

                DAGFactory.makeDAG(network, builder);
                return builder.Graph;

            }
        }, null);
    }
}
