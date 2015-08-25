package edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.util;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;

public class NetworkLengthApplier
{
    final private Network net;
    final private Map<String, List<String>> species2alleles;
    final double[] networkParameters;

    private int netParameterIndex = 0;

    public NetworkLengthApplier(Network net, double[] networkParameters, Map<String, List<String>> species2alleles)
    {
        this.net =net;
        this.species2alleles = species2alleles;
        this.networkParameters = networkParameters;
    }

    @SuppressWarnings("unchecked")
    private static <T> Network<T> uglyCast(Network net)
    {
        return net;
    }



    public void apply()
    {
        netParameterIndex = 0;
        Network<Double> foo = uglyCast(net);

        /*
        for(NetNode reticulation: foo.getNetworkNodes()){
            double inheritanceProb = networkParameters[netParameterIndex++];
            for(Object parent: reticulation.getParents()){
                reticulation.setParentProbability((NetNode)parent, inheritanceProb);
                inheritanceProb = 1 - inheritanceProb;
            }
        }
        */

        Map<NetNode, Set<String>> node2leaves = new HashMap<>();
        for(NetNode node: Networks.postTraversal(foo)){
            if(node.isRoot())break;
            Set<String> leaves = new HashSet<>();
            node2leaves.put(node, leaves);
            if(node.isLeaf()){
                leaves.addAll(species2alleles.get(node.getName()));
            }
            else{
                for(Object childO: node.getChildren()){
                    leaves.addAll(node2leaves.get(childO));
                }
            }
            if(leaves.size()>1){
                for(Object parentO: node.getParents()){
                    node.setParentDistance((NetNode)parentO, networkParameters[netParameterIndex++]);
                }
            }
        }

        if (netParameterIndex != networkParameters.length)
            throw new RuntimeException("FAIL adding netlengths");

    }
}
