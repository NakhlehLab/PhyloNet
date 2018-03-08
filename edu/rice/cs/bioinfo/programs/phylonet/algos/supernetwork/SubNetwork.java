package edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 2/17/18
 * Time: 12:02 PM
 * To change this template use File | Settings | File Templates.
 */
public class SubNetwork {

    public static void main(String[] args) {
        SuperNetwork superNetwork = new SuperNetwork(new ArrayList<>());
        Network trueNetwork = Networks.readNetworkWithRootPop("[0.01](((((C:0.0015:0.01)#H1:0.0045:0.01:0.5,B:0.006:0.01):0.003:0.01)#H2:0.0015:0.01:0.5,A:0.0105:0.01):0.0045:0.01,((D:0.003:0.01,#H1:0.0015:0.01:0.5):0.009:0.01,#H2:0.003:0.01:0.5):0.003:0.01);");
        System.out.println("True network: " + trueNetwork.toString());

        superNetwork.setTrueNetwork(trueNetwork);

        superNetwork.genAllSubNetworks(trueNetwork, 3);
        for (Network net : superNetwork.getSubNetworks()) {
            //System.out.println(net.toString());
            System.out.println("-truenet \""+ Networks.getFullString(net) + "\"");
        }
    }
}
