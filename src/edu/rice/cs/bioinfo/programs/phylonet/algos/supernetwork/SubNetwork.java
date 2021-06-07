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
        Network trueNetwork = Networks.readNetworkWithRootPop("(((((((P:4.723,(O:2.0)#H1:2.723::0.6):7.87,N:12.593):2.829,(M:6.976,(L:3.0)#H2:3.976::0.8):8.446):3.272,((K:11.378,(#H1:4.0::0.4,J:6.0):5.378):4.977,((H:10.99,((G:6.106,((I:1.0)#H3:1.106::0.7,F:2.106):4.0):3.462,E:9.568):1.422):2.74,(#H2:1.0::0.2,D:4.0):9.73):2.625):2.339):2.535,(#H3:6.229::0.3,C:7.229):14.0):6.136,B:27.365):12.913,A:40.278);");
        System.out.println("True network: " + trueNetwork.toString());

        superNetwork.setTrueNetwork(trueNetwork);
        List<String> leaves = new ArrayList<>();
        leaves.add("C");
        leaves.add("I");
        leaves.add("D");
        leaves.add("H");
        leaves.add("F");
        leaves.add("J");
        leaves.add("L");
        leaves.add("O");
        leaves.add("M");
        leaves.add("P");

        Network subnet = NetworkUtils.getSubNetwork(trueNetwork, leaves, true);
        System.out.println(Networks.getFullString(subnet));
    }
}
