package edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork;

import edu.rice.cs.bioinfo.library.programming.Tuple3;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 11/21/18
 * Time: 2:50 PM
 * To change this template use File | Settings | File Templates.
 */
public class SNProblem {
    protected List<Tuple3<Network, String, Double>> subnetworks = new ArrayList<>();
    protected Network trueNetwork = null;

    public void AddSubNetwork(Network network, String filename, double proportion) {
        subnetworks.add(new Tuple3<>(network, filename, proportion));
    }

    public void SetTrueNetwork(Network trueNetwork) {
        this.trueNetwork = trueNetwork.clone();
    }

    public int GetSize() {
        return subnetworks.size();
    }

    /*
     * @Author: Zhen Cao
     * @Date: 2019-10-22
     */
    public List<Tuple3<Network, String, Double>> getSubnetworks(){
        return subnetworks;
    }

//    @Override
//    public String toString() {
//        StringBuilder sb = new StringBuilder();
//        for (int i = 0; i < subnetworks.size(); i++) {
//            sb.append(subnetworks.get(i));
//            sb.append("\n");
//        }
//        return sb.toString();
//    }
}
