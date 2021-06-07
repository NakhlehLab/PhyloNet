package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.distribution;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.Splitting;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.BiAllelicGTR;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 2/17/19
 * Time: 11:23 PM
 * To change this template use File | Settings | File Templates.
 */
public class SNAPPLikelihoodSampling {

    static public double computeSNAPPLikelihoodST(Network network, List<Splitting> splittings, Map<RPattern, double[]> patterns, BiAllelicGTR BAGTRModel) {
        Algorithms.CORRECTION_AT_LEAVES = true;

        Networks.autoLabelNodes(network);

        Network net = Networks.readNetwork(network.toString());
        net.getRoot().setRootPopSize(network.getRoot().getRootPopSize());

        if(!Utils._ESTIMATE_POP_SIZE) {
            net.getRoot().setRootPopSize(Utils._POP_SIZE_MEAN);
        }

        Double theta = null;
        if(Utils._CONST_POP_SIZE) {
            // theta = network.getRoot().getRootPopSize();
            for(Object childObj : net.bfs()) {
                NetNode child = (NetNode) childObj;
                for(Object parentObj : child.getParents()) {
                    NetNode parent = (NetNode) parentObj;
                    child.setParentSupport(parent, net.getRoot().getRootPopSize());
                }
            }
        }

        String netstring = net.toString();

        double sum = 0.0;
        double sumMono = 0.0;
        double numsites = 0.0;
        int maxLineages = 0;

        for(RPattern pattern : patterns.keySet()) {
            maxLineages = Math.max(maxLineages, pattern.sumLineages());
        }

        if(Algorithms.HAS_DOMINANT_MARKER) {
            maxLineages = maxLineages * 2;
        }

        Network cloneNetwork = Networks.readNetwork(netstring);
        cloneNetwork.getRoot().setRootPopSize(network.getRoot().getRootPopSize());
        SNAPPAlgorithm run = new SNAPPAlgorithm(cloneNetwork, BAGTRModel, theta, maxLineages);
        long start = System.currentTimeMillis();

        for(int i = 0 ; i < splittings.size() ; i++) {
            double likelihood = 0;
            //long start = System.currentTimeMillis();
            try {
                likelihood = run.getProbability(splittings.get(i).getRPattern(), splittings.get(i));
                sum += Math.log(likelihood);
                numsites += 1;
                //System.out.println(represent + " " + likelihood  + " " + count );
            } catch (Exception e) {
                e.printStackTrace();
                System.out.println("Exceptional network" + netstring);
            }
        }

        return sum;
    }
}
