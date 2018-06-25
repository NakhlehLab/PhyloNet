package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.distribution;

import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.core.NetworkPrior;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.structs.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.structs.PopulationSize;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCPrior.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;

import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Created by wendingqiao on 2/25/16.
 */
public class SpeciesNetPriorDistribution {

    // topology
    private PoissonDistribution numReticulation;
    // node heights
    private ExponentialDistribution branchLength;
    private ExponentialDistribution nodeHeights;
    private ExponentialDistribution diameters;
    private Double maxTau;
    // pop size
    private PopulationSize popSize;

    private BetaDistribution inheritance;

    // valid test
    private Set<String> taxa = null;

    public SpeciesNetPriorDistribution(double poissonParam) {
        this(poissonParam, new PopulationSize());
    }

    public SpeciesNetPriorDistribution(double poissonParam, PopulationSize ps) {
        this.numReticulation = new PoissonDistribution(poissonParam);
        nodeHeights = Utils._TIMES_EXP_PRIOR ? new ExponentialDistribution(Utils.EXP_PARAM) : null;
        branchLength = Utils._BL_EXP_PRIOR ? new ExponentialDistribution(Utils.EXP_PARAM) : null;
        maxTau = Utils._BL_EXP_PRIOR ? Utils.EXP_CUTOFF : null;
        diameters = Utils._DIAMETER_PRIOR ? new ExponentialDistribution(Utils.EXP_PARAM) : null;
        this.popSize = ps;
        inheritance = new BetaDistribution(2, 2);
    }

    public double logPrior(Network<NetNodeInfo> net) {
        if(taxa == null) {
            taxa = new HashSet<>();
            for(NetNode<NetNodeInfo> node : net.getLeaves()) {
                taxa.add(node.getName());
            }
        }
        double topPrior = Utils.DISABLE_TOPOLOGY_MOVES ? 0.0 : getTopPrior(net);
        double diametersPrior = getDiametersPrior(net);
        double nodeHeightsPrior = getNodeHeightsPrior(net);
        double blPrior = getBLPrior(net);
        double inheritanceProbPrior = getInheritanceProbPrior(net);
        double popSizePrior = getPopSizePrior(net);
        double prior = topPrior + diametersPrior + nodeHeightsPrior
                + blPrior + inheritanceProbPrior
                + popSizePrior;
        return prior;
    }

    private double getTopPrior(Network<NetNodeInfo> net) {
        int numReti = net.getReticulationCount();
        int numEdges = net.getLeafCount() * 2 - 2;
        double networkSizePrior = Utils._NETWORK_SIZE_PRIOR ? getNetworkSizePrior(numReti, numEdges) : 0.0;
        return Math.log(numReticulation.probability(numReti)) + networkSizePrior;
    }

    private double getNetworkSizePrior(int reti, int edges) {
        double logP = 0;
        for(int i = 1; i <= reti; i++) {
            edges += 3;
            logP += - Math.log((double)(edges * (edges - 1)));
        }
        return logP;
    }

    private double getDiametersPrior(Network<NetNodeInfo> network) {
        if(diameters == null || network.getReticulationCount() == 0) return 0;
        Map<NetNode, Double> distMap = NetworkPrior.getDiameterMap(network.toString());
        double logP = 0.0;
        for(NetNode key: distMap.keySet()) {
            logP += Math.log(diameters.density(distMap.get(key)));
        }
        return logP;
    }

    private double getNodeHeightsPrior(Network<NetNodeInfo> net) {
        /* Uniform distribution */
        if (nodeHeights == null) return 0;
        /* Exponential distribution */
        double logP = 0;
        for(NetNode<NetNodeInfo> node : net.bfs()) {
            if(node.isLeaf()) continue;
            logP += Math.log(nodeHeights.density(node.getData().getHeight()));
        }
        return logP;
    }

    private double getBLPrior(Network<NetNodeInfo> net) {
        /* Uniform distribution */
        if (branchLength == null) return 0;
        /* Exponential distribution */
        double logP = 0;
        for(NetNode<NetNodeInfo> node : net.bfs()) {
            for(NetNode<NetNodeInfo> parent : node.getParents()) {
                double bl = node.getParentDistance(parent);
                if(maxTau != null) {
                    if(bl > maxTau) return Double.NEGATIVE_INFINITY;
                    logP += Math.log(branchLength.density(bl) / branchLength.cumulativeProbability(maxTau));
                } else {
                    logP += Math.log(branchLength.density(bl) );
                }
            }
        }
        return logP;
    }

    private double getInheritanceProbPrior(Network<NetNodeInfo> net) {
        /* Uniform distribution */
        //return 0;
        /* Beta distribution */
        double logP = 0;
        for(NetNode<NetNodeInfo> node : net.bfs()) {
            if(node.isNetworkNode()) {
                for(NetNode<NetNodeInfo> parent : node.getParents()) {
                    logP += Math.log(inheritance.density(node.getParentProbability(parent)));
                }
            }
        }
        return logP;
    }

    private double getPopSizePrior(Network<NetNodeInfo> net) {
        if(!Utils._ESTIMATE_POP_SIZE) {
            return 0;
        }
        double rootPopSize = net.getRoot().getRootPopSize();
        double logP = Math.log(popSize.density(rootPopSize));
        if(Double.isNaN(rootPopSize) || rootPopSize <= 0) {
            System.err.println("Wrong root pop size: " + rootPopSize);
        }
        if(Utils._CONST_POP_SIZE) {
            return logP + popSize.logDensity();
        }
        for(NetNode<NetNodeInfo> node : net.dfs()) {
            for(NetNode<NetNodeInfo> par : node.getParents()) {
                double ps = node.getParentSupport(par);
                if(Double.isNaN(ps) || ps <= 0) {
                    System.err.println("Wrong pop size: " + node.getName() + "->" + par.getName());
                }
                logP += Math.log(popSize.density(ps));
            }
        }
        return logP + popSize.logDensity(); // prior + hyper-prior
    }


    public boolean isValid(Network<NetNodeInfo> net) {
        if (Networks.hasCycle(net)) {
            if(Utils.DEBUG_MODE) System.err.println("has cycle");
            return false;
        }
        if (!Networks.isDisconnectedNetwork(net, null)) {
            if(Utils.DEBUG_MODE) System.err.println("disconnected network: " + net.toString());
            return false;
        }
        List<NetNode<NetNodeInfo>> networkLeafNode = IterableHelp.toList(net.getLeaves());
        // test leaves
        int count = 0;
        for (NetNode<NetNodeInfo> leaf : networkLeafNode) {
            if (!taxa.contains(leaf.getName())) {
                if(Utils.DEBUG_MODE) System.err.println(leaf.getName() + " is missing");
                return false;
            }
            count++;
        }
        if (count != taxa.size()) {
            if(Utils.DEBUG_MODE) System.err.println("taxa size doesn't match " + count);
            return false;
        }
        // test gamma
        for (NetNode<NetNodeInfo> node : net.dfs()) {
            if (node.isRoot()) continue;
            double totalProb = 0.0;
            StringBuilder log = new StringBuilder(node.getName() + " : ");
            for (NetNode parent : node.getParents()) {
                // distance
                if(node.getParentDistance(parent) < 0) {
                    if(Utils.DEBUG_MODE) {
                        System.err.println("negative branch length! " + node.getName() + "<-" + parent.getName());
                    }
                    return false;
                }
                // population size
                if(Utils.varyPopSizeAcrossBranches() && (
                        Double.isNaN(node.getParentSupport(parent)) || node.getParentSupport(parent) < 0)) {
                    if(Utils.DEBUG_MODE) {
                        System.err.println("Wrong pop size! " + node.getName() + "<-" + parent.getName());
                    }
                    return false;
                }
                // probability
                totalProb += node.getParentProbability(parent);
                log.append(node.getParentProbability(parent) + ", ");
            }
            if( (node.getChildCount() == 1 && node.getParentCount()==1) || (node.getParentCount()==0) ){
                throw new RuntimeException(node.getName() + " node.parentCount = " + node.getParentCount());
            }
            if (node.isNetworkNode() && Math.abs(totalProb - 1.0) > 0.00001) {
                throw new RuntimeException(log + " gamma != 1.0 " + net.toString());
            }
            if (node.isTreeNode() && !Double.isNaN(totalProb) && totalProb != node.NO_PROBABILITY) {
                throw new RuntimeException("gamma is invalid: " + totalProb + " " + net.toString());
            }
        }
        return true;
    }

}
