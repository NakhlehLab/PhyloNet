package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.EZONChapter;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.core.State;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.distribution.SNAPPLikelihood;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.alignment.MarkerSeq;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.counting.CoalescenceHistoriesCounting;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimSNPInNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.BiAllelicGTR;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 1/27/18
 * Time: 5:12 PM
 * To change this template use File | Settings | File Templates.
 */

public class FigureRunningTime {
    private static void initNetHeights(Network<NetNodeInfo> network, double popSize) {
        for(NetNode<NetNodeInfo> node : Networks.postTraversal(network)) {
            if(node.isLeaf()) {
                node.setData(new NetNodeInfo(Utils.DEFAULT_NET_LEAF_HEIGHT));
                continue;
            }
            double height = Double.MAX_VALUE;
            for(NetNode<NetNodeInfo> child : node.getChildren()) {
                height = Math.min(height, child.getParentDistance(node) + child.getData().getHeight());
            }
            node.setData(new NetNodeInfo(height));
        }
        boolean setPopSize = true;
        for(NetNode<NetNodeInfo> node : Networks.postTraversal(network)) {
            for(NetNode<NetNodeInfo> par : node.getParents()) {
                node.setParentDistance(par, par.getData().getHeight() - node.getData().getHeight());
                if(node.getParentSupport(par) == node.NO_POP_SIZE && setPopSize) {
                    node.setParentSupport(par, popSize);
                }
            }
        }
        network.getRoot().setRootPopSize(popSize);
    }

    private static void adopt(NetNode<NetNodeInfo> par, NetNode<NetNodeInfo> child, double[] params) {
        par.adoptChild(child, par.getData().getHeight() - child.getData().getHeight());
        child.setParentProbability(par, params[0]);
        child.setParentSupport(par, params[1]);
    }

    private static double[] getParameters(NetNode<NetNodeInfo> par, NetNode<NetNodeInfo> child) {
        return new double[] {child.getParentProbability(par), child.getParentSupport(par)};
    }

    private static void addReticulation(NetNode<NetNodeInfo> v1, NetNode<NetNodeInfo> v2,
                                     NetNode<NetNodeInfo> v3, NetNode<NetNodeInfo> v4,
                                     NetNode<NetNodeInfo> v5, NetNode<NetNodeInfo> v6,
                                     double gamma, double popSize) {

        double[] paramV3V4 = getParameters(v3, v4);
        double[] paramV5V6 = getParameters(v5, v6);

        double a = Math.min(paramV3V4[1], paramV5V6[1]);
        double b = Math.max(paramV3V4[1], paramV5V6[1]);

        v3.removeChild(v4);
        v5.removeChild(v6);

        adopt(v3, v1, new double[] {NetNode.NO_PROBABILITY, popSize} );
        adopt(v1, v4, paramV3V4);
        adopt(v5, v2, new double[] {1.0-gamma, popSize} );
        adopt(v2, v6, paramV5V6);

        adopt(v1, v2, new double[] {gamma, popSize} );
    }

    private static void addRandomReticulation(Network<NetNodeInfo> network) {
        NetNode<NetNodeInfo> _v1, _v2, _v3, _v4, _v5, _v6;
        _v1 = _v2 = _v3 = _v4 = _v5 = _v6 = null; // reset

        List<Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>>> edges = Networks.getAllEdges(network);
        int numEdges = edges.size();
        int numRetiNodes = network.getReticulationCount();

        Tuple<NetNode<NetNodeInfo>, NetNode<NetNodeInfo>> edge1, edge2;
        edge1 = edges.get(Randomizer.getRandomInt(numEdges));
        do {
            edge2 = edges.get(Randomizer.getRandomInt(numEdges));
        } while (edge1 == edge2);

        _v3 = edge1.Item2;
        _v4 = edge1.Item1;
        _v5 = edge2.Item2;
        _v6 = edge2.Item1;

        double t3 = _v3.getData().getHeight();
        double t4 = _v4.getData().getHeight();
        double l1 = t3 - t4;
        double t1 = t4 + 0.5 * l1;

        double t5 = _v5.getData().getHeight();
        double t6 = _v6.getData().getHeight();
        double l2 = t5 - t6;
        double t2 = t6 + 0.5 * l2;

        double gamma = 0.5;

        _v1 = new BniNetNode<>();
        _v1.setData(new NetNodeInfo(t1));
        _v2 = new BniNetNode<>();
        _v2.setData(new NetNodeInfo(t2));

        if(t1 > t2)
            addReticulation(_v1, _v2, _v3, _v4, _v5, _v6, gamma, network.getRoot().getRootPopSize()); // v3 v4 t1 v5 v6 t2
        else
            addReticulation(_v2, _v1, _v5, _v6, _v3, _v4, gamma, network.getRoot().getRootPopSize());
    }

    private static boolean isNetworkValid(Network<NetNodeInfo> network){
        Set<String> leaves = new HashSet<>();
        for(NetNode node : network.getLeaves()) {
            leaves.add(node.getName());
        }

        int count = 0;
        for(Object leaf: network.getLeaves()){
            if(leaves.contains(((NetNode)leaf).getName())){
                count++;
            }
            else{
                return false;
            }
        }
        if(count!=leaves.size())return false;
        if(!Networks.isDisconnectedNetwork(network,null))return false;
        for(Object node: Networks.postTraversal(network)){
            double totalProb = 0;
            for (Object parent : ((NetNode) node).getParents()) {
                totalProb += ((NetNode) node).getParentProbability((NetNode) parent);
            }
            if(((NetNode)node).getChildCount()==1 && ((NetNode)node).getParentCount()<2){
                return false;
            }
            if(totalProb!=NetNode.NO_PROBABILITY && ((NetNode)node).isNetworkNode()){
                if(Math.abs(totalProb - 1) > 0.00001) {
                    throw new RuntimeException(network.toString());
                }
            }
            else if(!((NetNode)node).isRoot()){
                if(totalProb != NetNode.NO_PROBABILITY){
                    throw new RuntimeException(network.toString());
                }
            }
        }
        return true;
    }

    public static void go() {
        Utils._NUM_THREADS = 8;

        int numSites = 10000;
        double pi0 = 0.5;
        double pi1 = 1- pi0;

        BiAllelicGTR BAGTRModel = new BiAllelicGTR(new double[] {pi0, pi1}, new double[] {1.0/ (2.0 * pi0)});
        boolean useOnlyPolymorphic = false;

        List<Network<NetNodeInfo>> starters = new ArrayList<>();
        starters.add(Networks.readNetwork("(1:22.48683929,((2:8.05934715,(((3:1.87509346,4:1.87509346):1.42588234,((5:1.20218468,6:1.20218468):0.04821014,(7:0.07337952,8:0.07337952):1.17701530):2.05058098):0.49372292,9:3.79469872):4.26464844):0.01494217,10:8.07428932):14.41254997);"));
        starters.add(Networks.readNetwork("((1:12.88111496,(2:5.18445969,((3:0.26446915,4:0.26446915):0.76512527,5:1.02959442):4.15486526):7.69665527):1.75988197,((6:1.86362267,7:1.86362267):1.99307442,(8:3.54789925,(9:2.77291870,10:2.77291870):0.77498055):0.30879784):10.78429985);"));
        starters.add(Networks.readNetwork("(((1:0.96506882,2:0.96506882):3.42176819,((3:0.68811798,4:0.68811798):1.97950363,5:2.66762161):1.71921539):6.10220146,(((6:0.07656860,7:0.07656860):1.15059662,8:1.22716522):8.26786232,(9:4.70946312,10:4.70946312):4.78556442):0.99401093);"));
        starters.add(Networks.readNetwork("(((((1:0.97398567,2:0.97398567):1.41451073,3:2.38849640):0.76133728,(4:0.30510712,5:0.30510712):2.84472656):3.37349415,(6:0.66458702,(7:0.13372993,8:0.13372993):0.53085709):5.85874081):1.54096222,(9:3.92808533,10:3.92808533):4.13620472);"));
        starters.add(Networks.readNetwork("(((1:1.38463020,2:1.38463020):7.23227310,(3:0.42516518,4:0.42516518):8.19173813):1.31324577,((5:4.79586792,(6:1.54055214,7:1.54055214):3.25531578):4.52950287,(8:4.13646507,(9:0.62176704,10:0.62176704):3.51469803):5.18890572):0.60477829);"));

        int numTaxa = 5;
        int maxReti = 1;
        int maxIndi = 4;
        int maxNumUnder = numTaxa;
        int maxDiameter = numTaxa + 1;
        boolean usePseudoLikelihood = true;

        double timeConsumption[][] = new double[maxNumUnder + 1][maxDiameter + 1];
        for(int i = 4 ; i <= maxIndi ; i++) {


            System.out.println("# of taxa = " + numTaxa);
            System.out.println("# of individuals/species = " + i);

            State state = null;


            for(int iter = 0 ; iter <= 150 ; iter++) {



                if(iter == 0) {


                    Network<NetNodeInfo> currentNetwork = Networks.readNetworkWithRootPop("[0.006](((((Q:0.004:0.006)I5#H1:0.002:0.005:0.7,A:0.006:0.006)I3:0.016:0.005,L:0.022:0.006)I2:0.02:0.005,(I5#H1:0.003:0.005:0.3,R:0.007:0.006)I4:0.035:0.005)I1:0.038:0.005,C:0.08:0.006);");
                    initNetHeights(currentNetwork, 0.01);

                    Map<String, String> alleles2species = new HashMap<>();
                    Map<String, List<String>> species2alleles = new HashMap<>();
                    for(NetNode node : currentNetwork.getLeaves()) {
                        String s = node.getName();

                        species2alleles.put(s, new ArrayList<>());
                        for(int k = 0 ; k < i ; k++) {
                            String s1 = String.format("%s_%d", s, k);
                            alleles2species.put(s1, s);
                            species2alleles.get(s).add(s1);
                        }
                    }

                    System.out.println("# of reticulations = " + currentNetwork.getReticulationCount());

                    SimSNPInNetwork simulator = new SimSNPInNetwork(BAGTRModel, 12345678L);
                    simulator._diploid = false;
                    Map<String, String> onesnp = simulator.generateSNPs(currentNetwork, species2alleles, numSites, !useOnlyPolymorphic);

                    List<MarkerSeq> alns = new ArrayList<>();
                    MarkerSeq aln = new MarkerSeq(onesnp);
                    aln._diploid = false;
                    alns.add(aln);
                    alns.get(0)._RPatterns = SNAPPLikelihood.haploidSequenceToPatterns(alleles2species, alns);

                    Utils._NET_MAX_RETI = 1;
                    Utils.DISABLE_PARAMETER_MOVES = true;

                    state = new State(
                            Networks.getFullString(currentNetwork),
                            Utils._START_GT_LIST,
                            alns,
                            Utils._POISSON_PARAM,
                            species2alleles,
                            BAGTRModel
                    );
                    state.calculatePrior();

                }
                state.propose();
                int numUnder = CoalescenceHistoriesCounting.getNumTaxaUnderReticulation(state.getNetworkObject());
                int diameter = CoalescenceHistoriesCounting.getDiscreteDiameter(state.getNetworkObject());


                if(!state.isValidState() || state.getNetworkObject().getReticulationCount() != 1 || numUnder < Math.min(iter / 50, 4)) {
                    state.undo(Utils.INVALID_MOVE);
                    iter--;
                    continue;
                }


                if(diameter > 5) {
                    System.out.println(Networks.getTopologyString(state.getNetworkObject()));
                }

                long start = System.currentTimeMillis();
                double ll;
                try {
                    ll = state.calculateLikelihood();
                }catch(Exception e) {
                    e.printStackTrace();
                    continue;
                }
                long end = System.currentTimeMillis();
                System.out.println("Time(s): " + (end - start) / 1000.0);


                if(iter < 10) {
                    continue;
                }

                timeConsumption[numUnder][diameter] = Math.max(timeConsumption[numUnder][diameter], (end - start) / 1000.0);

                System.out.println("Iter = " + iter);
                System.out.println("likelihood = " + ll + " numUnder = " + numUnder + " diameter = " + diameter);
            }
        }

        for(int i = 0 ; i <= maxNumUnder ; i++) {
            for(int j = 0 ; j <= maxDiameter ; j++) {
                System.out.print(timeConsumption[i][j] + " ");
            }
            System.out.println();
        }
    }

    public static void main(String[] args) {
        go();
    }
}

