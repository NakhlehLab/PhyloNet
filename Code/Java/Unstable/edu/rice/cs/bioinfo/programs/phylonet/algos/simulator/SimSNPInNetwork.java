package edu.rice.cs.bioinfo.programs.phylonet.algos.simulator;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.SNAPPAlgorithm;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.BiAllelicGTR;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations.OneNucleotideObservation;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.*;
import java.util.concurrent.TimeUnit;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 11/27/16
 * Time: 11:50 AM
 * To change this template use File | Settings | File Templates.
 */
public class SimSNPInNetwork {
    private Long _seed = null;
    private BiAllelicGTR _model;
    private SimGTInNetworkWithTheta _gtsim;
    private SimSNPInGT _snpsim;

    public SimSNPInNetwork(BiAllelicGTR model, Long seed) {
        _model = model;
        _seed = seed;
        _gtsim = new SimGTInNetworkWithTheta();
        _gtsim.setSeed(_seed);
        _snpsim = new SimSNPInGT(model);
        _snpsim.setSeed(_seed);
    }

    public Map<String, String> generateSNPs(Network network, Map<String, List<String>> species2alleles, int numGTs, boolean allowConstant){

        BiAllelicGTR model = _model;
        double pi0 = model.getEquilibriumVector().get(0, 0);
        double pi1 = model.getEquilibriumVector().get(1, 0);
        double u = 1.0 / (2.0 * pi0);
        double v = 1.0 / (2.0 * pi1);
        double mu = 2.0 * u * v / (u + v);
        Map<String, StringBuilder> res = new HashMap<>();
        for(int i = 0 ; i < numGTs ; i++) {
            List<Tree> gts = _gtsim.generateGTs(network, species2alleles, mu, 1);

            for (Tree gt : gts) {
                Map<String, String> snp = _snpsim.generateSingeSite(gt, allowConstant);
                for (String name : snp.keySet()) {
                    if (!res.containsKey(name))
                        res.put(name, new StringBuilder());
                    res.get(name).append(snp.get(name));
                    //res.put(name, res.get(name) + snp.get(name));
                }
            }
        }

        Map<String, String> ret = new HashMap<>();
        for(String s : res.keySet()) {
            ret.put(s, res.get(s).toString());
        }
        return ret;

    }

    public static void main(String []args) {
        double pi0 = 0.5;
        double pi1 = 1- pi0;

        BiAllelicGTR BAGTRModel = new BiAllelicGTR(new double[] {pi0, pi1}, new double[] {1.0/ (2.0 * pi0)});

        double gamma = 1 - 0.3;

        double y = 1.5;
        double x = 0.5;

        Network<Object> trueNetwork;
        trueNetwork = Networks.readNetwork("((((b:1.0,c:1.0)I4:" + y + ")I3#H1:0.1::" + (1 - gamma) + ",a:" + (1.1 + y) + ")I1:" + x + ",(I3#H1:0.1::" + gamma + ",d:" + (1.1 + y) + ")I2:"+ x +")I0;");
        trueNetwork = Networks.readNetwork("(((((A:0.7)I6#H1:1.3::0.8,Q:2.0)I4:1.0,L:3.0)I3:1.0,R:4.0)I2:1.0,(G:2.0,(I6#H1:0.7::0.2,C:1.4)I5:0.6)I1:3.0)I0;");
        trueNetwork = Networks.readNetwork("(((((Q:2.0,A:2.0)I4:1.0,L:3.0)I3:0.5)I8#H1:0.5::0.7,R:4.0)I2:1.0,(I8#H1:0.5::0.3,(G:2.0,C:2.0)I1:2.0)I7:1.0)I0;");
        //trueNetwork = Networks.readNetwork("(((((Q:0.5)I8#H1:0.5::0.7,A:1.0)I4:1.0,L:2.0)I3:2.0,(I8#H1:1.0::0.3,R:1.5)I7:2.5)I2:1.0,(G:1.0,C:1.0)I1:4.0)I0;");
        //trueNetwork = Networks.readNetwork("(((((Q:0.5)I8#H1:0.5::0.7,(A:0.5)I6#H2:0.5::0.8)I4:1.0,L:2.0)I3:2.0,(I8#H1:1.0::0.3,R:1.5)I7:2.5)I2:1.0,((I6#H2:0.5::0.2,C:1.0)I5:1.0,G:2.0)I1:3.0)I0;");

        System.out.println(trueNetwork.toString());

        double constTheta = 0.036;
        double mu = 1.0;
        for(Object nodeObject : Networks.postTraversal(trueNetwork)) {
            NetNode node = (NetNode) nodeObject;
            for(Object parentObject :  node.getParents()) {
                NetNode parent = (NetNode) parentObject;
                node.setParentSupport(parent, constTheta);
                node.setParentDistance(parent, node.getParentDistance(parent) * constTheta / 2.0);
            }
        }
        trueNetwork.getRoot().setRootPopSize(constTheta);


        /*trueNetwork = Networks.readNetwork("(((d:0.03312567269052971:0.05729812764763328)#H1:0.04040511494621882:0.01741043976665741:0.20775124134245948,(((b:0.028486216774249146:0.006372651326433015,c:0.028486216774249146:0.00978225510772307):0.022397718719795415:0.043337051035579346)#H2:0.005681182773572657:0.019241824721338403:0.7217411528293988,#H1:0.023439445577087506:0.07791563334890902:0.7922487586575405):0.01696566936913131:0.05670452561728086):0.5676859473369829:0.05332739995725745,(a:0.05264074053638943:0.003071386462697118,#H2:0.0017568050423448708:0.02272184117641578:0.27825884717060123):0.588575994437342:0.008691349345198683);");
        trueNetwork.getRoot().setRootPopSize(0.08023790737422791);
        int nameCount = 0;
        for(Object node : trueNetwork.dfs()) {
            NetNode mynode = (NetNode) node;
            if(mynode.getName().equals("")) {
                mynode.setName("I" + nameCount);
                nameCount++;
            }
        }
*/

        SimSNPInNetwork simulator = new SimSNPInNetwork(BAGTRModel, null);
        Map<String, String> onesnp = simulator.generateSNPs(trueNetwork, null, 1000000, true);

        List<Map<String, String>> snpdata = new ArrayList<>();
        snpdata.add(onesnp);
        List<Alignment> alns = new ArrayList<>();
        for(Map<String, String> input : snpdata) {
            for(String allele : input.keySet()) {
                System.out.println(allele + " " + input.get(allele));
            }
            System.out.println();
            Alignment aln = new Alignment(input);

            Map<Integer, Integer> cache = new HashMap<>();

            for(int i = 0 ; i < aln.getSiteCount() ; i++) {
                Map<String, Character> colorMap = new TreeMap<String, Character>();
                for(String taxon : aln.getAlignment().keySet()) {
                    colorMap.put(taxon, aln.getAlignment().get(taxon).charAt(i));
                }
                Integer represent = 0;
                for(String s : colorMap.keySet()) {
                    represent = (represent << 1) + (colorMap.get(s) == '0' ? 0 : 1);
                }
                if(!cache.containsKey(represent)) {
                    cache.put(represent, 0);
                }
                cache.put(represent, cache.get(represent) + 1);
            }

            aln.setCache(cache);
            alns.add(aln);


            for(Alignment alg : alns) {
                List<String> names = alg.getTaxaNames();
                double sum = 0;

                for(Integer represent : alg.getCache().keySet()) {
                    int cur = represent;
                    Integer count = alg.getCache().get(cur);
                    Map<String, Character> colorMap = new HashMap<String, Character>();
                    for(int i = names.size() - 1 ; i >= 0 ; i--) {
                        Character site = cur % 2 == 1 ? '1' : '0';
                        colorMap.put(names.get(i), site);
                        cur /= 2;
                    }
                    OneNucleotideObservation converter = new OneNucleotideObservation(colorMap);
                    Network cloneNetwork = Networks.readNetwork(trueNetwork.toString());
                    cloneNetwork.getRoot().setRootPopSize(trueNetwork.getRoot().getRootPopSize());
                    SNAPPAlgorithm run = new SNAPPAlgorithm(cloneNetwork, BAGTRModel, 1);
                    double likelihood = 0;
                    try {
                        likelihood = run.getProbability(converter);
                        sum += likelihood;
                        System.out.println(represent + " " + likelihood + " " + count);
                    } catch(Exception e) {
                        e.printStackTrace();
                        System.out.println("Exceptional network");
                    }


                }
                System.out.println("sum " + sum);


            }
        }
    }
}
