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
            Tree gt = gts.get(0);
            Map<String, String> snp = _snpsim.generateSingeSite(gt);

            boolean all1 = true;
            boolean all0 = true;
            for (String name : snp.keySet()) {
                if(snp.get(name).equals("0")) all1 = false;
                if(snp.get(name).equals("1")) all0 = false;
            }

            if(!allowConstant && (all1 || all0)) {
                i--;
                continue;
            }

            for (String name : snp.keySet()) {
                if (!res.containsKey(name))
                    res.put(name, new StringBuilder());
                res.get(name).append(snp.get(name));
                //res.put(name, res.get(name) + snp.get(name));
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
        trueNetwork = Networks.readNetwork("(((((Q:0.5)I8#H1:0.5::0.7,A:1.0)I4:1.0,L:2.0)I3:2.0,(I8#H1:1.0::0.3,R:1.5)I7:2.5)I2:1.0,(G:1.0,C:1.0)I1:4.0)I0;");
        trueNetwork = Networks.readNetwork("(((((Q:0.5)I8#H1:0.5::0.7,(A:0.5)I6#H2:0.5::0.8)I4:1.0,L:2.0)I3:2.0,(I8#H1:1.0::0.3,R:1.5)I7:2.5)I2:1.0,((I6#H2:0.5::0.2,C:1.0)I5:1.0,G:2.0)I1:3.0)I0;");

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


        trueNetwork = Networks.readNetwork("((G:0.018380841366896713,C:0.018380841366896713):0.07259787745284114,(((Q:0.019223930643805762)#H1:0.0079382721121616::0.3027236780708331,R:0.02716220275596736):0.04612672867584792,((L:0.02760957365672652)#H2:0.01933595272319038::0.07225913233066261,((#H1:3.112852137000452E-4::0.6972763219291669,A:0.019535215857505808):0.017171505371586866,#H2:0.009097147572366154::0.9277408676693374):0.010238805150824225):0.02634340505189838):0.01768978738792258);");
        trueNetwork.getRoot().setRootPopSize(0.03434332011291962);
        trueNetwork = Networks.readNetwork("(((((R:0.03692297297059491)#H1:0.026335083580234356::0.32659382420676164,(A:0.048296112704057884,(Q:0.0016061601679680012)#H2:0.04668995253608988::0.8824921981923378):0.01496194384677138):4.711142085103842E-4,((L:0.03025867321849643)#H3:0.03261745957608034::0.675454288991292,((G:0.0514837139295224,C:0.0514837139295224):6.708696808103431E-4)#H4:0.010721549184244027::0.16076035943533917):8.530379647628816E-4):0.04434852129641825,#H4:0.05592310844542516::0.8392396405646608):9.399206670118038E-4,((#H3:0.003242697377805359::0.324545711008708,#H2:0.031895210428333785::0.11750780180766218):0.058420991394148414,#H1:0.054999389019855294::0.6734061757932384):0.017095250732319503);");
        trueNetwork.getRoot().setRootPopSize(6.14445682191189E-4);
        //trueNetwork = Networks.readNetwork("((((C:0.07,(B:0.06)I7#H1:0.01::0.2)I4:0.02,(I7#H1:0.02::0.8)I6#H2:0.01::0.2)I3:0.02,(I6#H2:0.02::0.8)I5#H3:0.01::0.2)I2:0.07,(I5#H3:0.06::0.8,A:0.16)I1:0.02)I0;");
        //trueNetwork.getRoot().setRootPopSize(0.04);

        int nameCount = 0;
        for(Object nodeObject : Networks.postTraversal(trueNetwork)) {
            NetNode node = (NetNode) nodeObject;
            for(Object parentObject :  node.getParents()) {
                NetNode parent = (NetNode) parentObject;
                node.setParentSupport(parent, trueNetwork.getRoot().getRootPopSize());
            }
        }
        for(Object node : trueNetwork.dfs()) {
            NetNode mynode = (NetNode) node;
            if(mynode.getName().equals("")) {
                mynode.setName("II" + nameCount);
                nameCount++;
            }
        }


        SimSNPInNetwork simulator = new SimSNPInNetwork(BAGTRModel, 12345678L + 1000000);
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
                    SNAPPAlgorithm run = new SNAPPAlgorithm(cloneNetwork, BAGTRModel, null);
                    double likelihood = 0;
                    try {
                        likelihood = run.getProbability(converter);
                        sum += likelihood;
                        System.out.println(represent + " " + likelihood  + " " + count );
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
