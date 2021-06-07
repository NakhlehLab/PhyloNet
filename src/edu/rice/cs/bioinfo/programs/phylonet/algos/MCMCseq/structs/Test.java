package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.core.MC3Core;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.move.network.param.ChangeInheritance;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimGTInNetworkByMS;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimGTInNetworkWithTheta;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimSeqInGTBySeqGen;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;
import org.jblas.util.Random;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 9/1/18
 * Time: 5:28 PM
 * To change this template use File | Settings | File Templates.
 */
public class Test {
    static void testSamplingEmbedding() {
        double alpha1 = 0.00002;//0.5;
        double alpha2 = 0.00001;//0.4;
        double alpha3 = 0.00002;//0.6;
        double alpha4 = 0.00003;//0.8;
        double a = 0;
        double b = 0;
        double c = 0;
        double theta = 0.1;
        double gamma = 0.905;

        String netstring = String.format("[%f](((B:%f)I3#H1:%f::%f,A:%f)I1:%f,(C:%f,I3#H1:%f::%f)I2:%f)I0;", theta, alpha2, alpha1-alpha2, gamma, alpha1, alpha4-alpha1, alpha3, alpha3-alpha2, 1.0-gamma, alpha4-alpha3);
        //String netstring = "[0.02](((((C:0.01)#H1:0.02::0.6,B:0.03):0.01)#H2:0.02::0.7,A:0.06):0.03,((D:0.02,#H1:0.01::0.4):0.03,#H2:0.01::0.3):0.04);";
        //String netstring = "[0.03478971320729695](((((D:0.018773841591172784,C:0.018773841591172784):0.01665803647774957,B:0.035431878068922354):5.193091455444887E-4)#H1:0.042646185645307325::0.490471497495645,((#H1:0.01768017414602177::0.509528502504355,A:0.053631361360488614):2.981418077369069E-5)#H2:0.024936197318511863::0.28337456507014525):2.9602647955248584E-4,#H2:0.02523222379806435::0.7166254349298548);";
        //String netstring = "[0.07138571458425186]((((B:0.021545725769406308,(D:0.02072866724785186,C:0.02072866724785186):8.170585215544479E-4):0.012235537519323995)#H1:0.02015659627482603::0.17371727563554662,A:0.05393785956355633):0.005245366442141658,#H1:0.025401962716967687::0.8262827243644534);";
        Network<Object> trueNetwork;
        trueNetwork = Networks.readNetworkWithRootPop(netstring);
        System.out.println(Networks.getFullString(trueNetwork));

        String gtstring = "((A:1.0,B:1.0):0.3,C:1.3);";
        //String gtstring = "(((A:0.2,B:0.2):0.3,C:0.5):0.2,D:0.7);";
        //String gtstring = "((B:0.07,A:0.07)I1:0.03,(D:0.08,C:0.08)I2:0.02)I0;";
        String gtstring1 = "((B:0.12,A:0.12)I1:0.03,(D:0.13,C:0.13)I2:0.02)I0;";
        Tree gt0 = Trees.readTree(gtstring);
        List<UltrametricTree> gts = new ArrayList<>();
        for(int i = 0 ; i < 1 ; i++)
            gts.add(new UltrametricTree(Trees.readTree(gtstring)));
        //gts.add(new UltrametricTree(gt0));
        //gts.add(new UltrametricTree(gt0));

        //Utils._NUM_THREADS = 1;
        Utils._SEED = Random.nextInt(99999999);
        Utils._ESTIMATE_POP_SIZE = false;
        Utils._CONST_POP_SIZE = true;
        Utils.SAMPLE_EMBEDDINGS = true;
        UltrametricNetwork network = new UltrametricNetwork(trueNetwork, gts);
        network.setDirty(true);
        double prevLikelihood = network.logDensity();
        network._logGeneTreeNetwork = network.logDensity1();
        int acc = 0;
        TreeEmbedding face = null;
        ChangeInheritance changeInheritance = new ChangeInheritance(network);
        int probBins[] = new int[10];

        if(Utils.SAMPLE_EMBEDDINGS) {
            Map<Integer, Double> likelihood2 = new TreeMap<>();
            Map<Integer, Integer> embeddingCount2 = new TreeMap<>();

            Map<TreeEmbedding, Double> likelihood = new HashMap<>();
            Map<TreeEmbedding, Integer> embeddingCount = new HashMap<>();
            Map<TreeEmbedding, TreeEmbedding> embeddings = new HashMap<>();
            for (int i = 0; i < 10000; i++) {
                changeInheritance.propose();
                network.setDirty(true);
                //network._geneTrees.get(Randomizer.getRandomInt(gts.size()) ).setDirty(true);
                double embeddingLogHr = network.rebuildEmbeddings();
                //System.out.println(embeddingLogHr);
                double nextLikelihood = network.logDensity();
                //System.out.println(embeddingLogHr);
                double logAlpha = nextLikelihood - prevLikelihood + embeddingLogHr;
                System.out.println(embeddingLogHr + " " + nextLikelihood + " " + logAlpha);

                if (logAlpha >= Math.log(Randomizer.getRandomDouble())) {
                    network.accept();
                    for (UltrametricTree gt : network._geneTrees) gt.accept();
                    acc++;
                    prevLikelihood = nextLikelihood;
                } else {
                    network.undo();
                    network.reject();
                    //for (UltrametricTree gt : network._geneTrees) gt.undo();
                    for (UltrametricTree gt : network._geneTrees) gt.reject();
                }

                TreeEmbedding embedding = network._embeddings.get(0).clone();
                if (!embeddingCount.containsKey(embedding)) {
                    embeddingCount.put(embedding, 0);
                }
                embeddingCount.put(embedding, embeddingCount.get(embedding) + 1);
                likelihood.put(embedding, network._logGeneTreeNetwork[0]);

                if(gts.size() >= 2) {
                    if(face == null) face = network._embeddings.get(0).clone();

                    int s = 0;
                    for(int j = 0 ; j < gts.size() ; j++) {
                        if(network._embeddings.get(j).prob > 0.5) s++;
                    }

                    likelihood2.put(s, network._logGeneTreeNetwork[0]);
                    if (!embeddingCount2.containsKey(s)) {
                        embeddingCount2.put(s, 0);
                    }
                    embeddingCount2.put(s, embeddingCount2.get(s) + 1);
                }

                if(embeddings.containsKey(embedding) && embedding.probsum != embeddings.get(embedding).probsum) {
                    throw new RuntimeException("!!!");
                }
                if(embeddings.containsKey(embedding) && embedding.prob != embeddings.get(embedding).prob) {
                    //throw new RuntimeException("!!!");
                }
                embeddings.put(embedding, embedding);

                // System.out.println("!!!!! " + network.logDensity());
            }

            double sum = 0;
            for (TreeEmbedding embedding : embeddingCount.keySet()) {
                System.out.println(embeddingCount.get(embedding) + " " + likelihood.get(embedding) + " " + embeddings.get(embedding).prob + " " + embeddings.get(embedding).probsum);
                sum += Math.exp(likelihood.get(embedding));
            }

            System.out.println("Embedding2");
            for (int s : embeddingCount2.keySet()) {
                System.out.println(s + " " + embeddingCount2.get(s) + " " + likelihood2.get(s));
            }
            System.out.println("Embedding2");

            System.out.println(Math.log(sum));
            System.out.println("prev likelihood " + prevLikelihood);
            System.out.println("Acc " + acc);
        } else {
            System.out.println(prevLikelihood);
        }
    }

    static void testWhetherBiasedOrNot(String[] args) {
        double alpha1 = 0.5;
        double alpha2 = 0.4;
        double alpha3 = 0.6;
        double alpha4 = 0.8;
        double a = 0;
        double b = 0;
        double c = 0;
        double theta = 0.1;
        double gamma = 0.55;

        String netstring = String.format("[%f](((B:%f)I3#H1:%f::%f,A:%f)I1:%f,(C:%f,I3#H1:%f::%f)I2:%f)I0;", theta, alpha2, alpha1-alpha2, gamma, alpha1, alpha4-alpha1, alpha3, alpha3-alpha2, 1.0-gamma, alpha4-alpha3);
        //String netstring = "[0.02](((((C:0.01)#H1:0.02::0.6,B:0.03):0.01)#H2:0.02::0.7,A:0.06):0.03,((D:0.02,#H1:0.01::0.4):0.03,#H2:0.01::0.3):0.04);";
        //String netstring = "[0.03478971320729695](((((D:0.018773841591172784,C:0.018773841591172784):0.01665803647774957,B:0.035431878068922354):5.193091455444887E-4)#H1:0.042646185645307325::0.490471497495645,((#H1:0.01768017414602177::0.509528502504355,A:0.053631361360488614):2.981418077369069E-5)#H2:0.024936197318511863::0.28337456507014525):2.9602647955248584E-4,#H2:0.02523222379806435::0.7166254349298548);";
        //String netstring = "[0.07138571458425186]((((B:0.021545725769406308,(D:0.02072866724785186,C:0.02072866724785186):8.170585215544479E-4):0.012235537519323995)#H1:0.02015659627482603::0.17371727563554662,A:0.05393785956355633):0.005245366442141658,#H1:0.025401962716967687::0.8262827243644534);";
        Network<NetNodeInfo> trueNetwork = Networks.readNetwork("I0;");
        trueNetwork = Networks.readNetworkWithRootPop(netstring);
        System.out.println(Networks.getFullString(trueNetwork));
        Map<String, List<String>> species2alleles = new HashMap<>();
        species2alleles.put("I0", new ArrayList<>());
        species2alleles.get("I0").add("A");
        species2alleles.get("I0").add("B");
        species2alleles.get("I0").add("C");


        List<UltrametricTree> gts = new ArrayList<>();

        SimGTInNetworkWithTheta simulator = new SimGTInNetworkWithTheta();
        simulator.setSeed((long) Randomizer.getRandomInt(99999999));
        int numGT = 100;
        Network<NetNodeInfo> net = Networks.readNetworkWithRootPop(String.format("[%f](((B:%f)I3#H1:%f::%f,A:%f)I1:%f,(C:%f,I3#H1:%f::%f)I2:%f)I0;", theta, alpha2, alpha1-alpha2, 0.1, alpha1, alpha4-alpha1, alpha3, alpha3-alpha2, 0.9, alpha4-alpha3));//Networks.readNetworkWithRootPop("[0.1]I0;");

        for(NetNode<NetNodeInfo> node : net.dfs()){
            for(NetNode<NetNodeInfo> parent : node.getParents()) {
                if(node.getParentSupport(parent) == NetNode.NO_SUPPORT) {
                    node.setParentSupport(parent, net.getRoot().getRootPopSize());
                }
            }
        }

        List<Tree> simulatedGTs = new ArrayList<>();
        for(int i = 0 ; i < numGT ; i++) {
            List<Tree> onegt = simulator.generateGTs(net, null, 1.0, 1);
            Trees.autoLabelNodes((STITree)onegt.get(0));
            simulatedGTs.add(onegt.get(0));
        }


        for(int i = 0 ; i < numGT ; i++)
            gts.add(new UltrametricTree(simulatedGTs.get(i)));
        //gts.add(new UltrametricTree(gt0));
        //gts.add(new UltrametricTree(gt0));

        //Utils._NUM_THREADS = 1;
        Utils._SEED = Random.nextInt(99999999);
        Utils._ESTIMATE_POP_SIZE = false;
        Utils._CONST_POP_SIZE = true;
        Utils.SAMPLE_EMBEDDINGS = true;
        UltrametricNetwork network = new UltrametricNetwork(trueNetwork, gts);
        network.setDirty(true);
        double prevLikelihood = network.logDensity();
        network._logGeneTreeNetwork = network.logDensity1();
        int acc = 0;
        TreeEmbedding face = null;
        ChangeInheritance changeInheritance = new ChangeInheritance(network);
        int probBins[] = new int[10];
        double sumProb = 0.0;

        if(true) {
            Map<Integer, Double> likelihood2 = new TreeMap<>();
            Map<Integer, Integer> embeddingCount2 = new TreeMap<>();

            Map<TreeEmbedding, Double> likelihood = new HashMap<>();
            Map<TreeEmbedding, Integer> embeddingCount = new HashMap<>();
            Map<TreeEmbedding, TreeEmbedding> embeddings = new HashMap<>();
            for (int i = 0; i < 100000; i++) {
                changeInheritance.propose();
                network.setDirty(true);
                //network._geneTrees.get(Randomizer.getRandomInt(gts.size()) ).setDirty(true);
                double embeddingLogHr = network.rebuildEmbeddings();
                //System.out.println(embeddingLogHr);
                double nextLikelihood = network.logDensity();
                //System.out.println(embeddingLogHr);
                double logAlpha = nextLikelihood - prevLikelihood + embeddingLogHr;
                System.out.println(embeddingLogHr + " " + nextLikelihood + " " + logAlpha);

                if (logAlpha >= Math.log(Randomizer.getRandomDouble())) {
                    network.accept();
                    for (UltrametricTree gt : network._geneTrees) gt.accept();
                    acc++;
                    prevLikelihood = nextLikelihood;
                } else {
                    changeInheritance.undo();
                    network.undo();
                    network.reject();
                    //for (UltrametricTree gt : network._geneTrees) gt.undo();
                    for (UltrametricTree gt : network._geneTrees) gt.reject();
                }

                NetNode<NetNodeInfo> AParent = network.getNetwork().findNode("A").getParents().iterator().next();
                NetNode<NetNodeInfo> HNode = Networks.getOtherChild(AParent, network.getNetwork().findNode("A"));
                double iprob = HNode.getParentProbability(AParent);
                probBins[(int)(iprob / (1.0 / probBins.length))]++;

//                if(gts.size() >= 2 && Utils.SAMPLE_EMBEDDINGS) {
//                    if(face == null) face = network._embeddings.get(0).clone();
//
//                    int s = 0;
//                    for(int j = 0 ; j < gts.size() ; j++) {
//                        if(network._embeddings.get(j).prob > 0.5) s++;
//                    }
//
//                    likelihood2.put(s, network._logGeneTreeNetwork[0]);
//                    if (!embeddingCount2.containsKey(s)) {
//                        embeddingCount2.put(s, 0);
//                    }
//                    embeddingCount2.put(s, embeddingCount2.get(s) + 1);
//                }



                // System.out.println("!!!!! " + network.logDensity());
            }

            for(int i = 0 ; i < probBins.length ; i++) {
                System.out.println(i * (1.0 / probBins.length) + "->" + (i + 1) * (1.0 / probBins.length) + " " + probBins[i]);
            }

            double sum = 0;


            System.out.println("Embedding2");
            for (int s : embeddingCount2.keySet()) {
                System.out.println(s + " " + embeddingCount2.get(s) + " " + likelihood2.get(s));
            }
            System.out.println("Embedding2");

            System.out.println(Math.log(sum));
            System.out.println("prev likelihood " + prevLikelihood);
            System.out.println("Acc " + acc);
        } else {
            System.out.println(prevLikelihood);
        }
    }

    static void testSampler(String[] args) {
        double alpha1 = 3;
        double alpha2 = 2;
        double alpha3 = 4;
        double alpha4 = 8;
        double a = 0;
        double b = 0;
        double c = 0;
        double theta = 0.02;
        double gamma = 0.55;

        String netstring = String.format("[%f](((B:%f)I3#H1:%f::%f,A:%f)I1:%f,(C:%f,I3#H1:%f::%f)I2:%f)I0;", theta, alpha2, alpha1-alpha2, gamma, alpha1, alpha4-alpha1, alpha3, alpha3-alpha2, 1.0-gamma, alpha4-alpha3);
        //String netstring = "[0.02](((((C:0.01)#H1:0.02::0.6,B:0.03):0.01)#H2:0.02::0.7,A:0.06):0.03,((D:0.02,#H1:0.01::0.4):0.03,#H2:0.01::0.3):0.04);";
        //String netstring = "[0.03478971320729695](((((D:0.018773841591172784,C:0.018773841591172784):0.01665803647774957,B:0.035431878068922354):5.193091455444887E-4)#H1:0.042646185645307325::0.490471497495645,((#H1:0.01768017414602177::0.509528502504355,A:0.053631361360488614):2.981418077369069E-5)#H2:0.024936197318511863::0.28337456507014525):2.9602647955248584E-4,#H2:0.02523222379806435::0.7166254349298548);";
        //String netstring = "[0.07138571458425186]((((B:0.021545725769406308,(D:0.02072866724785186,C:0.02072866724785186):8.170585215544479E-4):0.012235537519323995)#H1:0.02015659627482603::0.17371727563554662,A:0.05393785956355633):0.005245366442141658,#H1:0.025401962716967687::0.8262827243644534);";
        Network<NetNodeInfo> trueNetwork = Networks.readNetwork("I0;");
        trueNetwork = Networks.readNetworkWithRootPop(netstring);
        System.out.println(Networks.getFullString(trueNetwork));
        Map<String, List<String>> species2alleles = new HashMap<>();
        species2alleles.put("A", new ArrayList<>());
        species2alleles.get("A").add("A_0");
        species2alleles.get("A").add("A_1");
        species2alleles.put("B", new ArrayList<>());
        species2alleles.get("B").add("B_0");
        species2alleles.get("B").add("B_1");
        species2alleles.put("C", new ArrayList<>());
        species2alleles.get("C").add("C_0");
        species2alleles.get("C").add("C_1");

        List<UltrametricTree> gts = new ArrayList<>();

        int numGT = 100;
        String MSPath = "/Users/zhujiafan/Documents/Luay/msdir/ms";
        Network<NetNodeInfo> net = trueNetwork.clone();//Networks.readNetworkWithRootPop("[0.1]I0;");

        SimGTInNetworkByMS simGT = new SimGTInNetworkByMS();
        List<Tree> simulatedGTs = simGT.generateGTs(net, species2alleles, numGT, MSPath);

        List<Map<String, String>> lociSeq = new ArrayList<>();
        double freq[] = new double[]{0.2112,0.2888,0.2896,0.2104};
        double rates[] = new double[]{0.2173,0.9798,0.2575,0.1038,1,0.2070};
        int lenPerGT = 100;
        String SeqGenPath = args.length > 0 ? args[0] : "/Users/zhujiafan/Documents/Luay/Seq-Gen.v1.3.3/source/seq-gen";
        for(Tree gt : simulatedGTs) {
            lociSeq.add(SimSeqInGTBySeqGen.execute(gt, theta, freq, rates, lenPerGT, null, SeqGenPath));
        }

        List<Alignment> alignments = new ArrayList<>();

        for(int locus_num = 0 ; locus_num < numGT ; locus_num++){
            alignments.add(new Alignment(lociSeq.get(locus_num), "locus_" + locus_num));
        }


        Map<String, String> geneTrees = new HashMap<>();
        for(int locus_num = 0 ; locus_num < numGT ; locus_num++){
            try {
                STITree<Object> tree = new STITree(simulatedGTs.get(locus_num).toNewick());
                for (Object t : tree.postTraverse()) {
                    STINode tt = (STINode) t;
                    if (!tt.isLeaf()) {
                        tt.setName("");
                    }

                    tt.setParentDistance(tt.getParentDistance() * theta / 2.0);
                }
                geneTrees.put("locus_" + locus_num, tree.toNewick());

            } catch(Exception e) {
                e.printStackTrace();
            }

        }


        List<Map<String, String>> startGTs = new ArrayList<>();

        startGTs.add(geneTrees);
        Utils._START_GT_LIST = startGTs;
        Utils._FIX_GENE_TREES = false;
        Utils._FIX_GENE_TREE_TOPOLOGIES = false;
        List<String> snets = new ArrayList<>();

        for (Object nodeObject : Networks.postTraversal(trueNetwork)) {
            NetNode node = (NetNode) nodeObject;
            for (Object parentObject : node.getParents()) {
                NetNode parent = (NetNode) parentObject;
                node.setParentSupport(parent, theta);
                node.setParentDistance(parent, node.getParentDistance(parent) * theta / 2.0);
            }
        }
        trueNetwork.getRoot().setRootPopSize(theta);

        snets.add(Networks.getFullString(trueNetwork));
        //Utils._START_NET = snets;

        Utils._TAXON_MAP = species2alleles;

        Utils._NET_MAX_RETI = 1;
        Utils._FIX_NET_TOPOLOGY = false;
        Utils._FIX_NET = false;

        Utils._ESTIMATE_POP_SIZE = false;
        Utils._ESTIMATE_POP_PARAM = false;
        Utils._POP_SIZE_MEAN = 0.02;

        Utils.SAMPLE_EMBEDDINGS = true;

        Utils._NUM_THREADS = 8;

        System.out.println("MCMC_SEQ");

        Collections.sort(alignments);
        MC3Core mc3 = new MC3Core(alignments);
        mc3.run();


    }

    static void testInheritanceProbMove() {
        double oldRate = 0.5;
        double _windowSize = 0.1;
        int probBins[] = new int[10];
        for(int i = 0 ; i < 100000 ; i++) {
            double newRate = oldRate + (Randomizer.getRandomDouble() - 0.50) * _windowSize;

            while (newRate > 1 || newRate < 0) {
                if (newRate > 1.0) newRate = 2.0 - newRate;
                if (newRate < 0.0) newRate = -newRate;
            }
            oldRate = newRate;
            probBins[(int)(newRate / (1.0 / probBins.length))]++;
        }

        for(int i = 0 ; i < probBins.length ; i++) {
            System.out.println(i * (1.0 / probBins.length) + "->" + (i + 1) * (1.0 / probBins.length) + " " + probBins[i]);
        }
    }

    public static void main(String[] args) {
        testWhetherBiasedOrNot(args);

    }
}
