package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.distribution;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.core.MC3Core;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.alignment.MarkerSeq;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.move.network.param.ChangeInheritance;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.move.network.param.ChangeTime;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.move.network.param.ScaleRootTime;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.UltrametricNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimSNPInNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.BiAllelicGTR;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.structs.NetNodeInfo;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 3/2/19
 * Time: 4:45 PM
 * To change this template use File | Settings | File Templates.
 */
public class Test {

    static void testSampler(String[] args) {
        double pi0 = 0.5;
        double pi1 = 1- pi0;

        BiAllelicGTR BAGTRModel = new BiAllelicGTR(new double[] {pi0, pi1}, new double[] {1.0/ (2.0 * pi0)});


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
        for(Object nodeObj : trueNetwork.dfs()) {
            NetNode node = (NetNode) nodeObj;
            for(Object parentObj : node.getParents()) {
                NetNode parent = (NetNode) parentObj;
                node.setParentSupport(parent, trueNetwork.getRoot().getRootPopSize());
            }
        }

        System.out.println(Networks.getFullString(trueNetwork));
        Map<String, List<String>> species2alleles = new HashMap<>();
        species2alleles.put("I0", new ArrayList<>());
        species2alleles.get("I0").add("A");
        species2alleles.get("I0").add("B");
        species2alleles.get("I0").add("C");


        int numSites = 10000;
        boolean useOnlyPolymorphic = false;
        SimSNPInNetwork simulator = new SimSNPInNetwork(BAGTRModel, 12345678L);
        Map<String, String> onesnp = simulator.generateSNPs(trueNetwork, null, numSites, !useOnlyPolymorphic);

        List<MarkerSeq> alns = new ArrayList<>();
        MarkerSeq aln = new MarkerSeq(onesnp);
        alns.add(aln);
        aln._RPatterns = SNAPPLikelihood.haploidSequenceToPatterns(null, alns);


        Utils._START_NET = Networks.getFullString(trueNetwork);
        Utils._POP_SIZE_MEAN = 0.01;
        Utils._ESTIMATE_POP_SIZE = false;
        Utils._ESTIMATE_POP_PARAM = false;

        Utils.DISABLE_TOPOLOGY_MOVES = true;

        Utils.SAMPLE_SPLITTING = true;

        Utils._NUM_THREADS = 2;

        System.out.println("MCMC_BiMarkers");

        MC3Core mc3 = new MC3Core(alns, BAGTRModel);
        mc3.run();

    }

    static void testBias(String[] args) {
        double pi0 = 0.5;
        double pi1 = 1- pi0;

        BiAllelicGTR BAGTRModel = new BiAllelicGTR(new double[] {pi0, pi1}, new double[] {1.0/ (2.0 * pi0)});


        double alpha1 = 0.5;
        double alpha2 = 0.4;
        double alpha3 = 0.6;
        double alpha4 = 0.8;
        double a = 0;
        double b = 0;
        double c = 0;
        double theta = 0.1;
        double gamma = 0.55;

        //String netstring = String.format("[%f](((B:%f)I3#H1:%f::%f,A:%f)I1:%f,(C:%f,I3#H1:%f::%f)I2:%f)I0;", theta, alpha2, alpha1-alpha2, gamma, alpha1, alpha4-alpha1, alpha3, alpha3-alpha2, 1.0-gamma, alpha4-alpha3);
        String netstring = String.format("[%f]((B:%f,A:%f)I1:%f,C:%f)I0;", theta, alpha2, alpha2, alpha4-alpha2, alpha4);

        //String netstring = "[0.02](((((C:0.01)#H1:0.02::0.6,B:0.03):0.01)#H2:0.02::0.7,A:0.06):0.03,((D:0.02,#H1:0.01::0.4):0.03,#H2:0.01::0.3):0.04);";
        //String netstring = "[0.03478971320729695](((((D:0.018773841591172784,C:0.018773841591172784):0.01665803647774957,B:0.035431878068922354):5.193091455444887E-4)#H1:0.042646185645307325::0.490471497495645,((#H1:0.01768017414602177::0.509528502504355,A:0.053631361360488614):2.981418077369069E-5)#H2:0.024936197318511863::0.28337456507014525):2.9602647955248584E-4,#H2:0.02523222379806435::0.7166254349298548);";
        //String netstring = "[0.07138571458425186]((((B:0.021545725769406308,(D:0.02072866724785186,C:0.02072866724785186):8.170585215544479E-4):0.012235537519323995)#H1:0.02015659627482603::0.17371727563554662,A:0.05393785956355633):0.005245366442141658,#H1:0.025401962716967687::0.8262827243644534);";
        Network<NetNodeInfo> trueNetwork = Networks.readNetwork("I0;");
        trueNetwork = Networks.readNetworkWithRootPop(netstring);
        for(Object nodeObj : trueNetwork.dfs()) {
            NetNode node = (NetNode) nodeObj;
            for(Object parentObj : node.getParents()) {
                NetNode parent = (NetNode) parentObj;
                node.setParentSupport(parent, trueNetwork.getRoot().getRootPopSize());
            }
        }

        System.out.println(Networks.getFullString(trueNetwork));
        Map<String, List<String>> species2alleles = new HashMap<>();
        species2alleles.put("I0", new ArrayList<>());
        species2alleles.get("I0").add("A");
        species2alleles.get("I0").add("B");
        species2alleles.get("I0").add("C");

        int numSites = 100;
        boolean useOnlyPolymorphic = false;
        SimSNPInNetwork simulator = new SimSNPInNetwork(BAGTRModel, 12345678L);
        Map<String, String> onesnp = simulator.generateSNPs(trueNetwork, null, numSites, !useOnlyPolymorphic);

        List<MarkerSeq> alns = new ArrayList<>();
        MarkerSeq aln = new MarkerSeq(onesnp);
        alns.add(aln);
        aln._RPatterns = SNAPPLikelihood.haploidSequenceToPatterns(null, alns);


        Utils._START_NET = Networks.getFullString(trueNetwork);
        Utils._POP_SIZE_MEAN = 0.01;
        Utils._ESTIMATE_POP_SIZE = false;
        Utils._ESTIMATE_POP_PARAM = false;

        Utils.DISABLE_TOPOLOGY_MOVES = true;

        Utils.SAMPLE_SPLITTING = true;

        Utils._NUM_THREADS = 2;

        UltrametricNetwork network = new UltrametricNetwork( Networks.getFullString(trueNetwork), alns, null, BAGTRModel);
        network.setDirty(true);
        double prevLikelihood = network.logDensity();
        ChangeInheritance changeInheritance = new ChangeInheritance(network);
        ChangeTime changeTime = new ChangeTime(network);
        int acc = 0;

        for(int i = 0 ; i < 100 ; i++) {
            network._operator = changeTime;
            network._operator.propose();
            network.setDirty(true);
            double nextLikelihood = network.logDensity();
            System.out.println(nextLikelihood);

            double logAlpha = nextLikelihood - prevLikelihood;
            //logAlpha += network.splitting.logWeight;
            //logAlpha -= network.nextSplitting.logWeight;
            System.out.println(prevLikelihood + " " + nextLikelihood + " " + logAlpha);

            if (logAlpha >= Math.log(Randomizer.getRandomDouble())) {
                network.accept();
                acc++;
                prevLikelihood = nextLikelihood;
                System.out.println("Accepted");
            } else {
                network.undo();
                network.reject();
            }
        }

        System.out.println("Acceptance: " + acc);

    }


    static void test2taxa(String[] args) {
        double pi0 = 0.5;
        double pi1 = 1- pi0;

        BiAllelicGTR BAGTRModel = new BiAllelicGTR(new double[] {pi0, pi1}, new double[] {1.0/ (2.0 * pi0)});


        double alpha1 = 0.5;
        double alpha2 = 0.4;
        double alpha3 = 0.6;
        double alpha4 = 0.8;
        double a = 0;
        double b = 0;
        double c = 0;
        double theta = 0.1;
        double gamma = 0.55;

        //String netstring = String.format("[%f](((B:%f)I3#H1:%f::%f,A:%f)I1:%f,(C:%f,I3#H1:%f::%f)I2:%f)I0;", theta, alpha2, alpha1-alpha2, gamma, alpha1, alpha4-alpha1, alpha3, alpha3-alpha2, 1.0-gamma, alpha4-alpha3);
        String netstring = String.format("[%f](B:%f,A:%f)I0;", theta, alpha2, alpha2);

        //String netstring = "[0.02](((((C:0.01)#H1:0.02::0.6,B:0.03):0.01)#H2:0.02::0.7,A:0.06):0.03,((D:0.02,#H1:0.01::0.4):0.03,#H2:0.01::0.3):0.04);";
        //String netstring = "[0.03478971320729695](((((D:0.018773841591172784,C:0.018773841591172784):0.01665803647774957,B:0.035431878068922354):5.193091455444887E-4)#H1:0.042646185645307325::0.490471497495645,((#H1:0.01768017414602177::0.509528502504355,A:0.053631361360488614):2.981418077369069E-5)#H2:0.024936197318511863::0.28337456507014525):2.9602647955248584E-4,#H2:0.02523222379806435::0.7166254349298548);";
        //String netstring = "[0.07138571458425186]((((B:0.021545725769406308,(D:0.02072866724785186,C:0.02072866724785186):8.170585215544479E-4):0.012235537519323995)#H1:0.02015659627482603::0.17371727563554662,A:0.05393785956355633):0.005245366442141658,#H1:0.025401962716967687::0.8262827243644534);";
        Network<NetNodeInfo> trueNetwork = Networks.readNetwork("I0;");
        trueNetwork = Networks.readNetworkWithRootPop(netstring);
        for(Object nodeObj : trueNetwork.dfs()) {
            NetNode node = (NetNode) nodeObj;
            for(Object parentObj : node.getParents()) {
                NetNode parent = (NetNode) parentObj;
                node.setParentSupport(parent, trueNetwork.getRoot().getRootPopSize());
            }
        }

        System.out.println(Networks.getFullString(trueNetwork));
        Map<String, List<String>> species2alleles = new HashMap<>();
        species2alleles.put("I0", new ArrayList<>());
        species2alleles.get("I0").add("A");
        species2alleles.get("I0").add("B");
        species2alleles.get("I0").add("C");

        int numSites = 1000;
        boolean useOnlyPolymorphic = false;
        SimSNPInNetwork simulator = new SimSNPInNetwork(BAGTRModel, 12345678L);
        simulator._diploid = false;
        Map<String, String> onesnp = simulator.generateSNPs(trueNetwork, null, numSites, !useOnlyPolymorphic);
        //onesnp.clear();
        //onesnp.put("A", "0");
        //onesnp.put("B", "0");
        //onesnp.put("A", "00000000000000000000000000000000000000000000000000000000000000000000000000000000");
        //onesnp.put("B", "00000000000000000000000000000000000000000000000000000000000000000000000000000000");

        List<MarkerSeq> alns = new ArrayList<>();
        MarkerSeq aln = new MarkerSeq(onesnp);
        alns.add(aln);
        aln._RPatterns = SNAPPLikelihood.haploidSequenceToPatterns(null, alns);


        Utils._START_NET = Networks.getFullString(trueNetwork);
        Utils._POP_SIZE_MEAN = 0.1;
        Utils._ESTIMATE_POP_SIZE = false;
        Utils._ESTIMATE_POP_PARAM = false;

        Utils.DISABLE_TOPOLOGY_MOVES = true;

        Utils.SAMPLE_SPLITTING = true;

        Utils._NUM_THREADS = 1;
        Randomizer.setSeed(new Random().nextLong());

        UltrametricNetwork network = new UltrametricNetwork( Networks.getFullString(trueNetwork), alns, null, BAGTRModel);
        network.setDirty(true);
        ChangeInheritance changeInheritance = new ChangeInheritance(network);
        ChangeTime changeTime = new ChangeTime(network);
        ScaleRootTime scaleRootTime = new ScaleRootTime(network);
        network._operator = null;
        network.propose();
        double prevLikelihood = network.logDensity();
        network.accept();

        int acc = 0;
        Mean heightMean = new Mean();

        for(int i = 0 ; i < 1000 ; i++) {
            network._operator = null;
            double logHR = network.propose();
            network.setDirty(true);
            double nextLikelihood = network.logDensity();
            System.out.println(nextLikelihood);

            double logAlpha = nextLikelihood - prevLikelihood + logHR;
            if(Utils.SAMPLE_SPLITTING) {
                //logAlpha += network.splitting.logWeight;
                //logAlpha -= network.nextSplitting.logWeight;
            }
            System.out.println(prevLikelihood + " " + nextLikelihood + " " + logAlpha);

            if (logAlpha >= Math.log(Randomizer.getRandomDouble())) {
                network.accept();
                acc++;
                prevLikelihood = nextLikelihood;
                System.out.println("Accepted");
            } else {
                network.undo();
                network.reject();
                nextLikelihood = network.logDensity();
                if(nextLikelihood != prevLikelihood) throw new RuntimeException("!!!!!");
            }
            heightMean.increment(network.getNetwork().getRoot().getData().getHeight());
        }

        System.out.println("Acceptance: " + acc);
        System.out.println("Mean height: " + heightMean.getResult());
    }

    public static void main(String args[]) {
        test2taxa(args);
    }
}
