package edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.Tests;

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.Tuple3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.distribution.SNAPPLikelihood;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.alignment.MarkerSeq;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimSNPInNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.BiAllelicGTR;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import org.junit.Test;

import java.util.*;

import static org.junit.Assert.assertEquals;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 11/14/17
 * Time: 2:43 PM
 * To change this template use File | Settings | File Templates.
 */
public class LikelihoodTest {
    @Test
    public void testNetwork1() throws Exception {
        Tuple<String, String> a = new Tuple<>("aa", "bb");
        Tuple<String, String> b = new Tuple<>("aa", "dd");
        System.out.println(a == b);
        System.out.println(a.equals(b));

        Tuple3<String, String, String> c = new Tuple3<>("aa", "bb", "cc");
        Tuple3<String, String, String> d = new Tuple3<>("aa", "bb", "dd");
        System.out.println(c == d);
        System.out.println(c.equals(d));
    }

    @Test
    public void testAnalytical() throws Exception {
        double pi0 = 0.5;
        double pi1 = 1- pi0;

        BiAllelicGTR BAGTRModel = new BiAllelicGTR(new double[] {pi0, pi1}, new double[] {1.0/ (2.0 * pi0)});
        boolean useOnlyPolymorphic = false;
        int numSites = 1000;

        double alpha1 = 0.5;
        double alpha2 = 0.4;
        double alpha3 = 0.6;
        double alpha4 = 0.8;
        double a = 0;
        double b = 0;
        double c = 0;
        double theta = 0.1;
        double gamma = 0.5;

        //String netstring = String.format("[%f]((((C:%f)I3#H1:%f::%f,B:%f)I2:%f,A:%f)I1:%f,I3#H1:%f::%f)I0;", theta,alpha3,alpha2-alpha3, gamma, alpha2, alpha1-alpha2, alpha1, alpha4-alpha1, alpha4-alpha3, 1.0-gamma);
        String netstring = String.format("[%f](((B:%f)I3#H1:%f::%f,A:%f)I1:%f,(C:%f,I3#H1:%f::%f)I2:%f)I0;", theta, alpha2, alpha1-alpha2, gamma, alpha1, alpha4-alpha1, alpha3, alpha3-alpha2, 1.0-gamma, alpha4-alpha3);

        Network<Object> trueNetwork;
        trueNetwork = Networks.readNetworkWithRootPop(netstring);

        for(Object nodeObj : trueNetwork.dfs()) {
            NetNode node = (NetNode) nodeObj;
            for(Object parentObj : node.getParents()) {
                NetNode parent = (NetNode) parentObj;
                node.setParentSupport(parent, trueNetwork.getRoot().getRootPopSize());
            }
        }

        SimSNPInNetwork simulator = new SimSNPInNetwork(BAGTRModel, new Random().nextLong());
        simulator._diploid = false;
        Map<String, String> onesnp = simulator.generateSNPs(trueNetwork, null, numSites, !useOnlyPolymorphic);

        List<MarkerSeq> alns = new ArrayList<>();
        MarkerSeq aln = new MarkerSeq(onesnp);
        alns.add(aln);
        aln._RPatterns = SNAPPLikelihood.haploidSequenceToPatterns(null, alns);

        double sum = 0.0;

        for(RPattern pattern : aln._RPatterns.keySet()) {
            Network cloneNetwork = trueNetwork.clone();
            R.maxLineages = aln._RPatterns.keySet().iterator().next().sumLineages();
            SNAPPAlgorithm run = new SNAPPAlgorithm(cloneNetwork, BAGTRModel, null);
            double likelihood = 0;
            try {
                long start = System.currentTimeMillis();
                likelihood = run.getProbability(pattern, null);
                sum += likelihood;
                System.out.println(pattern + " " + likelihood * Math.exp(aln._RPatterns.get(pattern)[1]) * numSites  + " " + aln._RPatterns.get(pattern)[0] );
                System.out.println("Time: " + (System.currentTimeMillis()-start)/1000.0);
            } catch(Exception e) {
                e.printStackTrace();
                System.out.println("Exceptional network");
            }
        }

        System.out.println("Sum = " + sum);

    }

    @Test
    public void testHaploid() throws Exception {
        double pi0 = 0.9;
        double pi1 = 1- pi0;

        BiAllelicGTR BAGTRModel = new BiAllelicGTR(new double[] {pi0, pi1}, new double[] {1.0/ (2.0 * pi0)});
        boolean useOnlyPolymorphic = false;
        int numSites = 1000;

        Network<Object> trueNetwork;
        trueNetwork = Networks.readNetwork("(((A:0.009342197568949029:0.023422567184402446,(R:1.9855469004628661E-4:0.023422567184402446)II0#H1:0.009143642878902743:0.023422567184402446:0.09961893620579232)I0:0.024426689411379966:0.023422567184402446,II0#H1:0.03357033229028271:0.023422567184402446:0.9003810637942077)I0:0.04212588078922362:0.023422567184402446,C:0.07589476776955262:0.023422567184402446)I1;");
        trueNetwork.getRoot().setRootPopSize(0.023422567184402446);

        int nameCount = 0;
        for(Object node : trueNetwork.dfs()) {
            NetNode mynode = (NetNode) node;
            if(mynode.getName().equals("")) {
                mynode.setName("I" + nameCount);
                nameCount++;
            }
        }

        Utils._NUM_THREADS = 1;
        Utils._CONST_POP_SIZE = false;
        SimSNPInNetwork simulator = new SimSNPInNetwork(BAGTRModel, new Random().nextLong());
        simulator._diploid = false;
        Map<String, String> onesnp = simulator.generateSNPs(trueNetwork, null, numSites, !useOnlyPolymorphic);
        //onesnp = new HashMap<>();
        //onesnp.put("R", "1001011010101011001101010001010111001110011001100101011111001011110000001010101001100000100100001001011000111011011101001010111101000010100100011100011001100011101111000001101111100101100001100101001100100010010000011100010100000111101101110100100010001100011111000101101101110000101100011101111000001110011100001101100101111010000001111000111010011010001110100011000101110111101100110011101010100100011100000011110101001111111000010011011011111011011000100100101000100001010101101111011011011001001011011111100110001111111111001011000010010110111100010001101001001001010101110110110000010100110010001101001110111001010010101100101000100100111101110000010010110100001000010001011101010110010010010011001010000011011001011111000100000010111100101000000000111000101110111110110101110010000101000011010000110011011111001001100101110100001110000100010101111101010000101010111011000101010010100001110101101100111010001000011011000100010101111000011000101010111000000100011011100100100100111110000100000110");
        //onesnp.put("Q", "1001011010101011001001010101010111001010011001100101111111001011110000001010001001100000110100001011011000111001011101101010111101000010100100011100011001100111101111000001100111101101000001001101001100110010010100001100011100000111001101110100100010001101011111000101101101110000101100011101110010001110011100001101100101111011000001111010111010011010001110100011000101110111101100110111101010100100011100000111110101001111111000011011011010111011011000000100101000100001010101101101011111111001001011011111100110001111111110001011000010010110111100011001101011001001010110111110111010010100110010001101001110011000010010101100101000100100111111110000010010110100001000000001011101000110010010010000011000001011011001011111000100000010111100101000000000111001101010111110100111110010000101100011010000110011010111001001100101110100001110000100010101111101010000101010101011000101000010100001110101101000111011000000011011000100010001111000011000101011111000000100011011100100100100111110000100000110");
        //onesnp.put("A", "1001011010101011001000010101010111001010011001100101111011000011111000001010001001100000110100001011001000111001011101101010111001000010100100011110011001100111101111000001101111101101000001001101001100110010010100001100011100000111001101110100100010001101011111000101100101110000101100010101110010001110011100001101100101111010000001111010111010011110001110100011000111110111101100110111101000101100011100000011110101001111111000010011011010111011011000000100101000100001010101101101010111101001001011011111100110001111111110001011000011010110111100010001101011001001010111111110111010010100010010001101001110011000110010001100101000100100111111110000010010110100001000000001011101010110000010010001011000001011011001011111000100000010111100101000000000111001101010111110100111110010000101100011010000110011010101001001100101110100001110010100010101111101010000101010101011000101000010100001110101101000111011000000011011000100010001111000001000101011111000000100011011100100100100111110000100000110");

        List<MarkerSeq> alns = new ArrayList<>();
        MarkerSeq aln = new MarkerSeq(onesnp);
        alns.add(aln);
        aln._RPatterns = SNAPPLikelihood.haploidSequenceToPatterns(null, alns);

        System.out.println(SNAPPLikelihood.computeSNAPPLikelihood(trueNetwork, aln._RPatterns, BAGTRModel));

    }

    @Test
    public void testDiploid() throws Exception {
        double pi0 = 0.9;
        double pi1 = 1- pi0;

        BiAllelicGTR BAGTRModel = new BiAllelicGTR(new double[] {pi0, pi1}, new double[] {1.0/ (2.0 * pi0)});
        boolean useOnlyPolymorphic = false;
        int numSites = 1000000;

        Network<Object> trueNetwork;

        //trueNetwork = Networks.readNetwork("((C:0.09604283044864592:0.006,G:0.09604283044864592:0.036)II1:0.08821096143809048:0.006,(R:0.16796441222747166:0.036,(L:0.13918714276096708:0.036,(A:0.10079837440250675:0.006,Q:0.10079837440250675:0.006)II4:0.038388768358460335:0.036)II3:0.02877726946650458:0.036)II2:0.016289379659264747:0.036)II0;");
        //trueNetwork.getRoot().setRootPopSize(0.036);
        //trueNetwork = Networks.readNetwork("(((G:0.03271098551352612,(C:0.014581262376039935,(A:0.007848436317467452)#H2:0.006732826058572483::0.20092242303662247):0.01812972313748619):6.589767080962955E-4)#H1:0.056974310762192545::0.973670280875714,(((Q:0.01575520654534508)#H3:0.008437888036381732::0.29205403332422386,R:0.02419309458172681):0.04760724352481251,(#H1:0.012622246801820981::0.026329719124286055,((#H3:1.8671708033869347E-4::0.7079459666757761,#H2:0.00809348730821632::0.7990775769633776):0.01963878547064636,L:0.03558070909633013):0.010411499927113266):0.025808129083095925):0.01854393487727564);");
        //trueNetwork.getRoot().setRootPopSize(0.036);
        trueNetwork = Networks.readNetwork("((((((((Q:0.007706894443703847:0.036)#H3:0.011286322972969953:0.036:0.23227570802381947,R:0.0189932174166738:0.036):4.832779408627669E-6:0.036)#H2:0.006601622173161768:0.036:0.8787824364801154,L:0.025599672369244195:0.036):0.005238671994843097:0.036,(#H3:0.016537156753233577:0.036:0.7677242919761805,A:0.024244051196937424:0.036):0.006594293167149868:0.036):0.004030047262403039:0.036)#H1:0.04499887208586771:0.036:0.6757746276842883,(#H2:0.03928130470143773:0.036:0.12121756351988455,#H1:0.023410963271029823:0.036:0.3242253723157117):0.021587908814837888:0.036):0.007571037854375656:0.036,(G:0.03494992023933711:0.036,C:0.03494992023933711:0.036):0.05248838132739659:0.036);");
        trueNetwork.getRoot().setRootPopSize(0.036);

        int nameCount = 0;
        for(Object node : trueNetwork.dfs()) {
            NetNode mynode = (NetNode) node;
            if(mynode.getName().equals("")) {
                mynode.setName("I" + nameCount);
                nameCount++;
            }
        }

        SimSNPInNetwork simulator = new SimSNPInNetwork(BAGTRModel, new Random().nextLong());
        simulator._diploid = true;
        Map<String, String> onesnp = simulator.generateSNPs(trueNetwork, null, numSites, !useOnlyPolymorphic);

        List<MarkerSeq> alns = new ArrayList<>();
        MarkerSeq aln = new MarkerSeq(onesnp);
        alns.add(aln);
        aln._RPatterns = SNAPPLikelihood.diploidSequenceToPatterns(null, alns);

        Algorithms.SWITCH_APPROX_SPLIT = false;

        int count = 0;
        double sum = 0.0;

        for(RPattern pattern : aln._RPatterns.keySet()) {
            //if(count > 10) break;
            count++;
            Network cloneNetwork = Networks.readNetwork(trueNetwork.toString());
            cloneNetwork.getRoot().setRootPopSize(trueNetwork.getRoot().getRootPopSize());
            R.maxLineages = aln._RPatterns.keySet().iterator().next().sumLineages();
            SNAPPAlgorithm run = new SNAPPAlgorithm(cloneNetwork, BAGTRModel, null);
            double likelihood = 0;
            try {
                long start = System.currentTimeMillis();
                likelihood = run.getProbability(pattern, null);
                sum += likelihood;
                System.out.println(pattern + " " + likelihood * Math.exp(aln._RPatterns.get(pattern)[1]) * numSites  + " " + aln._RPatterns.get(pattern)[0] );
                System.out.println("Time: " + (System.currentTimeMillis()-start)/1000.0);
            } catch(Exception e) {
                e.printStackTrace();
                System.out.println("Exceptional network");
            }
        }
        System.out.println("Sum " + sum);

        assertEquals(sum, 1.0, 1e-4);

        return;
    }

    @Test
    public void testPseudoLikelihood() throws Exception {
        double pi0 = 0.9;
        double pi1 = 1- pi0;

        BiAllelicGTR BAGTRModel = new BiAllelicGTR(new double[] {pi0, pi1}, new double[] {1.0/ (2.0 * pi0)});
        boolean useOnlyPolymorphic = false;
        int numSites = 1000000;

        Network<Object> trueNetwork;

        trueNetwork = Networks.readNetwork("((C:0.09604283044864592:0.036,G:0.09604283044864592:0.036)II1:0.08821096143809048:0.036,(R:0.16796441222747166:0.036,(L:0.13918714276096708:0.036,(A:0.10079837440250675:0.036,Q:0.10079837440250675:0.036)II4:0.038388768358460335:0.036)II3:0.02877726946650458:0.036)II2:0.016289379659264747:0.036)II0;");
        trueNetwork.getRoot().setRootPopSize(0.036);
        //trueNetwork = Networks.readNetwork("(((G:0.03271098551352612,(C:0.014581262376039935,(A:0.007848436317467452)#H2:0.006732826058572483::0.20092242303662247):0.01812972313748619):6.589767080962955E-4)#H1:0.056974310762192545::0.973670280875714,(((Q:0.01575520654534508)#H3:0.008437888036381732::0.29205403332422386,R:0.02419309458172681):0.04760724352481251,(#H1:0.012622246801820981::0.026329719124286055,((#H3:1.8671708033869347E-4::0.7079459666757761,#H2:0.00809348730821632::0.7990775769633776):0.01963878547064636,L:0.03558070909633013):0.010411499927113266):0.025808129083095925):0.01854393487727564);");
        //trueNetwork.getRoot().setRootPopSize(0.036);
        //trueNetwork = Networks.readNetwork("((((((((Q:0.007706894443703847)#H3:0.011286322972969953::0.23227570802381947,R:0.0189932174166738):4.832779408627669E-6)#H2:0.006601622173161768::0.8787824364801154,L:0.025599672369244195):0.005238671994843097,(#H3:0.016537156753233577::0.7677242919761805,A:0.024244051196937424):0.006594293167149868):0.004030047262403039)#H1:0.04499887208586771::0.6757746276842883,(#H2:0.03928130470143773::0.12121756351988455,#H1:0.023410963271029823::0.3242253723157117):0.021587908814837888):0.007571037854375656,(G:0.03494992023933711,C:0.03494992023933711):0.05248838132739659);");
        //trueNetwork.getRoot().setRootPopSize(0.036);
        trueNetwork = Networks.readNetwork("(((((G:3.723217345458331:0.036,H:3.723217345458331:0.036):5.131391932345352:0.036)#H1:23.578573233736876:0.036:0.5,((K:29.30932864678949:0.036,C:29.30932864678949:0.036):1.842046384485137:0.036,((((I:2.2260388160275455:0.036)#H2:14.637597073878071:0.036:0.5,J:16.863635889905616:0.036):6.487530888001629:0.036,(F:20.767587832150696:0.036,A:20.767587832150696:0.036):2.5835789457565497:0.036):2.5939521119047697:0.036,M:25.945118889812015:0.036):5.20625614146261:0.036):1.2818074802659325:0.036):0.6945919360420447:0.036,(((((#H2:6.433379566176971:0.036:0.5,L:8.659418382204516:0.036):0.9589721989891125:0.036)#H4:0.40313261644486253:0.036:0.5)#H3:3.185882965808421:0.036:0.5,#H4:3.5890155822532837:0.036:0.5):4.756014611940664:0.036,E:17.963420775387576:0.036):15.164353672195027:0.036):2.2165211933435387:0.036,(B:17.465139861448325:0.036,(#H1:7.473121360018542:0.036:0.5,((#H3:0.8821470719568634:0.036:0.5,N:10.903670269595354:0.036):2.7154673994911676:0.036,(D:5.4062591350430385:0.036,O:5.4062591350430385:0.036):8.212878534043483:0.036):2.7085929687357027:0.036):1.1374092236261006:0.036):17.879155779477816:0.036);");
        trueNetwork.getRoot().setRootPopSize(0.036);

        //trueNetwork = (new SuperNetwork(new ArrayList<>())).genRandomNetwork(15, 4);

        int nameCount = 0;
        for(Object node : trueNetwork.dfs()) {
            NetNode mynode = (NetNode) node;
            if(mynode.getName().equals("")) {
                mynode.setName("I" + nameCount);
                nameCount++;
            }
        }

        Utils._NUM_THREADS = 8;
        Utils._CONST_POP_SIZE = false;
        Algorithms.SWITCH_APPROX_SPLIT = true;
        SimSNPInNetwork simulator = new SimSNPInNetwork(BAGTRModel, 12345678L);
        simulator._diploid = true;
        Map<String, String> onesnp = simulator.generateSNPs(trueNetwork, null, numSites, !useOnlyPolymorphic);

        List<MarkerSeq> alns = new ArrayList<>();
        MarkerSeq aln = new MarkerSeq(onesnp);
        alns.add(aln);
        aln._RPatterns = SNAPPLikelihood.diploidSequenceToPatterns(null, alns);

        long start = System.currentTimeMillis();
        //System.out.println(SNAPPLikelihood.computeSNAPPLikelihood(trueNetwork, aln._RPatterns, BAGTRModel));
        System.out.println("Time of computing full: " + (System.currentTimeMillis()-start)/1000.0);


        start = System.currentTimeMillis();
        SNAPPPseudoLikelihood pseudoLikelihood = new SNAPPPseudoLikelihood(null, alns, simulator._diploid);
        System.out.println("Time of constructor: " + (System.currentTimeMillis()-start)/1000.0);
        start = System.currentTimeMillis();
        double l = pseudoLikelihood.computeSNAPPPseudoLogLikelihood(trueNetwork, null, BAGTRModel);
        System.out.println(l);
        System.out.println("Time of computing: " + (System.currentTimeMillis()-start)/1000.0);

        int count = 0;
        double sum = 0.0;

        for(RPattern pattern : aln._RPatterns.keySet()) {
            //if(count > 10) break;
            count++;
            Network cloneNetwork = Networks.readNetwork(trueNetwork.toString());
            cloneNetwork.getRoot().setRootPopSize(trueNetwork.getRoot().getRootPopSize());
            R.maxLineages = aln._RPatterns.keySet().iterator().next().sumLineages();
            SNAPPPseudoLikelihood run = new SNAPPPseudoLikelihood(cloneNetwork, BAGTRModel, null);
            double likelihood = 0;
            try {
                start = System.currentTimeMillis();
                likelihood = run.getProbability(pattern);
                sum += likelihood;
                System.out.println(pattern + " " + likelihood * Math.exp(aln._RPatterns.get(pattern)[1]) * numSites  + " " + aln._RPatterns.get(pattern)[0] );
                System.out.println("Time: " + (System.currentTimeMillis()-start)/1000.0);
            } catch(Exception e) {
                e.printStackTrace();
                System.out.println("Exceptional network");
            }
        }
        System.out.println("Sum " + sum);

        assertEquals(sum, 1.0, 1e-4);

        return;
    }

    @Test
    public void testLikelihood() {
        double alpha1 = 0.005;
        double alpha2 = 0.004;
        double alpha3 = 0.006;
        double alpha4 = 0.008;
        double theta = 0.01;
        double gamma = 0.6;
        double pi0 = 0.5;
        double pi1 = 1 - pi0;

        Map<String, List<String>> species2alleles = new HashMap<>();
        Map<String, String> alleles2species = new HashMap<>();

        species2alleles.put("A", new ArrayList<>());
        species2alleles.put("B", new ArrayList<>());
        species2alleles.put("C", new ArrayList<>());
        for (int i = 0; i < 1; i++) {
            species2alleles.get("A").add("A_" + i);
            species2alleles.get("B").add("B_" + i);
            species2alleles.get("C").add("C_" + i);
            alleles2species.put("A_" + i, "A");
            alleles2species.put("B_" + i, "B");
            alleles2species.put("C_" + i, "C");
        }

        int numSites = 10;

        BiAllelicGTR BAGTRModel = new BiAllelicGTR(new double[]{pi0, pi1}, new double[]{1.0 / (2.0 * pi0)});

        String netstring = String.format("[%f](((B:%f)I3#H1:%f::%f,A:%f)I1:%f,(C:%f,I3#H1:%f::%f)I2:%f)I0;", theta, alpha2, alpha1 - alpha2, gamma, alpha1, alpha4 - alpha1, alpha3, alpha3 - alpha2, 1.0 - gamma, alpha4 - alpha3);

        Network trueNetwork = Networks.readNetworkWithRootPop(netstring);
        for (Object nodeObj : Networks.postTraversal(trueNetwork)) {
            NetNode node = (NetNode) nodeObj;
            for (Object parentObj : node.getParents()) {
                NetNode parent = (NetNode) parentObj;
                node.setParentSupport(parent, theta);
            }
        }

        SimSNPInNetwork simulator = new SimSNPInNetwork(BAGTRModel, 12345678L);
        Map<String, String> snps = simulator.generateSNPs(trueNetwork, null, numSites, true);

        List<MarkerSeq> alns = new ArrayList<>();
        MarkerSeq aln = new MarkerSeq(snps);
        alns.add(aln);
        aln._RPatterns = SNAPPLikelihood.haploidSequenceToPatterns(null, alns);

        System.out.println(SNAPPLikelihood.computeSNAPPLikelihood(trueNetwork, aln._RPatterns, BAGTRModel));

    }
}
