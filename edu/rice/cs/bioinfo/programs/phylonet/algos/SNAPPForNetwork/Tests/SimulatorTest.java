package edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.Tests;

import edu.rice.cs.bioinfo.library.programming.Tuple3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.distribution.SNAPPLikelihood;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.R;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.RPattern;
import edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPPForNetwork.SNAPPAlgorithm;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimSNPInNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.model.BiAllelicGTR;
import edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork.SuperNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.junit.Test;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 12/15/17
 * Time: 2:45 PM
 * To change this template use File | Settings | File Templates.
 */
public class SimulatorTest {
    @Test
    public void testPolyploid() {
        double pi0 = 0.9;
        double pi1 = 1- pi0;

        BiAllelicGTR BAGTRModel = new BiAllelicGTR(new double[] {pi0, pi1}, new double[] {1.0/ (2.0 * pi0)});
        boolean useOnlyPolymorphic = false;
        int numSites = 10000;

        Network trueNetwork = Networks.readNetworkWithRootPop("[0.01](((((C:0.0015:0.01)#H1:0.0045:0.01:0.5,B:0.006:0.01):0.003:0.01)#H2:0.0015:0.01:0.5,A:0.0105:0.01):0.0045:0.01,((D:0.003:0.01,#H1:0.0015:0.01:0.5):0.009:0.01,#H2:0.003:0.01:0.5):0.003:0.01);");

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
                mynode.setName("I" + nameCount);
                nameCount++;
            }
        }

        Map<Tuple3<String, String, String>, Map<RPattern, double[]>> triplets2patterns = new HashMap<>();

        SimSNPInNetwork simulator = new SimSNPInNetwork(BAGTRModel, 12345678L);
        simulator._diploid = false;
        simulator._polyploid = 3;
        Map<String, String> onesnp = simulator.generateSNPs(trueNetwork, null, numSites, !useOnlyPolymorphic);

        List<Alignment> alnwarp = new ArrayList<>();
        alnwarp.add(new Alignment(onesnp));

        Map<RPattern, double[]> patterns;
        patterns = SNAPPLikelihood.polyploidSequenceToPatterns(null, alnwarp, 3);

        for(String species : onesnp.keySet()) {
            System.out.println(species +  " " + onesnp.get(species));
        }

        for(RPattern pattern : patterns.keySet()) {
            //System.out.println(pattern + " " + patterns.get(pattern)[0] + " " + patterns.get(pattern)[1]);
            R.maxLineages = patterns.keySet().iterator().next().sumLineages();
            SNAPPAlgorithm run = new SNAPPAlgorithm(trueNetwork.clone(), BAGTRModel, null);
            double likelihood = 0;
            try {
                long start = System.currentTimeMillis();
                likelihood = run.getProbability(pattern);
            } catch(Exception e) {
                e.printStackTrace();
                System.out.println("Exceptional network");
            }
            System.out.println(pattern + " " + patterns.get(pattern)[0] + " " + likelihood * Math.exp(patterns.get(pattern)[1]) * numSites);
        }
    }

    @Test
    public void subnet() throws Exception {
        double pi0 = 0.9;
        double pi1 = 1- pi0;

        BiAllelicGTR BAGTRModel = new BiAllelicGTR(new double[] {pi0, pi1}, new double[] {1.0/ (2.0 * pi0)});
        boolean useOnlyPolymorphic = false;
        int numSites = 1000000;

        //Network trueNetwork = Networks.readNetwork("(((A:0.01:0.006,B:0.01:0.006):0.01:0.006,C:0.02:0.006):0.01:0.006,D:0.03:0.006);");
        //trueNetwork.getRoot().setRootPopSize(0.006);
        //Network trueNetwork = Networks.readNetwork("(((((Q:0.004:0.006)I5#H1:0.002:0.005:0.7,(A:0.003:0.006)I6#H2:0.003:0.005:0.6)I3:0.016:0.005,L:0.022:0.006)I2:0.02:0.005,(I5#H1:0.003:0.005:0.3,R:0.014:0.006)I4:0.028:0.005)I1:0.038:0.005,(C:0.005:0.006,I6#H2:0.002:0.005:0.4)I7:0.075:0.005);");
        //trueNetwork.getRoot().setRootPopSize(0.006);
        //Network trueNetwork = Networks.readNetwork("(((((Q:0.004:0.006)I5#H1:0.002:0.005:0.0,(A:0.003:0.006)I6#H2:0.003:0.005:0.0)I3:0.016:0.005,L:0.022:0.006)I2:0.02:0.005,(I5#H1:0.003:0.005:1.0,R:0.014:0.006)I4:0.028:0.005)I1:0.038:0.005,(C:0.005:0.006,I6#H2:0.002:0.005:1.0)I7:0.075:0.005);");
        //trueNetwork.getRoot().setRootPopSize(0.006);
        //Network trueNetwork = Networks.readNetwork("(((L:0.022:0.006)I2:0.02:0.005,((Q:0.004:0.006)I5:0.003:0.005,R:0.014:0.006)I4:0.028:0.005)I1:0.038:0.005,(C:0.005:0.006,(A:0.003:0.006)I6:0.002:0.005)I7:0.075:0.005);");
        //trueNetwork.getRoot().setRootPopSize(0.006);
        //Network trueNetwork = Networks.readNetwork("(((A:0.003:0.006)I6#H1:0.03900000000000001:0.006:0.0,R:0.042:0.006)I1:0.038:0.006,(C:0.005:0.006,I6#H1:0.002:0.006:1.0)I7:0.075:0.006)I0;");
        //trueNetwork.getRoot().setRootPopSize(0.006);
        //Network trueNetwork = Networks.readNetwork("((R:0.042:0.006)I1:0.038:0.006,(C:0.005:0.006,(A:0.003:0.006)I6:0.002:0.006)I7:0.075:0.006)I0;");
        //trueNetwork.getRoot().setRootPopSize(0.006);
        //Network trueNetwork = Networks.readNetworkWithRootPop("[0.01](A:0.19138876485031905:0.01,(B:0.12682431493516344:0.01,(C:0.09614402921753337:0.01,(((D:0.06864873769305778:0.01,((E:0.04784096272944867:0.01,(F:0.03053172961628988:0.01,G:0.03053172961628988:0.01):0.017309233113158788:0.01):0.007107765870675464:0.01,(H:0.04381981472730136:0.01,I:0.04381981472730136:0.01):0.011128913872822777:0.01):0.013700009092933646:0.01):0.003126652986989553:0.01,(J:0.06689143414599366:0.01,K:0.06689143414599366:0.01):0.004883956534053671:0.01):0.0016928336785998616:0.01,((L:0.044882393305892485:0.01,M:0.044882393305892485:0.01):0.022226283989763562:0.01,(N:0.06296395354420453:0.01,(O:0.03361286516308759:0.01,P:0.03361286516308759:0.01):0.02935108838111694:0.01):0.004144723751451515:0.01):0.006359547062991147:0.01):0.022675804858886178:0.01):0.030680285717630068:0.01):0.06456444991515561:0.01);");
        Network trueNetwork = Networks.readNetworkWithRootPop("[0.01](A:0.20138876485031905:0.01,(B:0.13682431493516344:0.01,((C:0.03614402921753337:0.01,(I:0.005:0.01)#H1:0.03114402921753337:0.01:0.3):0.07:0.01,((((D:0.02:0.01,(L:0.015:0.01)#H3:0.005:0.01:0.2):0.04864873769305778:0.01,((E:0.04784096272944867:0.01,((F:0.01053172961628988:0.01,#H1:0.00553172961628988:0.01:0.7):0.02:0.01,G:0.03053172961628988:0.01):0.017309233113158788:0.01):0.007107765870675464:0.01,H: 0.0549487286:0.01):0.013700009092933646:0.01):0.013126652986989553:0.01,((J:0.03:0.01, (O:0.01:0.01)#H2:0.02:0.01:0.4):0.02689143414599366:0.01,K:0.05689143414599366:0.01):0.024883956534053671:0.01):0.0116928336785998616:0.01,((#H3:0.019882393305892485:0.01:0.8,M:0.044882393305892485:0.01):0.032226283989763562:0.01,(N:0.06296395354420453:0.01,(#H2:0.01361286516308759:0.01:0.6,P:0.02361286516308759:0.01):0.03935108838111694:0.01):0.014144723751451515:0.01):0.016359547062991147:0.01):0.012675804858886178:0.01):0.030680285717630068:0.01):0.06456444991515561:0.01);");

        Map<Tuple3<String, String, String>, Map<RPattern, double[]>> triplets2patterns = new HashMap<>();

        SimSNPInNetwork simulator = new SimSNPInNetwork(BAGTRModel, 12345678L);
        simulator._diploid = false;
        Map<String, String> onesnp = simulator.generateSNPs(trueNetwork, null, numSites, !useOnlyPolymorphic);

        List<Alignment> alns = new ArrayList<>();
        Alignment aln = new Alignment(onesnp);
        alns.add(aln);
        List<Tuple3<String, String, String>> triplets = new ArrayList<>();

        List<String> taxa = alns.get(0).getTaxaNames();
        for(int i = 0 ; i < taxa.size() ; i++) {
            for(int j = i + 1 ; j < taxa.size() ; j++) {
                for(int k = j + 1 ; k < taxa.size() ; k++) {
                    Tuple3<String, String, String> triplet = new Tuple3<>(taxa.get(i), taxa.get(j), taxa.get(k));
                    triplets.add(triplet);

                    Map<String, String> sequence = new HashMap<>();
                    sequence.put(triplet.Item1, alns.get(0).getAlignment().get(triplet.Item1));
                    sequence.put(triplet.Item2, alns.get(0).getAlignment().get(triplet.Item2));
                    sequence.put(triplet.Item3, alns.get(0).getAlignment().get(triplet.Item3));
                    List<Alignment> alnwarp = new ArrayList<>();
                    alnwarp.add(new Alignment(sequence));

                    Map<RPattern, double[]> patterns;
                    patterns = SNAPPLikelihood.haploidSequenceToPatterns(null, alnwarp);

                    triplets2patterns.put(triplet, patterns);

                    System.out.println(triplet.Item1 + " " + triplet.Item2 + " " + triplet.Item3);

                }
            }
        }

        Map<Tuple3<String, String, String>, Network> triplets2subnetworks = new HashMap<>();

        for(Tuple3<String, String, String> triplet : triplets) {
            SuperNetwork superNetwork = new SuperNetwork(new ArrayList<>());

            List<String> leaves = new ArrayList<>();
            leaves.add(triplet.Item1);
            leaves.add(triplet.Item2);
            leaves.add(triplet.Item3);

            Network subNetwork = superNetwork.getSubNetwork(trueNetwork, leaves, true, false);
            Networks.autoLabelNodes(subNetwork);
            triplets2subnetworks.put(triplet, subNetwork);
        }

        Mean mean = new Mean();
        StandardDeviation standardDeviation = new StandardDeviation();

        for(Tuple3<String, String, String> triplet : triplets) {
            Network network = triplets2subnetworks.get(triplet);
            System.out.println(network);

            int nameCount = 0;
            for (Object node : network.dfs()) {
                NetNode mynode = (NetNode) node;
                if (mynode.getName().equals("")) {
                    mynode.setName("I" + nameCount);
                    nameCount++;
                }
            }

            // simulate sites for a subnet
            // compare to the portion of pattern with full data
            /*Map<String, String> onesnp1 = simulator.generateSNPs(network, null, numSites, !useOnlyPolymorphic);

            List<Alignment> alns1 = new ArrayList<>();
            Alignment aln1 = new Alignment(onesnp1);
            alns1.add(aln1);
            aln1._RPatterns = SNAPPLikelihood.haploidSequenceToPatterns(null, alns1);

            for (RPattern pattern : aln1._RPatterns.keySet()) {
                System.out.println(pattern + " " + aln1._RPatterns.get(pattern)[0] + "     " + triplets2patterns.get(triplet).get(pattern)[0]);
            }

            System.out.println();*/

            // compute expected number of this pattern for a subnet
            // compare to the portion of pattern with full data
            for (RPattern pattern : triplets2patterns.get(triplet).keySet()) {
                R.maxLineages = triplets2patterns.get(triplet).keySet().iterator().next().sumLineages();
                SNAPPAlgorithm run = new SNAPPAlgorithm(triplets2subnetworks.get(triplet), BAGTRModel, null);
                double likelihood = 0;
                try {
                    long start = System.currentTimeMillis();
                    likelihood = run.getProbability(pattern);
                } catch(Exception e) {
                    e.printStackTrace();
                    System.out.println("Exceptional network");
                }
                System.out.println(pattern + " " + triplets2patterns.get(triplet).get(pattern)[0] + " " + likelihood);
                double difference = triplets2patterns.get(triplet).get(pattern)[0] / numSites - likelihood;
                mean.increment(difference);
                standardDeviation.increment(difference);
            }

        }
        System.out.println("Mean: " + mean.getResult() + " SD: " + standardDeviation.getResult());
    }
}
