package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.ISMB2018;

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
 * Date: 4/29/18
 * Time: 12:19 PM
 * To change this template use File | Settings | File Templates.
 */
public class FigureConvergencePlot {
    public static void go() {
        double pi0 = 0.5;
        double pi1 = 1- pi0;

        BiAllelicGTR BAGTRModel = new BiAllelicGTR(new double[] {pi0, pi1}, new double[] {1.0/ (2.0 * pi0)});
        boolean useOnlyPolymorphic = false;
        int sizes[] = new int[]{100, 1000, 10000, 100000, 1000000};

        try {
            PrintWriter out = new PrintWriter("ConvergencePlot.txt");

            for(int numSites : sizes) {
                out.println("START");
                out.println(numSites);
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
                        out.println(difference);
                        mean.increment(difference);
                        standardDeviation.increment(difference);
                    }

                }
                out.println("END");
            }

            out.close();
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }

    }

    public static void main(String[] args) {
        go();
    }
}
