package edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork;

import com.google.gson.Gson;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.Tuple3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.core.MC3Core;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.NetNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.structs.UltrametricNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.summary.SummaryBL;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.*;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 5/26/18
 * Time: 9:41 AM
 * To change this template use File | Settings | File Templates.
 */
public class Pipeline {
    public static final boolean printDetails_ = true;

    public static void initNetHeights(Network<NetNodeInfo> network) {
        for(NetNode<NetNodeInfo> node : Networks.postTraversal(network)) {
            if(node.getData() == null) {
                node.setData(new NetNodeInfo(0.0));
            }
            for(NetNode<NetNodeInfo> par : node.getParents()) {
                double dist = node.getParentDistance(par);
                if(par.getData() == null) {
                    par.setData(new NetNodeInfo(node.getData().getHeight() + dist));
                }
            }
        }
    }

    public static Map<String, Map<String, String>> randomSelectLoci(Map<String, Map<String, String>> multilocusData, int subsetLociNumber, int curN) {
        Random random = new java.util.Random(Utils._SEED/* + curN*/);
        List<String> names = new ArrayList<>(multilocusData.keySet());
        Collections.shuffle(names, random);

        Map<String, Map<String, String>> newdata = new HashMap<>();
        for(int i = 0 ; i < subsetLociNumber ; i++) {
            //System.out.println(names.get(i));
            newdata.put(names.get(i), multilocusData.get(names.get(i)));
        }

        return newdata;
    }

    // stage 1: divide the input data into combinations of 3-species
    // and run on desired combination
    public static void stage1(Map<String, Map<String, String>> multilocusData, Map<String, List<String>> taxonMap, int number, int subsetLociNumber, String outgroup, String gtoutgroup) {
        int nTaxa = taxonMap.size();
        int index = 0;
        List<String> taxa = new ArrayList<>();
        List<String> selectedTaxa = new ArrayList<>();
        taxa.addAll(taxonMap.keySet());

        if(outgroup != null) {
            nTaxa--;
            taxa.remove(outgroup);
        }

        Collections.sort(taxa);
        for(int i1 = 0 ; i1 < nTaxa ; i1++) {
            for(int i2 = i1 + 1 ; i2 < nTaxa ; i2++) {
                for(int i3 = i2 + 1 ; i3 < nTaxa ; i3++) {
                    if(index == number) {
                        selectedTaxa.add(taxa.get(i1));
                        selectedTaxa.add(taxa.get(i2));
                        selectedTaxa.add(taxa.get(i3));
                        i1 = nTaxa;
                        i2 = nTaxa;
                        break;
                    }
                    index++;
                }
            }
        }

        if(selectedTaxa.size() != 3) {
            throw new RuntimeException("Wrong number: " + number);
        }

        if(outgroup != null) {
            selectedTaxa.add(outgroup);
        }

        Map<String, Map<String, String>> subsetData = new HashMap<>();
        for(String locus : multilocusData.keySet()) {
            subsetData.put(locus, new HashMap<>());
            for(String taxon : selectedTaxa) {
                for(String allele : taxonMap.get(taxon)) {
                    subsetData.get(locus).put(allele, multilocusData.get(locus).get(allele));
                }
            }
        }

        if(subsetLociNumber > 0 && subsetLociNumber < subsetData.size())
            subsetData = randomSelectLoci(subsetData, subsetLociNumber, number);

        Map<String, List<String>> subsetTaxonMap = new HashMap<>();
        for(String taxon : selectedTaxa) {
            subsetTaxonMap.put(taxon, taxonMap.get(taxon));
        }


        if(Utils._START_GT_LIST != null) {

            for(String species : subsetTaxonMap.keySet()) {
                if(subsetTaxonMap.get(species).contains(gtoutgroup)) {
                    subsetTaxonMap.get(species).remove(gtoutgroup);
                    for (String locus : subsetData.keySet()) {
                        subsetData.get(locus).remove(gtoutgroup);
                    }
                }
            }

            for(int i = 0 ; i < Utils._START_GT_LIST.size() ; i++) {
                List<String> toRemove = new ArrayList<>();
                for (String locus : Utils._START_GT_LIST.get(i).keySet()) {
                    if (!subsetData.containsKey(locus)) {
                        toRemove.add(locus);
                    }
                }

                for(String locus : toRemove) {
                    Utils._START_GT_LIST.get(i).remove(locus);
                }
            }

            for(int i = 0 ; i < Utils._START_GT_LIST.size() ; i++) {
                for (String locus : Utils._START_GT_LIST.get(i).keySet()) {
                    STITree stitree = new STITree(Trees.readTree(Utils._START_GT_LIST.get(i).get(locus)));
                    List<String> leaves = new ArrayList<>();
                    for(String species : subsetTaxonMap.keySet()) {
                        for(String allele : subsetTaxonMap.get(species)) {
                            leaves.add(allele);
                        }
                    }

                    stitree.rerootTreeAtEdge(gtoutgroup);
                    stitree.removeNode(gtoutgroup);
                    Trees.removeBinaryNodes(stitree);
                    stitree.constrainByLeaves(leaves);

                    for(Object nodeObj : stitree.postTraverse()) {
                        STINode node = (STINode) nodeObj;
                        node.setParentDistance(TNode.NO_DISTANCE);
                    }

                    Utils._START_GT_LIST.get(i).put(locus, stitree.toNewick());
                }
            }

        }

//        if(Utils.ONLY_BACKBONE_OP) {
//            int maxReti = Utils._NET_MAX_RETI;
//            long burnin = Utils._BURNIN_LEN;
//            long chainlen = Utils._CHAIN_LEN;
//            List<Double> mc3_chains = Utils._MC3_CHAINS;
//
//            Utils.ONLY_BACKBONE_OP = false;
//            Utils._BURNIN_LEN = 2 * Utils._SAMPLE_FREQUENCY;
//            Utils._CHAIN_LEN = 5 * Utils._SAMPLE_FREQUENCY;
//            Utils._MC3_CHAINS = null;
//            Utils._NET_MAX_RETI = 0;
//
//            MC3Core mc3 = new MC3Core(alignments);
//            mc3.run();
//
//            Utils._START_NET = new ArrayList<>();
//            Utils._START_NET.add(mc3.getLastNetwork());
//
//            Utils.ONLY_BACKBONE_OP = true;
//            Utils._BURNIN_LEN = burnin;
//            Utils._CHAIN_LEN = chainlen;
//            Utils._MC3_CHAINS = mc3_chains;
//            Utils._NET_MAX_RETI = maxReti;
//        }
        Utils._TAXON_MAP = subsetTaxonMap;
        List<Alignment> alignments = new ArrayList<>();
        for(String key : subsetData.keySet()) {
            alignments.add(new Alignment(subsetData.get(key), key));
        }

        Collections.sort(alignments);
        MC3Core mc3 = new MC3Core(alignments);
        mc3.run();
    }

    public static void stage1_5(Map<String, Map<String, String>> multilocusData, Map<String, List<String>> taxonMap, int number) {
        int nTaxa = taxonMap.size();
        int index = 0;
        List<String> taxa = new ArrayList<>();
        List<String> selectedTaxa = new ArrayList<>();
        taxa.addAll(taxonMap.keySet());
        Collections.sort(taxa);
        for(int i1 = 0 ; i1 < nTaxa ; i1++) {
            for(int i2 = i1 + 1 ; i2 < nTaxa ; i2++) {
                for(int i3 = i2 + 1 ; i3 < nTaxa ; i3++) {
                    for(int i4 = i3 + 1 ; i4 < nTaxa ; i4++) {
                        for(int i5 = i4 + 1 ; i5 < nTaxa ; i5++) {
                            if (index == number) {
                                selectedTaxa.add(taxa.get(i1));
                                selectedTaxa.add(taxa.get(i2));
                                selectedTaxa.add(taxa.get(i3));
                                selectedTaxa.add(taxa.get(i4));
                                selectedTaxa.add(taxa.get(i5));
                                i1 = nTaxa;
                                i2 = nTaxa;
                                i3 = nTaxa;
                                i4 = nTaxa;
                                break;
                            }
                            index++;
                        }
                    }
                }
            }
        }

        if(selectedTaxa.size() != 5) {
            throw new RuntimeException("Wrong number: " + number);
        }
        System.out.println();
        for(int i = 0 ; i < selectedTaxa.size() ; i++) {
            System.out.print(selectedTaxa.get(i) + "\t");
        }
        System.out.println();

        Map<String, Map<String, String>> subsetData = new HashMap<>();
        for(String locus : multilocusData.keySet()) {
            subsetData.put(locus, new HashMap<>());
            for(String taxon : selectedTaxa) {
                for(String allele : taxonMap.get(taxon)) {
                    subsetData.get(locus).put(allele, multilocusData.get(locus).get(allele));
                }
            }
        }

        Map<String, List<String>> subsetTaxonMap = new HashMap<>();
        for(String taxon : selectedTaxa) {
            subsetTaxonMap.put(taxon, taxonMap.get(taxon));
        }
        Utils._TAXON_MAP = subsetTaxonMap;

        List<Alignment> alignments = new ArrayList<>();
        for(String key : subsetData.keySet()) {
            alignments.add(new Alignment(subsetData.get(key), key));
        }

        Collections.sort(alignments);
        MC3Core mc3 = new MC3Core(alignments);
        mc3.run();
    }

    public static SNSummary stage2(List<String> filenames, int chainlen, int burnin, int sample_freq, SNOptions options) {
        SNProblem problem = new SNProblem();

        class Interceptor extends PrintStream
        {
            public Interceptor(OutputStream out)
            {
                super(out, true);
            }
            @Override
            public void println(String s)
            {//do what ever you like
                //super.print(s);
            }
            @Override
            public void println(char c) {

            }
            @Override
            public void println(int i) {

            }

        }

        PrintStream origOut = System.out;
        PrintStream interceptor = new Interceptor(origOut);

        if(!printDetails_) {
            System.setOut(interceptor);
        }

        int start = (int)(burnin / sample_freq) + 1;
        int end = (int)(chainlen / sample_freq);

        for(String filename : filenames) {

            try {
                BufferedReader in = new BufferedReader(new FileReader(filename));
                String s;
                int index = 0;
                String topSample = null;
                Map<Network, Integer> count = new HashMap<>();
                boolean begin = false;
                int total = 0;
                double lastESS = 0;

                while((s = in.readLine()) != null) {

                    if(begin) {
                        if (s.startsWith("[")) {
                            Network curSample = Networks.readNetworkWithRootPop(s);
                            total++;
                            if(total >= start) {
                                boolean exist = false;
                                for (Network net : count.keySet()) {
                                    if (Networks.hasTheSameTopology(net, curSample)) {
                                        count.put(net, count.get(net) + 1);
                                        exist = true;
                                        break;
                                    }
                                }
                                if (!exist) {
                                    count.put(curSample, 1);
                                }
                            }
                        } else {
                            String ss[] = s.split("\\s+");
                            if(ss.length == 7 && Character.isDigit(ss[2].charAt(0)) && ss[2].charAt(ss[2].length() - 1) == ';')
                                lastESS = Double.parseDouble(ss[2].substring(0, ss[2].length() - 1));
                        }
                    } else {
                        if(s.contains("Logger")) {
                            begin = true;
                        }
                    }
                    index++;
                }

                int totalValue = 0;
                int maxValue = 0;
                for(Network net : count.keySet()) {
                    if(maxValue < count.get(net)) {
                        topSample = Networks.getFullString(net);
                        maxValue = count.get(net);
                    }
                    totalValue += count.get(net);
                }


                in.close();


                if(topSample == null) {
                    continue;
                }

                if(lastESS < 20 || 1.0 * maxValue / totalValue < 0.5) {
                    continue;
                }

                SummaryBL sbl = new SummaryBL(topSample);
                sbl.addFile(filename, true, start, end);
                sbl.report(1.0, 1.0);
                Network meanNet = sbl.getMeanNetwork();
                problem.AddSubNetwork(meanNet, filename, sbl.getPercentage());
                //problem.AddSubNetwork(Networks.readNetworkWithRootPop(topSample), filename, 1.0 * maxValue / totalValue);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        if(printDetails_) {
            System.out.println("Subnetworks: " + problem.GetSize());
            System.out.println(problem.toString());
        }

        if(!printDetails_) {
            System.setOut(origOut);
        }

        return SNSolver.Solve(problem, options) ;
    }

    public static SNSummary stage2_ideal(SNOptions options) {
        Network trueNetwork = options.trueNetwork.clone();
        System.out.println("True network: " + Networks.getDendroscopeCompatibleString(trueNetwork));

        List<Network> netlist = NetworkUtils.genAllSubNetworks(trueNetwork, 3);

        SNProblem problem = new SNProblem();
        Random random = new Random(12345678L);
        for(Network network : netlist) {
            Network newnet = network.clone();
            NetworkUtils.alterHeights(newnet, random);
            problem.AddSubNetwork(newnet, "", 1.0);
        }

        SNSummary summary = SNSolver.Solve(problem, options);

        return summary;
    }

    private static void GetInputFromFolder(List<String> filenames, int start, int end, Map<String, List<String>> file2samples, Map<String, Tuple<String, Double>> file2topsample) {
        for(String filename : filenames) {
            file2samples.put(filename, new ArrayList<>());

            try {
                BufferedReader in = new BufferedReader(new FileReader(filename));
                String s;
                int index = 0;
                String topSample = null;
                Map<Network, Integer> count = new HashMap<>();
                boolean begin = false;
                int total = 0;
                double lastESS = 0.0;

                while((s = in.readLine()) != null) {

                    if(begin) {
                        if (s.startsWith("[")) {
                            Network curSample = Networks.readNetworkWithRootPop(s);
                            total++;
                            if(total >= start) {
                                if(total > end) break; // !!!!!

                                file2samples.get(filename).add(s);
                                boolean exist = false;
                                for (Network net : count.keySet()) {
                                    if (Networks.hasTheSameTopology(net, curSample)) {
                                        count.put(net, count.get(net) + 1);
                                        exist = true;
                                        break;
                                    }
                                }
                                if (!exist) {
                                    count.put(curSample, 1);
                                }
                            }
                        } else {
                            String ss[] = s.split("\\s+");
                            if(ss.length == 7 && Character.isDigit(ss[2].charAt(0)) && ss[2].charAt(ss[2].length() - 1) == ';')
                                lastESS = Double.parseDouble(ss[2].substring(0, ss[2].length() - 1));

                        }
                    } else {
                        if(s.contains("Logger")) {
                            begin = true;
                        }
                    }
                    index++;
                }

                int totalValue = 0;
                int maxValue = 0;
                for(Network net : count.keySet()) {
                    if(maxValue < count.get(net)) {
                        topSample = Networks.getFullString(net);
                        maxValue = count.get(net);
                    }
                    totalValue += count.get(net);
                }


                in.close();


                if(topSample == null) {
                    file2samples.remove(filename);
                    continue;
                }

                if(file2samples.get(filename).size() == 0) {
                    file2samples.remove(filename);
                    continue;
                }
//                if(lastESS < 20 || 1.0 * maxValue / totalValue < 0.5) {
//                    file2samples.remove(filename);
//                    continue;
//                }

//                List<String> toRemove = new ArrayList<>();
//                for(String sample : file2samples.get(filename)) {
//                    if(!Networks.hasTheSameTopology(Networks.readNetworkWithRootPop(sample), Networks.readNetworkWithRootPop(topSample))) {
//                        toRemove.add(sample);
//                    }
//                }
//                file2samples.get(filename).removeAll(toRemove);

                file2topsample.put(filename, new Tuple<>(topSample, 1.0 * maxValue / totalValue));
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    private static void GetInputFromTrueNetwork(List<String> filenames, int num, Network trueNetworkInput, Random random, Map<String, List<String>> file2samples, Map<String, Tuple<String, Double>> file2topsample) {
        Network trueNetwork = trueNetworkInput.clone();

        List<Network> netlist = NetworkUtils.genAllSubNetworks(trueNetwork, 3);
        Map<String, Network> filename2truesubnet = new HashMap<>();
        for(Network network : netlist) {
            Network newnet = network.clone();
            String filename = "";
            for(Object leafObj : newnet.getLeaves()) {
                NetNode leaf = (NetNode) leafObj;
                filename += leaf.getName() + " ";
            }
            filename2truesubnet.put(filename, newnet);
            file2topsample.put(filename, new Tuple<>(newnet.toString(), 1.0));
            file2samples.put(filename, new ArrayList<>());
        }

        for(int i = 0 ; i < num ; i++) {
            for(String filename : filename2truesubnet.keySet()) {
                Network newnet = filename2truesubnet.get(filename).clone();
                NetworkUtils.alterHeights(newnet, random);
                file2samples.get(filename).add(newnet.toString());
            }
        }
    }

    public static SNSummary stage2_1(List<String> filenames, int chainlen, int burnin, int sample_freq, SNOptions options) {

        long starttime = System.currentTimeMillis();

        class Interceptor extends PrintStream
        {
            public Interceptor(OutputStream out)
            {
                super(out, true);
            }
            @Override
            public void println(String s)
            {//do what ever you like
                //super.print(s);
            }
            @Override
            public void println(char c) {

            }
            @Override
            public void println(int i) {

            }

        }

        PrintStream origOut = System.out;
        PrintStream interceptor = new Interceptor(origOut);

        if(!printDetails_) {
            System.setOut(interceptor);
        }

        int start = (int)(burnin / sample_freq) + 1;
        int end = (int)(chainlen / sample_freq);

        Map<String, List<String>> file2samples = new HashMap<>();
        Map<String, Tuple<String, Double>> file2topsample = new HashMap<>();

        if(filenames != null) {
            GetInputFromFolder(filenames, start, end, file2samples, file2topsample);
        } else {
            filenames = new ArrayList<>();
            GetInputFromTrueNetwork(filenames, 100, options.trueNetwork, new Random(12345678L), file2samples, file2topsample);
        }

        SNProblem problem0 = new SNProblem();
        for(String filename : file2topsample.keySet()) {
            Tuple<String, Double> tuple = file2topsample.get(filename);
            SummaryBL summaryBL = new SummaryBL(tuple.Item1);
            for(String sample : file2samples.get(filename)) {
                Network curnet = Networks.readNetworkWithRootPop(sample);
                if(Networks.hasTheSameTopology(Networks.readNetworkWithRootPop(tuple.Item1), curnet)) {
                    summaryBL.addNetwork(curnet);
                }
            }
            summaryBL.computeMean(1,1);
            Network meannet = summaryBL.getMeanNetwork();

            problem0.AddSubNetwork(meannet, filename, tuple.Item2);
        }
        SuperNetwork3 sn0 = new SuperNetwork3(problem0);
        List<String> buildOrder = null;
        Set<String> backboneLeaves = null;
        if(options.tripletFilename != null) {
            sn0.ReduceTrinets(options.allele2species, options.tripletFilename);
            SuperNetwork3.eps = 0.001;
            List<Set<String>> requiredTrinets = sn0.FindMoreRequiredTrinets();
            buildOrder = sn0.GetBuildOrder();
            backboneLeaves = sn0.GetBackboneLeaves();

            while(requiredTrinets.size() > 0) {

                System.out.println("Need more trinets: " + requiredTrinets.size());
                for (Set<String> requiredTaxa : requiredTrinets) {
                    List<String> taxa = new ArrayList<>(requiredTaxa);
                    System.out.println(taxa.get(0) + " " + taxa.get(1) + " " + taxa.get(2));
                }

                sn0.AddBackReduceTrinet(requiredTrinets);
                requiredTrinets = sn0.FindMoreRequiredTrinets();
                buildOrder = sn0.GetBuildOrder();
                backboneLeaves = sn0.GetBackboneLeaves();
            }

            sn0.SetBackboneLeaves(backboneLeaves);
            sn0.SetBuildOrder(buildOrder);

            List<String> reducedFilenames = new ArrayList<>();
            List<String> removedFilenames = new ArrayList<>();
            for (SuperNetwork3.NetworkWithInfo netinfo : sn0.subnetworks_) {
                reducedFilenames.add(netinfo.filename);
            }

            removedFilenames.addAll(filenames);
            removedFilenames.removeAll(reducedFilenames);

            for(String filename : removedFilenames) {
                file2samples.remove(filename);
                file2topsample.remove(filename);
            }
            SuperNetwork3.eps = 0.01;
            options.tripletFilename = null;
            options.buildOrder = buildOrder;
            options.backboneLeaves = backboneLeaves;
        }

//        if(options.tripletFilename != null) {
//            sn0.ReduceTrinets(options.allele2species, options.tripletFilename);
//            sn0.CheckReducedTrinets();
//
//            List<String> reducedFilenames = new ArrayList<>();
//            List<String> removedFilenames = new ArrayList<>();
//            for (SuperNetwork3.NetworkWithInfo netinfo : sn0.subnetworks_) {
//                reducedFilenames.add(netinfo.filename);
//            }
//
//            removedFilenames.addAll(filenames);
//            removedFilenames.removeAll(reducedFilenames);
//
//            for(String filename : removedFilenames) {
//                file2samples.remove(filename);
//                file2topsample.remove(filename);
//            }
//
//            options.tripletFilename = null;
//        }


        int i = 0;
        int totalNumber = 0;
        int goodNumber = 0;
        boolean finished = false;
        SNSummary allsummary = new SNSummary();
        Map<Network, Integer> count = new HashMap<>();
        Map<Network, SummaryBL> summaryBL = new HashMap<>();
        Map<Network, List<Integer>> indexList = new HashMap<>();

        while(!finished && i < 100) {
            SNProblem problem = new SNProblem();
            for(String filename : file2samples.keySet()) {
//                if(i >= file2samples.get(filename).size()) {
//                    finished = true;
//                    System.out.println(filename);
//                    break;
//                }

                //problem.AddSubNetwork(Networks.readNetworkWithRootPop(file2samples.get(filename).get(i)), filename, 1.0);
                problem.AddSubNetwork(Networks.readNetworkWithRootPop(file2samples.get(filename).get(Randomizer.getRandomInt(file2samples.get(filename).size()))), filename, 1.0);
            }
            if(finished) break;

            //if(i < 18) {i++;continue;}

            SNSummary summary = SNSolver.Solve(problem, options);
            System.out.println(i);
            totalNumber++;

            if(summary.inferredNetwork == null) {i++;continue;}

            if(options.trueNetwork != null) {
                if (Networks.hasTheSameTopology(summary.inferredNetwork, options.trueNetwork.clone())) {
                    System.out.println("Good");
                    goodNumber++;
                } else {
                    System.out.println("Not Good");
                }
            }

            boolean exist = false;
            for (Network net : count.keySet()) {
                if (Networks.hasTheSameTopology(net, summary.inferredNetwork)) {
                    count.put(net, count.get(net) + 1);
                    summaryBL.get(net).addNetwork(summary.inferredNetwork);
                    indexList.get(net).add(i);
                    exist = true;
                    break;
                }
            }
            if (!exist) {
                count.put(summary.inferredNetwork, 1);
                summaryBL.put(summary.inferredNetwork, new SummaryBL(Networks.getFullString(summary.inferredNetwork)));
                summaryBL.get(summary.inferredNetwork).addNetwork(summary.inferredNetwork.clone());
                indexList.put(summary.inferredNetwork,  new ArrayList<>());
                indexList.get(summary.inferredNetwork).add(i);
            }

            System.out.println(Networks.getDendroscopeCompatibleString(summary.inferredNetwork));

            //System.out.println(NetworkUtils.ComputeScore(summary.inferredNetwork, problem));
            System.out.println(NetworkUtils.ComputeScore(summary.inferredNetwork, sn0));

            i++;
            if(i == 100) break;

            //break;
            //return summary;

        }

        SuperNetwork3.printDetails_ = false;
        int bestCount = 0;
        double bestScore = 0;
        Network bestMeanNet = null;
        for (Network net : count.keySet()) {
            summaryBL.get(net).computeMean(1,1);
            Network meannet = summaryBL.get(net).getMeanNetwork();
            if(Double.isNaN(meannet.findNode(SuperNetwork3.outgroup).getParentDistance((NetNode) meannet.findNode(SuperNetwork3.outgroup).getParents().iterator().next()))) continue;
            meannet = SuperNetwork3.getSubNetwork(meannet, sn0.getTaxaNames(), true).Item1;
            System.out.println(Networks.getFullString(meannet));
            System.out.println(Networks.getCoalUnitString(meannet));
            System.out.println(Networks.getDendroscopeCompatibleString(meannet));
            System.out.println(count.get(net));
            double score = NetworkUtils.ComputeScore(net, sn0);
            System.out.println(score);
            System.out.print("Indices: ");
            for(Integer index : indexList.get(net)) {
                System.out.print(index + " ");
            }
            System.out.println();
            System.out.println();

            //if(count.get(net) <= 1) continue;
            if(1.0 * count.get(net) / totalNumber > 2.0/3.0) {
                bestCount = count.get(net);
                bestScore = score;
                bestMeanNet = meannet;
                break;
            }

            //if(bestCount < count.get(net) || (bestCount == count.get(net) && bestScore < score)) {
            if(bestScore < score || (bestScore == score && bestCount < count.get(net))) {
                bestCount = count.get(net);
                bestScore = score;
                bestMeanNet = meannet;
            }
        }

        if(!printDetails_) {
            System.setOut(origOut);
        }

        System.out.println("Good number: " + goodNumber);
        System.out.println("Total number: " + totalNumber);
        System.out.println();

        System.out.println("Final number: " + bestCount);
        System.out.println("Final score: " + bestScore);
        System.out.println("Final result: ");
        System.out.println(bestMeanNet);


        SNSummary summary = new SNSummary();
        summary.netinfos = sn0.subnetworks_;
        summary.inferredNetwork = bestMeanNet;
        summary.taxaNames = sn0.getTaxaNames();

        System.out.println("Time (s): " + (System.currentTimeMillis()-starttime)/1000.0);

        return summary ;
    }

    public static void test1(String[] args) {
        Network<Double> net = Networks.readNetwork("[0.01](A:0.20138876485031905:0.01,(B:0.13682431493516344:0.01,((C:0.03614402921753337:0.01,(I:0.005:0.01)#H1:0.03114402921753337:0.01:0.3):0.07:0.01,((((D:0.02:0.01,(L:0.015:0.01)#H3:0.005:0.01:0.2):0.04864873769305778:0.01,((E:0.04784096272944867:0.01,((F:0.01053172961628988:0.01,#H1:0.00553172961628988:0.01:0.7):0.02:0.01,G:0.03053172961628988:0.01):0.017309233113158788:0.01):0.007107765870675464:0.01,H: 0.0549487286:0.01):0.013700009092933646:0.01):0.013126652986989553:0.01,((J:0.03:0.01, (O:0.01:0.01)#H2:0.02:0.01:0.4):0.02689143414599366:0.01,K:0.05689143414599366:0.01):0.024883956534053671:0.01):0.0116928336785998616:0.01,((#H3:0.019882393305892485:0.01:0.8,M:0.044882393305892485:0.01):0.032226283989763562:0.01,(N:0.06296395354420453:0.01,(#H2:0.01361286516308759:0.01:0.6,P:0.02361286516308759:0.01):0.03935108838111694:0.01):0.014144723751451515:0.01):0.016359547062991147:0.01):0.012675804858886178:0.01):0.030680285717630068:0.01):0.06456444991515561:0.01);");
        for(NetNode<Double> node : Networks.postTraversal(net))  {
            for(NetNode<Double> parent : node.getParents()) {
                DecimalFormat df = new DecimalFormat("#.###");
                df.setRoundingMode(RoundingMode.HALF_UP);
                node.setParentDistance(parent, node.getParentDistance(parent) / node.getParentSupport(parent) * 2.0);
                node.setParentSupport(parent, NetNode.NO_SUPPORT);
                node.setParentDistance(parent, Double.parseDouble(df.format(node.getParentDistance(parent))));
            }
        }

        System.out.println(net.toString());
        UltrametricNetwork unet = new UltrametricNetwork(net.toString());
        System.out.println(unet.isUltrametric());
    }

    public static Map<String, String> parseBeastFields(String input) {
        Map<String, String> result = new HashMap<>();
        input += ",";
        StringBuilder key = new StringBuilder();
        StringBuilder value = new StringBuilder();
        boolean readingKey = true;
        boolean inBracket = false;
        for(int i = 0 ; i < input.length() ; i++) {
            if(input.charAt(i) == ',' && !inBracket) {
                result.put(key.toString(), value.toString());
                key = new StringBuilder();
                value = new StringBuilder();
                readingKey = true;
            } else if(input.charAt(i) == '=') {
                readingKey = false;
            } else if(input.charAt(i) == '{') {
                inBracket = true;

            } else if(input.charAt(i) == '}') {
                inBracket = false;
            }
            else {
                if(readingKey) {
                    key.append(input.charAt(i));
                } else {
                    value.append(input.charAt(i));
                }
            }
        }

        return result;
    }

    public static String convertBeastNetworkString(String input) {
        StringBuilder sb = new StringBuilder();
        Map<String, Double> probs = new HashMap<>();
        String netnode = null;

        for(int i = 0 ; i < input.length() ; i++) {
            if(input.charAt(i) == '#') {
                int left = i;
                while(input.charAt(left) != '[') left++;
                netnode = input.substring(i, left);
                sb.append(netnode);
                i = left - 1;
            } else if(input.charAt(i) == '[') {
                int left = i + 2;
                int right = i + 2;
                while(input.charAt(right) != ']') right++;
                Map<String, String> fields = parseBeastFields(input.substring(left, right));
                int colon = right + 1;

                int endOfNumber = colon + 1;
                if(input.charAt(colon) != ';') {
                    while (Character.isDigit(input.charAt(endOfNumber)) || input.charAt(endOfNumber) == '.' || input.charAt(endOfNumber) == 'E' || input.charAt(endOfNumber) == '-')
                        endOfNumber++;
                    fields.put("BranchLength", input.substring(colon + 1, endOfNumber));
                } else {
                    endOfNumber = colon;
                }

                if(fields.containsKey("BranchLength")) {
                    sb.append(":");
                    sb.append(fields.get("BranchLength"));
                }
                if(netnode != null) {
                    if(false) {
                        if (fields.containsKey("gamma")) {
                            sb.append("::");
                            sb.append(fields.get("gamma"));

                            probs.put(netnode, Double.parseDouble(fields.get("gamma")));
                        } else {
                            sb.append("::");
                            sb.append(1.0 - probs.get(netnode));
                        }
                    }
                }

                netnode = null;
                i = endOfNumber - 1;
            } else {
                sb.append(input.charAt(i));
            }
        }

        return sb.toString();
    }

    public static boolean isBackboneOf(Network net1, Network net2) {
        List<Network> backbones = getAllBackboneNets(net2, Integer.MAX_VALUE);
        for(Network backbone : backbones) {
            if(Networks.hasTheSameTopology(net1, backbone)) {
                return true;
            }
        }
        return false;
    }

    static double ComputeNetworkDifficulty(Network net) {
        double score = 0;

        for(Object nodeObj : net.dfs()) {
            NetNode node = (NetNode) nodeObj;
            if(node.isNetworkNode()) {
                int diameter = SuperNetwork3.getReticulationNodeDiameter(node);
                int leavesUnder = SuperNetwork3.getTaxaNamesUnderReticulation(node).size();
                int depend = 0;
                Set<NetNode> ancestors = SuperNetwork3.getAllAncestors(node);
                for(NetNode ancestor : ancestors) {
                    if(ancestor.isNetworkNode()) {
                        depend++;
                    }
                }

                Iterator it = node.getParents().iterator();
                NetNode parent1 = (NetNode) it.next();
                NetNode parent2 = (NetNode) it.next();

                int leavesUnder1 = SuperNetwork3.getTaxaNamesUnderReticulation(parent1).size();
                int leavesUnder2 = SuperNetwork3.getTaxaNamesUnderReticulation(parent2).size();
                score += leavesUnder1 + leavesUnder2 + leavesUnder + depend * net.getLeafCount();
            }
        }

        return score;
    }

    public static void getAllBackbonesDfs(int index, int retiLimit, Network cur, List<NetNode> curReticulations, List<Network> results) {
        if(index >= curReticulations.size()) return;

        NetNode curnode = curReticulations.get(index);
        List<NetNode> parents = new ArrayList<>();
        for(Object parentObj : curnode.getParents()) {
            NetNode parent = (NetNode) parentObj;
            parents.add(parent);
        }

        for(NetNode parent : parents) {
            double distanceBackup = curnode.getParentDistance(parent);
            parent.removeChild(curnode);
            if(cur.getReticulationCount() <= retiLimit) {
                Network temp = cur.clone();
                Networks.removeBinaryNodes(temp);
                results.add(temp);
            }
            getAllBackbonesDfs(index + 1, retiLimit, cur, curReticulations, results);

            parent.adoptChild(curnode, distanceBackup);
        }

        getAllBackbonesDfs(index + 1, retiLimit, cur, curReticulations, results);
    }

    public static List<Network> getBackboneNetsWithNumReti(Network network, int numReti) {
        Network clonedNetwork = network.clone();

        numReti = clonedNetwork.getReticulationCount() - numReti;
        List<NetNode> trueReticulations = new ArrayList<>();
        for(Object nodeObj : clonedNetwork.dfs()) {
            NetNode node = (NetNode) nodeObj;
            if(node.isNetworkNode()) {
                trueReticulations.add(node);
                for(Object parentObj : node.getParents()) {
                    NetNode parent = (NetNode) parentObj;
                    node.setParentProbability(parent, NetNode.NO_PROBABILITY);
                }
            }
        }

        List<Network> results = new ArrayList<>();

        List<int[]> subsets = new ArrayList<>();

        int[] s = new int[numReti];

        if (numReti <= trueReticulations.size()) {
            // first index sequence: 0, 1, 2, ...
            for (int i = 0; (s[i] = i) < numReti - 1; i++);
            subsets.add(s.clone());
            for(;;) {
                int i;
                // find position of item that can be incremented
                for (i = numReti - 1; i >= 0 && s[i] == trueReticulations.size() - numReti + i; i--);
                if (i < 0) {
                    break;
                }
                s[i]++;                    // increment this item
                for (++i; i < numReti; i++) {    // fill up remaining items
                    s[i] = s[i - 1] + 1;
                }
                subsets.add(s.clone());
            }
        }

        for(int[] indices : subsets) {
            getBackbonesWithNumRetiDfs(0, clonedNetwork, trueReticulations, indices, results);
        }

        return results;
    }

    public static void getBackbonesWithNumRetiDfs(int index, Network cur, List<NetNode> curReticulations, int[] indices, List<Network> results) {
        if(index == indices.length) {
            Network temp = cur.clone();
            Networks.removeBinaryNodes(temp);
            results.add(temp);
        }

        if(index >= indices.length) return;

        NetNode curnode = curReticulations.get(indices[index]);
        List<NetNode> parents = new ArrayList<>();
        for(Object parentObj : curnode.getParents()) {
            NetNode parent = (NetNode) parentObj;
            parents.add(parent);
        }

        if(parents.size() < 2) return;

        for(NetNode parent : parents) {
            double distanceBackup = curnode.getParentDistance(parent);
            parent.removeChild(curnode);

            getBackbonesWithNumRetiDfs(index + 1, cur, curReticulations, indices, results);

            parent.adoptChild(curnode, distanceBackup);
        }

    }

    public static List<Network> getAllBackboneNets(Network network, int retiLimit) {
        Network clonedNetwork = network.clone();

        List<NetNode> trueReticulations = new ArrayList<>();
        for(Object nodeObj : clonedNetwork.dfs()) {
            NetNode node = (NetNode) nodeObj;
            if(node.isNetworkNode()) {
                trueReticulations.add(node);
                for(Object parentObj : node.getParents()) {
                    NetNode parent = (NetNode) parentObj;
                    node.setParentProbability(parent, NetNode.NO_PROBABILITY);
                }
            }
        }

        List<Network> results = new ArrayList<>();

        if(network.getReticulationCount() > retiLimit) {

        }

        getAllBackbonesDfs(0, retiLimit, clonedNetwork, trueReticulations, results);
        return results;
    }

    public static Tuple3<Network, Network, Double> CheckWithTrueBackbone(Network inferredNetwork, Network trueNetwork) {
        List<Network> trueBackboneNets = getAllBackboneNets(trueNetwork, Integer.MAX_VALUE);
        trueBackboneNets.add(0, trueNetwork.clone());

        List<Network> inferredBackboneNets = new ArrayList<>();
        inferredBackboneNets.add(0, inferredNetwork);

        Tuple3<Network, Network, Double> closest = null;
        for(Network inferredBackboneNet : inferredBackboneNets) {
            for (Network trueBackboneNet : trueBackboneNets) {
                double dist = Networks.computeDistanceBetweenTwoNetworks(trueBackboneNet, inferredBackboneNet);
                if (closest == null || closest.Item3 > dist || (closest.Item3 == dist && inferredBackboneNet.getReticulationCount() > closest.Item1.getReticulationCount())) {
                    closest = new Tuple3<>(inferredBackboneNet, trueBackboneNet, dist);
                }
            }
        }
        return closest;
    }

    public static Tuple3<Network, Network, Double> CheckWithTrueNetwork(Network inferredNetwork, Network trueNetwork) {
        List<String> leaves = new ArrayList<>();
        for(Object leafObj : trueNetwork.getLeaves()) {
            NetNode leaf = (NetNode) leafObj;
            leaves.add(leaf.getName());
        }

        List<Network> trueBackboneNets = getAllBackboneNets(trueNetwork, Integer.MAX_VALUE);
        trueBackboneNets.add(0, trueNetwork.clone());

        List<Network> inferredBackboneNets = getAllBackboneNets(inferredNetwork, Integer.MAX_VALUE);
        inferredBackboneNets.add(0, inferredNetwork);

        Tuple3<Network, Network, Double> closest = null;
        for(Network inferredBackboneNet : inferredBackboneNets) {
            if(inferredBackboneNet.getLeafCount() != trueNetwork.getLeafCount())  {
                inferredBackboneNet = SuperNetwork3.getSubNetwork(inferredBackboneNet, leaves, true).Item1;
            }

            for (Network trueBackboneNet : trueBackboneNets) {
                if(inferredBackboneNet.getReticulationCount() != trueBackboneNet.getReticulationCount()) continue;
                double dist = Networks.computeDistanceBetweenTwoNetworks(trueBackboneNet, inferredBackboneNet);
                if (closest == null || closest.Item3 > dist || (closest.Item3 == dist && inferredBackboneNet.getReticulationCount() > closest.Item1.getReticulationCount())) {
                    closest = new Tuple3<>(inferredBackboneNet, trueBackboneNet, dist);
                }
            }
        }
        return closest;
    }

    public static int CompareNodes(Network inferredNetwork, Network trueNetwork) {

        boolean finished = false;
        List<String> currentSubset = new ArrayList<>();

        for(Object leafObj : trueNetwork.getLeaves()) {
            NetNode leaf = (NetNode) leafObj;
            currentSubset.add(leaf.getName());
        }


        Network curINet = SuperNetwork3.getSubNetwork(inferredNetwork, currentSubset, true).Item1;
        Network curTNet = trueNetwork.clone();

        if(Networks.hasTheSameTopology(curINet, curTNet)) {
            return curINet.getLeafCount() * 2 - 1 + curINet.getReticulationCount() * 2;
        }


        while(!finished) {
            Network nextINet = null;
            Network nextTNet = null;
            int nextNodeCountI = 0;
            String leaf2Delete = null;
            double nextD = Double.MAX_VALUE;
            for(String targetName : currentSubset) {
                List<String> selectedLeaves = new ArrayList<>(currentSubset);
                selectedLeaves.remove(targetName);
                Tuple<Network, Map<NetNode, NetNode>> tupleI = SuperNetwork3.getSubNetwork(curINet, selectedLeaves, true);
                Tuple<Network, Map<NetNode, NetNode>> tupleT = SuperNetwork3.getSubNetwork(curTNet, selectedLeaves, true);
                int nodeCountI = tupleI.Item1.getLeafCount() * 2 - 1 + tupleI.Item1.getReticulationCount();
                int nodeCountT = tupleT.Item1.getLeafCount() * 2 - 1 + tupleT.Item1.getReticulationCount();

                double dist =  Networks.computeDistanceBetweenTwoNetworks(tupleI.Item1, tupleT.Item1);

                if(nextD > dist || (nextD == dist && (nodeCountI > nextNodeCountI))) {
                    nextD = dist;
                    nextINet = tupleI.Item1;
                    nextTNet = tupleT.Item1;
                    nextNodeCountI = nodeCountI;
                    leaf2Delete = targetName;
                }
            }

            System.out.println(leaf2Delete);
            currentSubset.remove(leaf2Delete);
            curINet = nextINet;
            curTNet = nextTNet;

            if(Networks.hasTheSameTopology(curINet, curTNet)) {
                System.out.println("Inferred network after deletion: ");
                System.out.println(curINet);
                System.out.println("True network after deletion: ");
                System.out.println(curTNet);

                finished = true;
                return nextNodeCountI;
            }

        }

        return 0;
    }

    static public List<Set<String>> convertAlleleTripletsToSpeciesTriplets(Map<String, List<String>> species2alleles, String outgroup, Set<String[]> alleleTriplets) {
        List<Set<String>> result = new ArrayList<>();
        if(species2alleles == null) {
            for(String[] triplet : alleleTriplets) {
                Set<String> cur_alleles = new HashSet<>();
                cur_alleles.add(triplet[0]);
                cur_alleles.add(triplet[1]);
                cur_alleles.add(triplet[2]);
                result.add(cur_alleles);
            }

            return result;
        }

        Map<String, String> allele2species = new HashMap<>();
        for(String species : species2alleles.keySet()) {
            for(String allele : species2alleles.get(species)) {
                allele2species.put(allele, species);
            }
        }

        for(String[] triplet : alleleTriplets) {
            List<String> cur_alleles = new ArrayList<>();
            cur_alleles.add(triplet[0]);
            cur_alleles.add(triplet[1]);
            cur_alleles.add(triplet[2]);

            Set<String> cur_species = new HashSet<>();

            for(String allele : cur_alleles) {
                cur_species.add(allele2species.get(allele));
            }

            if(cur_species.size() == 1) {
                continue;
            } else if(cur_species.size() == 2) {
                cur_species.add(outgroup);
            }

            if(!result.contains(cur_species)) {
                result.add(cur_species);
            }
        }

        return result;
    }

    public static SNSummary mergeBayesianResults(List<String> filenames, int chainlen, int burnin, int sample_freq, SNOptions options) {

        long starttime = System.currentTimeMillis();

        class Interceptor extends PrintStream
        {
            public Interceptor(OutputStream out)
            {
                super(out, true);
            }
            @Override
            public void println(String s)
            {//do what ever you like
                //super.print(s);
            }
            @Override
            public void println(char c) {

            }
            @Override
            public void println(int i) {

            }

        }

        PrintStream origOut = System.out;
        PrintStream interceptor = new Interceptor(origOut);

        if(!printDetails_) {
            System.setOut(interceptor);
        }

        int start = (int)(burnin / sample_freq) + 1;
        int end = (int)(chainlen / sample_freq);

        Map<String, List<String>> file2samples = new HashMap<>();
        Map<String, Tuple<String, Double>> file2topsample = new HashMap<>();

        GetInputFromFolder(filenames, start, end, file2samples, file2topsample);

        SNProblem problem0 = new SNProblem();
        for(String filename : file2topsample.keySet()) {
            Tuple<String, Double> tuple = file2topsample.get(filename);
            SummaryBL summaryBL = new SummaryBL(tuple.Item1);
            for(String sample : file2samples.get(filename)) {
                Network curnet = Networks.readNetworkWithRootPop(sample);
                if(Networks.hasTheSameTopology(Networks.readNetworkWithRootPop(tuple.Item1), curnet)) {
                    summaryBL.addNetwork(curnet);
                }
            }
            summaryBL.computeMean(1,1);
            Network meannet = summaryBL.getMeanNetwork();

            problem0.AddSubNetwork(meannet, filename, tuple.Item2);
        }
        SuperNetwork3 sn0 = new SuperNetwork3(problem0);
        List<String> buildOrder = null;
        Set<String> backboneLeaves = null;


        // Check whether trinets are enough
        SuperNetwork3.eps = options.eps / 10.;
        List<Set<String>> requiredTrinets = sn0.FindMoreRequiredTrinets();
        buildOrder = sn0.GetBuildOrder();
        backboneLeaves = sn0.GetBackboneLeaves();

        if(requiredTrinets.size() > 0) {

            System.out.println("Need more trinets: " + requiredTrinets.size());
            for (Set<String> requiredTaxa : requiredTrinets) {
                List<String> taxa = new ArrayList<>(requiredTaxa);
                System.out.println(taxa.get(0) + " " + taxa.get(1) + " " + taxa.get(2));
            }

            System.out.println("Order:");
            for(String s : buildOrder) {
                System.out.print(s + " ");
            }

            System.out.println("Backbone:");
            for(String s : backboneLeaves) {
                System.out.print(s + " ");
            }
            System.out.println();

            return null;
        }

        SuperNetwork3.eps = options.eps;
        options.tripletFilename = null;
        options.buildOrder = buildOrder;
        options.backboneLeaves = backboneLeaves;



        int i = 0;
        int totalNumber = 0;
        int goodNumber = 0;
        boolean finished = false;
        SNSummary allsummary = new SNSummary();
        Map<Network, Integer> count = new HashMap<>();
        Map<Network, SummaryBL> summaryBL = new HashMap<>();
        Map<Network, List<Integer>> indexList = new HashMap<>();

        while(!finished && i < 100) {
            SNProblem problem = new SNProblem();
            for(String filename : file2samples.keySet()) {
//                if(i >= file2samples.get(filename).size()) {
//                    finished = true;
//                    System.out.println(filename);
//                    break;
//                }

                //problem.AddSubNetwork(Networks.readNetworkWithRootPop(file2samples.get(filename).get(i)), filename, 1.0);
                problem.AddSubNetwork(Networks.readNetworkWithRootPop(file2samples.get(filename).get(Randomizer.getRandomInt(file2samples.get(filename).size()))), filename, 1.0);
            }
            if(finished) break;

            //if(i < 18) {i++;continue;}

            SNSummary summary = SNSolver.Solve(problem, options);
            System.out.println(i);
            totalNumber++;

            if(summary.inferredNetwork == null) {i++;continue;}

            if(options.trueNetwork != null) {
                if (Networks.hasTheSameTopology(summary.inferredNetwork, options.trueNetwork.clone())) {
                    System.out.println("Good");
                    goodNumber++;
                } else {
                    System.out.println("Not Good");
                }
            }

            boolean exist = false;
            for (Network net : count.keySet()) {
                if (Networks.hasTheSameTopology(net, summary.inferredNetwork)) {
                    count.put(net, count.get(net) + 1);
                    summaryBL.get(net).addNetwork(summary.inferredNetwork);
                    indexList.get(net).add(i);
                    exist = true;
                    break;
                }
            }
            if (!exist) {
                count.put(summary.inferredNetwork, 1);
                summaryBL.put(summary.inferredNetwork, new SummaryBL(Networks.getFullString(summary.inferredNetwork)));
                summaryBL.get(summary.inferredNetwork).addNetwork(summary.inferredNetwork.clone());
                indexList.put(summary.inferredNetwork,  new ArrayList<>());
                indexList.get(summary.inferredNetwork).add(i);
            }

            System.out.println(Networks.getDendroscopeCompatibleString(summary.inferredNetwork));

            //System.out.println(NetworkUtils.ComputeScore(summary.inferredNetwork, problem));
            System.out.println(NetworkUtils.ComputeScore(summary.inferredNetwork, sn0));

            i++;
            if(i == 100) break;

            //break;
            //return summary;

        }

        SuperNetwork3.printDetails_ = false;
        int bestCount = 0;
        double bestScore = 0;
        Network bestMeanNet = null;
        for (Network net : count.keySet()) {
            summaryBL.get(net).computeMean(1,1);
            Network meannet = summaryBL.get(net).getMeanNetwork();
            if(Double.isNaN(meannet.findNode(SuperNetwork3.outgroup).getParentDistance((NetNode) meannet.findNode(SuperNetwork3.outgroup).getParents().iterator().next()))) continue;
            meannet = SuperNetwork3.getSubNetwork(meannet, sn0.getTaxaNames(), true).Item1;
            System.out.println(Networks.getFullString(meannet));
            System.out.println(Networks.getCoalUnitString(meannet));
            System.out.println(Networks.getDendroscopeCompatibleString(meannet));
            System.out.println(count.get(net));
            double score = NetworkUtils.ComputeScore(net, sn0);
            System.out.println(score);
            System.out.print("Indices: ");
            for(Integer index : indexList.get(net)) {
                System.out.print(index + " ");
            }
            System.out.println();
            System.out.println();

            //if(count.get(net) <= 1) continue;
            if(1.0 * count.get(net) / totalNumber > 2.0/3.0) {
                bestCount = count.get(net);
                bestScore = score;
                bestMeanNet = meannet;
                break;
            }

            //if(bestCount < count.get(net) || (bestCount == count.get(net) && bestScore < score)) {
            if(bestScore < score || (bestScore == score && bestCount < count.get(net))) {
                bestCount = count.get(net);
                bestScore = score;
                bestMeanNet = meannet;
            }
        }

        if(!printDetails_) {
            System.setOut(origOut);
        }

        System.out.println("Good number: " + goodNumber);
        System.out.println("Total number: " + totalNumber);
        System.out.println();

        System.out.println("Final number: " + bestCount);
        System.out.println("Final score: " + bestScore);
        System.out.println("Final result: ");
        System.out.println(bestMeanNet);


        SNSummary summary = new SNSummary();
        summary.netinfos = sn0.subnetworks_;
        summary.inferredNetwork = bestMeanNet;
        summary.taxaNames = sn0.getTaxaNames();

        System.out.println("Time (s): " + (System.currentTimeMillis()-starttime)/1000.0);

        return summary ;
    }

    static public void main(String []args) {
        //Network trueNetwork = Networks.readNetwork("(Z:100.0,(((((D:7.430083706879999,((N:1.2)#H2:2.3831808)#H1:3.846902906879999)S20:3.2692368310271975,(I:2.9859839999999997,H:2.9859839999999997)S19:7.713336537907196)S18:4.707701036679165,#H1:11.823840774586362)S17:22.93057834988837,((E:6.191736422399999)#H4:25.756263514662276)#H3:6.389599987412456)S16:41.159247278916055,((A:55.206143891243606,(((F:5.159780351999999,(J:2.48832,((M:1.44)#H5:0.6335999999999999,K:2.0736)S15:0.41472)S14:2.6714603519999995)S13:7.679404293488635,(#H2:3.099816959999999,G:4.299816959999999)S12:8.539367685488635)S11:33.16593526388104,((#H4:12.296689467103633,(B:8.916100448255998,C:8.916100448255998)S10:9.572325441247633)S9:8.1349073913816,((#H5:0.28800000000000003,L:1.728)S8:20.458111067404356,(O:1.0,P:1.0)S7:21.186111067404358)S6:4.437222213480872)S5:19.381786628484445)S4:9.20102398187393)S3:11.041228778248716,#H3:34.299372732430044)S2:13.249474533898464)S1:20.503152796609214);");
        //Network inferredNetwork = Networks.readNetwork("(Z:0.7096995771149552,(((((H:0.029034319213268003,I:0.029034319213268003)I12:0.07571284181423195,(D:0.07054975280170007,((N:0.01677046154182873)I13#H2:0.03547962985339541::0.29734002016051975)I8#H1:0.018299661406475938::0.45197743073632335)I11:0.03419740822579988)I7:0.04215771020355635,I8#H1:0.09465477983583218::0.5480225692636767)I4:0.21268141868533277,(E:0.14551081013277983)I5#H3:0.21407547978360925::0.2643208014344447)I2:0.23803919763961257,(((((P:0.009021966743936025,O:0.009021966743936025)I20:0.194039105912683,(L:0.015768583495269863,(M:0.015668583495269863)I22#H4:9.99999999999994E-5::0.7841851661421151)I19:0.18729248916134916)I15:0.04056689988800191,((B:0.08646536587599503,C:0.08646536587599503)I18:0.0852119469162843,I5#H3:0.02616650265949949::0.7356791985655553)I14:0.0719506597523416)I9:0.14730339578831697,((G:0.039965889280701807,I13#H2:0.023195427738873075::0.7026599798394803)I16:0.0830744882850622,(F:0.050398625425341764,(J:0.024499579087946654,(I22#H4:0.003481035489195875::0.2158148338578849,K:0.01914961898446574)I23:0.005349960103480916)I21:0.02589904633739511)I17:0.07264175214042223)I10:0.26789099076717393)I6:0.06327200779650383,A:0.4542033761294417)I3:0.14342211142655992)I1:0.1120740895589536)I0;");

        Network trueNetwork = Networks.readNetwork("(((((I:2.0736,H:2.0736)S9:0.9123839999999999,F:2.9859839999999997)S8:5.9301164482559985,(((L:1.44,K:1.44)S7:0.28800000000000003,J:1.728)S6:0.7603199999999999,G:2.48832)S5:6.427780448255998)S4:9.572325441247633,(((M:1.2,N:1.2)S18:4.991736422399999,A:6.191736422399999)S17:9.215285152186361,(((E:3.5831807999999996,(O:1.0,P:1.0)S14:2.5831807999999996)S13:0.7166361599999997,D:4.299816959999999)S12:8.539367685488635,(C:5.159780351999999,B:5.159780351999999)S15:7.679404293488635)S11:2.5678369290977265)S10:3.08140431491727)S3:81.51157411049637,Z:100.0);");
        Network inferredNetwork = Networks.readNetwork("(Z:0.7098935326187776,(((((P:0.00934308447367168,O:0.00934308447367168)I13:0.02678146831338358,E:0.03612455278705526)I8:0.005426082481262834,D:0.04155063526831809)I4:0.0808717831323044,(C:0.049349787332531966,B:0.049349787332531966)I5:0.07307263106809053)I2:0.043025824302415736,(((N:0.012697205524323789,M:0.012697205524323789)I14:0.04801980392309516,A:0.06071700944741895)I9:0.1046312332556193,((((L:0.013593236724210435,K:0.013593236724210435)I17:0.004007457531740718,J:0.017600694255951153)I16:0.007174349510463831,G:0.024775043766414984)I12:0.06306047173805217,(F:0.02934114268517848,(I:0.019812564232474315,H:0.019812564232474315)I15:0.009528578452704165)I11:0.05849437281928867)I7:0.0775127271985711)I3:9.999999999998899E-5)I1:0.5444452899157394)I0;");

        Tuple3<Network, Network, Double> closest = Pipeline.CheckWithTrueNetwork(inferredNetwork, trueNetwork);
        System.out.println("Closest true # Reti: " + closest.Item1.getReticulationCount());
        System.out.println(closest.Item2);
        System.out.println("Closest inferred # Reti: " + closest.Item2.getReticulationCount());
        System.out.println(closest.Item1);
        System.out.println("Distance: " + closest.Item3);

        System.out.println(Networks.hasTheSameTopology(inferredNetwork, trueNetwork));
        System.out.println(CompareNodes(inferredNetwork, trueNetwork));
    }

}
