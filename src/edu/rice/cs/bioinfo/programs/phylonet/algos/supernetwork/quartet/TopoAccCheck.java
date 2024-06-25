package edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork.quartet;
/*
 * @ClassName:   TopoAccCheck
 * @Description:
 * @Author:      Zhen Cao
 * @Date:        3/30/23 4:15 PM
 */

import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.Tuple3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCsnapp.summary.SummaryBL;
import edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork.*;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.*;
import java.util.*;

import static edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork.Pipeline.getAllBackboneNets;
import static edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork.Pipeline.test;

public class TopoAccCheck {
    /* Constructor */
    public TopoAccCheck() {

    }

//    public static Map<Set<String>, Network> inferred_quartets(String quartetdir){
//        int chainlen = 80000000;
//        int burnin = 8000000;
//        int sf = 5000;
//
//        File path = new File(quartetdir);
//        List<String> filenames = new ArrayList<>();
//        File[] files = path.listFiles();
//        for (int i = 0; i < files.length; i++) {
//            if (files[i].isFile()) { //this line weeds out other directories/folders
//                if (files[i].toString().endsWith(".out")) {
//                    //System.out.println(files[i]);
//                    filenames.add(files[i].toString());
//                }
//            }
//        }
//        Map<Set<String>, Network> networkmap = new HashMap();
//
//        Collections.sort(filenames);
//        int start = (int) (burnin*1.0/sf);
//        int end = (int) (chainlen*1.0/sf);
//
//        for(String filename : filenames) {
//
//            try {
//                BufferedReader in = new BufferedReader(new FileReader(filename));
//                String s;
//                int index = 0;
//                String topSample = null;
//                Map<Network, Integer> count = new HashMap<>();
//                boolean begin = false;
//                int total = 0;
//                double lastESS = 0;
//
//                while((s = in.readLine()) != null) {
//
//                    if(begin) {
//                        if (s.startsWith("[")) {
//                            Network curSample = Networks.readNetworkWithRootPop(s);
//                            total++;
//                            if(total >= start) {
//                                boolean exist = false;
//                                for (Network net : count.keySet()) {
//                                    if (Networks.hasTheSameTopology(net, curSample)) {
//                                        count.put(net, count.get(net) + 1);
//                                        exist = true;
//                                        break;
//                                    }
//                                }
//                                if (!exist) {
//                                    count.put(curSample, 1);
//                                }
//                            }
//                        } else {
//                            String ss[] = s.split("\\s+");
//                            if(ss.length == 7 && Character.isDigit(ss[2].charAt(0)) && ss[2].charAt(ss[2].length() - 1) == ';')
//                                lastESS = Double.parseDouble(ss[2].substring(0, ss[2].length() - 1));
//                        }
//                    } else {
//                        if(s.contains("Logger")) {
//                            begin = true;
//                        }
//                    }
//                    index++;
//                }
//
//                int totalValue = 0;
//                int maxValue = 0;
//                for(Network net : count.keySet()) {
//                    if(maxValue < count.get(net)) {
//                        topSample = Networks.getFullString(net);
//                        maxValue = count.get(net);
//                    }
//                    totalValue += count.get(net);
//                }
//
//
//                in.close();
//
//
//                if(topSample == null) {
//                    continue;
//                }
//
////                if(lastESS < 20 || 1.0 * maxValue / totalValue < 0.5) {
////                    continue;
////                }
//
////                SummaryBL sbl = new SummaryBL(topSample);
////                sbl.addFile(filename, true, start, end);
////                sbl.report(1.0, 1.0);
////                Network meanNet = sbl.getMeanNetwork();
//                List<String> leafset = new ArrayList<>();
//                Network net = Networks.readNetwork(topSample);
//                for(Object leafObj : net.getLeaves()) {
//                    NetNode leaf = (NetNode) leafObj;
//                    leafset.add(leaf.getName());
//                }
//                System.out.println(topSample);
//                networkmap.put(new HashSet<>(leafset), net);
//            } catch (IOException e) {
//                e.printStackTrace();
//            }
//        }
//        return networkmap;
//    }



    public static void summarize_inferred_quartets(String quartetdir, String output_path){
        int chainlen = 80000000;
        int burnin = 10000000;
        int sf = 5000;
        try {
            File path = new File(quartetdir);
            List<String> filenames = new ArrayList<>();
            File[] files = path.listFiles();
            BufferedWriter writer = new BufferedWriter(new FileWriter(output_path));


            for (int i = 0; i < files.length; i++) {
                if (files[i].isFile()) { //this line weeds out other directories/folders
                    if (files[i].toString().endsWith(".out")) {
                        //System.out.println(files[i]);
                        filenames.add(files[i].toString());
                    }
                }
            }
            Map<Set<String>, Network> networkmap = new HashMap();

            Collections.sort(filenames);
            int start = (int) (burnin*1.0/sf);
            int end = (int) (chainlen*1.0/sf);

            for(String filename : filenames) {


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

//                if(lastESS < 20 || 1.0 * maxValue / totalValue < 0.5) {
//                    continue;
//                }

//                SummaryBL sbl = new SummaryBL(topSample);
//                sbl.addFile(filename, true, start, end);
//                sbl.report(1.0, 1.0);
//                Network meanNet = sbl.getMeanNetwork();
                List<String> leafset = new ArrayList<>();
//                Network net = Networks.readNetwork(topSample);
//                for(Object leafObj : net.getLeaves()) {
//                    NetNode leaf = (NetNode) leafObj;
//                    leafset.add(leaf.getName());
//                }
                System.out.println(topSample);
                writer.write(topSample);

            }

        }catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static Map<Set<String>, Network> inferred_quartets(String filename){
        Map<Set<String>, Network> res = new HashMap<>();
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
                Network net = Networks.readNetworkWithRootPop(s);
                List<String> leafset = new ArrayList<>();
                for(Object leafObj : net.getLeaves()) {
                    NetNode leaf = (NetNode) leafObj;
                    leafset.add(leaf.getName());
                }
                res.put(new HashSet<>(leafset), net);
            }

        }catch (IOException e){

        }
        return res;
    }



    public static List<List<String>> generateCombinations(List<String> list, int k) {
        List<List<String>> result = new ArrayList<>();
        if (k == 0) {
            result.add(new ArrayList<>());
            return result;
        }
        if (list.isEmpty()) {
            return result;
        }
        String firstElement = list.get(0);
        List<List<String>> combinationsWithoutFirst = generateCombinations(list.subList(1, list.size()), k);
        List<List<String>> combinationsWithFirst = generateCombinations(list.subList(1, list.size()), k - 1);
        for (List<String> combination : combinationsWithFirst) {
            combination.add(0, firstElement);
        }
        result.addAll(combinationsWithoutFirst);
        result.addAll(combinationsWithFirst);
        return result;
    }



    public static Map<Set<String>, Network> true_quartets(String truestring){
        Network truenet = Networks.readNetwork(truestring);
        List<String> taxalist = new ArrayList<>();
        // get taxa list
        for (Object o: truenet.getLeaves()){
            NetNode n = (NetNode) o;
            taxalist.add(n.getName());
        }
        List<Network> network_list = new ArrayList<>();
        Map<Set<String>, Network> hashtable = new HashMap<>();

        // enumrate all 4 combinations
        List<List<String>> result = generateCombinations(taxalist, 4);
        for (List<String> selected :result){
            Tuple<Network, Map<NetNode, NetNode>> tuple = SuperNetwork3.getSubNetwork(truenet, selected, true);
            network_list.add(tuple.Item1);
            hashtable.put(new HashSet<>(selected), tuple.Item1);
        }
        return hashtable;

    }

    public static Tuple3<Network, Network, Double> CheckInsideTrueBackbone(Network inferredNetwork, Network trueNetwork) {
        List<Network> trueBackboneNets = getAllBackboneNets(trueNetwork, Integer.MAX_VALUE);
        trueBackboneNets.add(0, trueNetwork.clone());



        Tuple3<Network, Network, Double> closest = null;
        for (Network trueBackboneNet : trueBackboneNets) {
            double dist = Networks.computeDistanceBetweenTwoNetworks(trueBackboneNet, inferredNetwork);
            if (closest == null || closest.Item3 > dist || (closest.Item3 == dist && inferredNetwork.getReticulationCount() > closest.Item1.getReticulationCount())) {
                closest = new Tuple3<>(inferredNetwork, trueBackboneNet, dist);
            }
        }

        return closest;
    }

    public static void compare_true_infer(){
        String truestring = "(Z:100.0,((((F:4.299816959999999)#H2:14.188608929503632::0.65)#H1:19.8491740349711::0.4,(((B:12.839184645488634,((((N:1.2,O:1.2)S19:0.8735999999999999,K:2.0736)S18:3.0861803519999995,#H2:0.859963392::0.35)S17:1.0319560703999997,E:6.191736422399999)S16:6.647448223088635)S15:2.5678369290977265,(H:2.9859839999999997,(J:2.48832,I:2.48832)S14:0.4976639999999999)S13:12.42103757458636)S12:16.54097836247591,#H1:13.459574047558643::0.6)S11:6.389599987412456)S10:16.868543966768875,((A:26.62333328088523,((((M:1.44)#H3:0.28800000000000003::0.29999999999999993,L:1.728)S9:8.971320537907197,(#H3:5.990083706879998::0.7000000000000001,(P:1.0)#H4:6.430083706879999::0.2499999999999999)S8:3.2692368310271975)S7:11.486790529497162,(D:8.916100448255998,C:8.916100448255998)S6:13.27001061914836)S5:4.437222213480872)S4:19.381786628484445,(G:3.5831807999999996,#H4:2.5831807999999996::0.7500000000000001)S3:42.421939109369674)S2:9.20102398187393)S1:44.793856108756394);";
        String quartet_dir = "/Users/zhen/Desktop/Zhen/research/phylogenetics/mcmc_speedup/data/quartet/Reti4_B_startiq/";
        String quartet_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/mcmc_speedup/data/quartet/inferred_quartets_startiq.txt";

//        summarize_inferred_quartets(quartet_dir, quartet_path);
        Map<Set<String>, Network> true_q = true_quartets(truestring);
        Map<Set<String>, Network> infer_q = inferred_quartets(quartet_path);

        int identical_cnt = 0;
        int inside_cnt = 0;
        int uninferred = 0;

        for (Set<String> leafset : true_q.keySet()){
            Network inferred_net = infer_q.get(leafset);
            Network true_net = true_q.get(leafset);
            List<String> testset = new ArrayList<>();
            testset.add("A");
            testset.add("B");
            testset.add("C");
            testset.add("F");
            if (leafset.containsAll(testset)){
                System.out.println(inferred_net.toString());
                System.out.println(true_net.toString());
            }

            if (inferred_net == null){
                uninferred += 1;
                continue;
            }
            if (Networks.hasTheSameTopology(inferred_net, true_net)){
                identical_cnt += 1;
            }
            else{
//                Tuple3<Network, Network, Double> tuple3 = Pipeline.CheckWithTrueNetwork(inferred_net, true_net);
                Tuple3<Network, Network, Double> tuple3 = CheckInsideTrueBackbone(inferred_net, true_net);
                if (tuple3.Item3 < 0.1){
                    inside_cnt += 1;
                }
            }


//            boolean inside_cnt = display
        }
        System.out.println(identical_cnt);
        System.out.println(inside_cnt);
        System.out.println(2380-uninferred-identical_cnt-inside_cnt);
        System.out.println(uninferred);
    }


    public static void check_stage() {
        String resultFolder = "/Users/zhen/Desktop/Zhen/research/phylogenetics/mcmc_speedup/data/quartet/Reti4_B_slurms/";
        File path = new File(resultFolder);

        List<String> filenames = new ArrayList<>();

        File [] files = path.listFiles();
        for (int i = 0; i < files.length; i++){
            if (files[i].isFile()){ //this line weeds out other directories/folders
                if(files[i].toString().endsWith(".out")) {
                    //System.out.println(files[i]);
                    filenames.add(files[i].toString());
                }
            }
        }

        Collections.sort(filenames);
        SNOptions options = new SNOptions();
        options.eps = 0.01;

        SuperNetwork3.printDetails_ = true;
        SNSummary summary = Pipeline.stage2_1(filenames, 80000000, 8000000, 5000, options);

        Network inferred = summary.inferredNetwork;
        System.out.println(inferred);
    }

    public static void check_compute(){
        String resultPath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/mcmc_speedup/data/quartet/inferred_quartets.txt";
        String trueString = "(Z:100.0,((((F:4.299816959999999)#H2:14.188608929503632::0.65)#H1:19.8491740349711::0.4,(((B:12.839184645488634,((((N:1.2,O:1.2)S19:0.8735999999999999,K:2.0736)S18:3.0861803519999995,#H2:0.859963392::0.35)S17:1.0319560703999997,E:6.191736422399999)S16:6.647448223088635)S15:2.5678369290977265,(H:2.9859839999999997,(J:2.48832,I:2.48832)S14:0.4976639999999999)S13:12.42103757458636)S12:16.54097836247591,#H1:13.459574047558643::0.6)S11:6.389599987412456)S10:16.868543966768875,((A:26.62333328088523,((((M:1.44)#H3:0.28800000000000003::0.29999999999999993,L:1.728)S9:8.971320537907197,(#H3:5.990083706879998::0.7000000000000001,(P:1.0)#H4:6.430083706879999::0.2499999999999999)S8:3.2692368310271975)S7:11.486790529497162,(D:8.916100448255998,C:8.916100448255998)S6:13.27001061914836)S5:4.437222213480872)S4:19.381786628484445,(G:3.5831807999999996,#H4:2.5831807999999996::0.7500000000000001)S3:42.421939109369674)S2:9.20102398187393)S1:44.793856108756394);";
        Network truenet = Networks.readNetwork(trueString);

//        Map<Set<String>, Network> inferred_map = inferred_quartets(resultPath);
//        List<Tuple3<Network,String,Double>> newnetlist = new ArrayList<>();
//
//        for (Network net : inferred_map.values()) {
//            newnetlist.add(new Tuple3<>(net, "", 1.0));
////                if(net.findNode("I") == null) continue;
////                if(net.findNode("T") == null) continue;
////                if(net.findNode("B") == null) continue;
//            System.out.println(net.toString());
//        }
//        SuperNetwork4 superNetwork4 = new SuperNetwork4(newnetlist);
//        Network result = superNetwork4.compute();
//        System.out.println(result.toString());
//        System.out.println(Networks.hasTheSameTopology(result, truenet));
//        Tuple3<Network, Network, Double> tuple3 = Pipeline.CheckWithTrueNetwork(result, truenet);
//        System.out.println(tuple3.Item3);

        String mergedstring = "(((K:2.119170788125523E-4)I4#H1:0.00290095019804064::0.5617656110869003,((G:3.7134281459115504E-4)I8#H2:6.987055190779999E-4::0.003768494201859096,((Z:3.366647038764408E-4,(P:1.7391959263406277E-5)I22#H3:3.192727446130345E-4::0.9722551367036831)I13:9.589955848104313E-5,(E:4.149324658270239E-5)I14#H4:3.9107101577478154E-4::0.60034828111886)I7:6.37484071311671E-4)I3:0.0020428189431840374)I1:0.7095673650611073,(((I:0.0010163295698608933,((J:1.7439924534806613E-4)I28#H5:4.6387266727052583E-4::0.013692670976666517,((I22#H3:5.048286850778905E-5::0.02774486329631687,H:6.787482777119533E-5)I32:2.4182642398402085E-4,I14#H4:2.682080051725138E-4::0.39965171888114004)I27:3.285706608633758E-4)I19:3.780576572423014E-4)I11:0.143451166452831,((((O:0.0015729999933138763,N:0.0015729999933138763)I33:1.145329908659972E-4,I28#H5:0.0015131337388318074::0.9863073290233335)I29:0.0019484551013015962,((I8#H2:9.212882133219404E-4::0.9962315057981409,(M:7.42161800639618E-4,F:7.42161800639618E-4)I31:5.504692272734774E-4)I25:5.086950076606191E-4)I16#H6:0.0018346620499077552::0.2880109191533754)I21:0.002529526871677449,(I4#H1:4.84412966824446E-4::0.43823438891309974,B:6.963300456369983E-4)I20:0.005469184911521921)I12:0.138301981065533)I6:0.3086535051609408,((I16#H6:0.006552553860610178::0.7119890808466246,(((D:4.7808152693661825E-4)I17#H7:1.701022494851514E-5::0.9879713938586498,(C:2.1375097633185785E-4)I30#H8:2.813407755532755E-4::0.9944407651431842)I23:0.003822068901056106,(A:1.1065279166858547E-4)I24#H9:0.004206507861272654::0.5992780998320164)I15:0.004036719243242653)I9:0.4081835175513336,(I17#H7:0.20556419174851184::0.012028606141350207,((L:6.143927428186824E-4,I30#H8:4.006417664868245E-4::0.00555923485681573)I26:0.0025058757931726654,I24#H9:0.0030096157443227625::0.40072190016798354)I18:0.20292200473945712)I10:0.21049512417206903)I5:0.03658360373611519)I2:0.25955923115432783)I0;";
        Network mergednet = Networks.readNetwork(mergedstring);
//        System.out.println(Networks.hasTheSameTopology(mergednet, truenet));
        Tuple3<Network, Network, Double> tuple3 = Pipeline.CheckWithTrueNetwork(mergednet, truenet);
        System.out.println(tuple3.Item3);

    }


    public static void main(String[] args) {
//        Network truenet = Networks.readNetwork(truestring);
//
//        List<String> taxalist = new ArrayList<>();
//        // get taxa list
//        for (Object o: truenet.getLeaves()){
//            NetNode n = (NetNode) o;
//            taxalist.add(n.getName());
//        }
//        List<List<String>> result = generateCombinations(taxalist, 4);
//        for (List<String> l:result){
//            System.out.println(l);
//        }
//        System.out.println(result.size());

        compare_true_infer();
//        check_compute();

//        Network truenet = Networks.readNetwork("(Z:100.0,((((F:4.299816959999999)#H2:14.188608929503632::0.65)#H1:19.8491740349711::0.4,(((B:12.839184645488634,((((N:1.2,O:1.2)S19:0.8735999999999999,K:2.0736)S18:3.0861803519999995,#H2:0.859963392::0.35)S17:1.0319560703999997,E:6.191736422399999)S16:6.647448223088635)S15:2.5678369290977265,(H:2.9859839999999997,(J:2.48832,I:2.48832)S14:0.4976639999999999)S13:12.42103757458636)S12:16.54097836247591,#H1:13.459574047558643::0.6)S11:6.389599987412456)S10:16.868543966768875,((A:26.62333328088523,((((M:1.44)#H3:0.28800000000000003::0.29999999999999993,L:1.728)S9:8.971320537907197,(#H3:5.990083706879998::0.7000000000000001,(P:1.0)#H4:6.430083706879999::0.2499999999999999)S8:3.2692368310271975)S7:11.486790529497162,(D:8.916100448255998,C:8.916100448255998)S6:13.27001061914836)S5:4.437222213480872)S4:19.381786628484445,(G:3.5831807999999996,#H4:2.5831807999999996::0.7500000000000001)S3:42.421939109369674)S2:9.20102398187393)S1:44.793856108756394);");
//        Network infernet = Networks.readNetwork("((((((J:0.022477947716178236,I:0.022477947716178236):0.005397749363972762,H:0.027875697080151):0.1108981039133852,(((((O:0.010740898031913573,N:0.010740898031913573):0.0069016503369853845,K:0.017642548368898957):0.028895188592877704,(F:0.038180488124454626)#H1:0.008357248837322036::0.2903500929622906):0.011285498661766688,E:0.05782323562354335):0.059237849493109763,B:0.11706108511665311):0.021712715876883087):0.010285812151233698,#H1:0.11087912502031527::0.7096499070377094):0.2943976037547026,((G:0.029966324100279308,(P:0.004935246053082393)#H2:0.025031078047196915::0.7676546650822365):0.3428876231473879,((((#H2:0.06450189212344853::0.2323453349177635,(M:0.014111055616117188)#H3:0.055326082560413735::0.6697157428923726):0.029620525833882458,(L:0.0164615151767947,#H3:0.0023504595606775123::0.33028425710762743):0.08259614883361868):0.09732613395363696,(D:0.08308034932794546,C:0.08308034932794546):0.11330344863610488):0.044048918615180344,A:0.24043271657923068):0.13242123066843656):0.07060326965180524):0.28541539571358054,Z:0.728872612613053);");
//        Tuple3<Network, Network, Double> tuple = Pipeline.CheckWithTrueNetwork(infernet, truenet);
//        Tuple3<Network, Network, Double> tuple2 = Pipeline.CheckWithTrueBackbone(infernet, truenet);
//
//        System.out.println(tuple.Item3);
//        System.out.println(tuple2.Item3);
    }
}
