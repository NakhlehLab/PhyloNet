package edu.rice.cs.bioinfo.programs.phylonet.algos.chi2;
/*
 * @ClassName:   GTProbCalculator
 * @Description:
 * @Author:      Zhen Cao
 * @Date:        1/5/22 1:47 PM
 */

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.*;

public class GTProbCalculator {
    /* Constructor */
    public GTProbCalculator() {

    }

    public static String read_ML_species_net(String path, String outgroup){
        try {
            BufferedReader br = new BufferedReader(new FileReader(path));
            String line;
            Boolean begin = false;
            while((line = br.readLine()) != null) {
                if (line.startsWith("Inferred Network #1")) {
                    begin = true;
                }
                else if (begin) {

                    Network net = Networks.readNetwork(line.trim());
                    System.out.println(net.toString());
                    NetNode root = net.getRoot();
//                    root.
                    NetNode out_node = net.findNode(outgroup);
                    root.removeChild(out_node);
                    NetNode new_root = (NetNode) root.getChildren().iterator().next();
                    net.resetRoot(new_root);
                    out_node.removeItself();
                    Networks.removeBinaryNodes(net);
                    System.out.println(net.toString());
//                    STITree tree = new STITree(line.trim());
//                    tree.removeNode(outgroup);
//                    Trees.removeBinaryNodes(tree);
//                    tree.getRoot().setParentDistance(0);
                    return net.toString();
                }
            }
        }catch (Exception e){
            System.err.println(e);
        }
        return null;
    }


    public static void computeGTProb(String snet, String output_path){
        HashMap<String, List<String>> taxonmap = new HashMap<>();

        Network net = Networks.readNetwork(snet);
        int i = 1;
        for(Object o: net.getLeaves()){
            NetNode node = (NetNode) o;

            List<String> list = new ArrayList<>();
            list.add(String.valueOf(i));
            taxonmap.put(node.getName(), list);
            i += 1;
        }
        Enumerator enu = new Enumerator(3, taxonmap);
        StringBuilder sb = enu.computeLikelihood(net);
        Enumerator.outputResult(output_path, sb);
    }

    //diff = (inferred - true)/true
    public static Map<String, Double> compute_brl_diff(Network net, Map<NetNode, NetNode> map) {
        Map<String, Double> id_brldiff = new HashMap<>();
//        Networks.autoLabelNodes(net);
        for(Object n1 : Networks.postTraversal(net)) {
            NetNode child = (NetNode) n1;
            for(Object n2 : child.getParents()) {
                NetNode parent = (NetNode) n2;
                String id = parent.getName() + "-" + child.getName();
                NetNode ch = map.get(child), par = map.get(parent);
                double br1 = child.getParentDistance(parent); //true
                double br2 = ch.getParentDistance(par);  //inferred
                if (!child.isLeaf()){
                    id_brldiff.put(id, (br2-br1)/br1); // (inferred-true)/inferred
                }


            }
        }
        return id_brldiff;
    }

//    public static Map<String, Double> compute_brl_diff(Network net, Map<NetNode, NetNode> map) {
//        Map<String, Double> id_brldiff = new HashMap<>();
////        Networks.autoLabelNodes(net);
//        for(Object n1 : Networks.postTraversal(net)) {
//            NetNode child = (NetNode) n1;
//            for(Object n2 : child.getParents()) {
//                NetNode parent = (NetNode) n2;
//                String id = parent.getName() + "-" + child.getName();
//                NetNode ch = map.get(child), par = map.get(parent);
//                double br1 = child.getParentDistance(parent); //true
//                double br2 = ch.getParentDistance(par);  //inferred
//                id_brldiff.put(id, br1 - br2);
//
//            }
//        }
//        return id_brldiff;
//    }

    public static void compute_all_gt_prob(String root_path, String scale, boolean obs_truegt,int num_taxa, boolean compute_gtprob, String base_rate, int re_start, int re_end, int num_reti){
//        String root_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/taxa3_tall/036/";
        Network trunet = null;
        Map<String, Network> TRUENETMAP = new HashMap<>();

        if (num_taxa == 5){
            TRUENETMAP.put("short", Networks.readNetwork("(((Q:0.3,R:0.3)I3:0.5,L:0.8)I1:0.2,(G:0.4,C:0.4)I2:0.6);"));
            TRUENETMAP.put("medium", Networks.readNetwork("(((Q:0.6,R:0.6):1,L:1.6):0.4,(G:0.8,C:0.8):1.2);"));
            TRUENETMAP.put("long", Networks.readNetwork("(((Q:1.5,R:1.5):2.5,L:4):1,(G:2,C:2):3);"));
            if (num_reti == 0){
                trunet = TRUENETMAP.get(scale);
            }
            else if(num_reti == 1){
                trunet = Networks.readNetwork("((((Q:1)#H1:0.5::0.7,R:1.5)I3:2.5,(L:2,#H1:1::0.3):2)I1:1,(G:2,C:2)I2:3);");

            }
            else if (num_reti == 2){
                trunet = Networks.readNetwork("((((Q:1)#H1:0.5::0.7,R:1.5)I3:2.5,((L:2,#H1:1::0.3):1)#H2:1.0::0.6):1.0,((G:2,C:2)I2:2,#H2:1.0::0.4):1.0);");

            }
        }
        else if(num_taxa == 3){
            trunet = Networks.readNetwork("(L:2.0,(Q:1.0,A:1.0):1.0);");
        }
        StringBuilder sb_brl = new StringBuilder();
        StringBuilder sb_st_topo = new StringBuilder();
        sb_st_topo.append("replicate,newick,correct\n");
        String outgroup = "Z";
        Map<String, List<Double>> id_brl = new HashMap<>();
        Networks.autoLabelNodes(trunet);
        System.out.println(trunet.toString());
        int cnt = 0;

        // compute expected gene tree probabilities
        for (int i = re_start; i <= re_end; i++){
            String inpath = root_path + scale + "/" + String.valueOf(i) + "/heter/";
            String outpath = root_path + scale + "/" + String.valueOf(i) + "/heter/";

            if (base_rate.equals("")){
                if(obs_truegt){
                    inpath += "ML_true_"+num_reti+".out";
                    outpath += "gt_prob_exp_true_"+num_reti+".csv";

                }
                else{
                    inpath += "ML_iq_"+num_reti+".out";
                    outpath += "gt_prob_exp_iq_"+num_reti+".csv";
                }
            }

            else{
                if(obs_truegt) {
                    inpath += base_rate + "/ML_true_"+num_reti+".out";
                    outpath += base_rate + "/gt_prob_exp_true_"+num_reti+".csv";
                }
                else{
                    inpath += base_rate + "/ML_iq_"+num_reti+".out";
                    outpath += base_rate + "/gt_prob_exp_iq_"+num_reti+".csv";
                }
            }


            String snet = read_ML_species_net(inpath, outgroup);
            Network net = Networks.readNetwork(snet);
            System.out.println(net);
            boolean correct_topology = Networks.hasTheSameTopology(net, trunet);

            Map<NetNode, NetNode> nodemap = Networks.mapTwoNetworks(trunet, net);
            if (nodemap != null){
                cnt += 1;
                Map<String, Double> branch_lengths = compute_brl_diff(trunet, nodemap);
                for (String name: branch_lengths.keySet()){
                    id_brl.putIfAbsent(name, new ArrayList<>());
                    id_brl.get(name).add(branch_lengths.get(name));
                }
            }
            sb_st_topo.append(i+",\""+net.toString()+"\","+correct_topology+"\n");
            System.out.println(i+",\""+net.toString()+"\","+correct_topology);
            if (compute_gtprob && false){
                computeGTProb(snet, outpath);
            }
        }
//        System.out.println(id_brl);
        String st_topo = root_path+scale+"/";

        if (obs_truegt){
            st_topo += "st_topo_truegt_"+base_rate+"_"+num_reti+".csv";
        }
        else if (base_rate.equals("")){
            st_topo += "st_topo_infergt_"+num_reti+".csv";

        }
        else{
            st_topo += "/st_topo_infergt_"+base_rate+"_"+num_reti+".csv";
        }
        if (compute_gtprob) {
            String st_brl_path = root_path+scale + "/";
            if (obs_truegt){
                st_brl_path += "st_brl_truegt_"+base_rate+"_"+num_reti+".csv";
            }
            else if (base_rate.equals("")){
                st_brl_path += "st_brl_infergt_"+num_reti+".csv";

            }
            else{
                st_brl_path += "/st_brl_infergt_"+base_rate+"_"+num_reti+".csv";
            }
            List<String> ids = new ArrayList<>(id_brl.keySet());

            for (String id : ids){
                sb_brl.append(id);
                sb_brl.append(",");
            }
            sb_brl.append("\n");

            for (int i = 0; i < cnt; i++){
                for (String id : ids){
                    sb_brl.append(id_brl.get(id).get(i));
                    sb_brl.append(",");
                }
                sb_brl.append("\n");
            }

            Enumerator.outputResult(st_brl_path, sb_brl);
            Enumerator.outputResult(st_topo, sb_st_topo);
        }
    }

    public static void run_all(String path, int re_end, int num_reti){
//        String path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/compare/net0/";
        String[] scale_array = {"short", "medium", "long"};
//        String[] scale_array = {"short"};
        int[] locus_length_array = {500, 2000};
        int re_start = 1;
//        int re_end = 10;
        int num_taxa = 5;
        System.out.println(path);
        System.out.println(re_end);
        boolean[] gt_list = {true, false};
        for (String scale: scale_array){
            for(int locus_length : locus_length_array){
                for(boolean obs_truegt : gt_list){
                    compute_all_gt_prob(path, scale, obs_truegt, num_taxa, true, String.valueOf(locus_length), re_start, re_end, num_reti);
                }
//                compute_all_gt_prob(path+scale+"/", obs_truegt, 5, true, String.valueOf(locus_length), num_replicate);
            }
        }
    }



    public static void main(String[] args) {

//        compute_all_gt_prob(args[0], Boolean.parseBoolean(args[1]), Integer.parseInt(args[2]), Boolean.parseBoolean(args[3]), args[4], Integer.parseInt(args[5]));

//        String snet = "(((Q:0.3,R:0.3)I3:0.5,L:0.8)I1:0.2,(G:0.4,C:0.4)I2:0.6);";
//        String snet = "(((Q:0.6,R:0.6):1,L:1.6):0.4,(G:0.8,C:0.8):1.2);";
//        String snet = "(((Q:1.5,R:1.5):2.5,L:4):1,(G:2,C:2):3);";
//        String outputpath = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/long/gt_prob_exp_true_st.csv";
//        computeGTProb("(L:1.0,(Q:1.0,(A:1.0,(G:1.0,C:1.0):0.5407287513326626):0.4172513268835793):0.10057388068846826);", "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/mosquito/MAF/network_inference/scripts/consensus/6taxa/gt_probs.csv");
//        computeGTProb("((R:1.0,L:1.0):0.04248124417722194,((A:1.0,(C:1.0,G:1.0):0.5036894677522945):0.4638026225852166,Q:1.0):0.08659831598939814);", "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/mosquito/MAF/network_inference/scripts/consensus/6taxa/gt_probs.csv");
//        computeGTProb("(HmelRef:1.0,(Hhsa:1.0,HeraRef:1.0):1.334189880324706);", "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/butterfly/3_3/gt_probs_exp.csv");
//        computeGTProb("(Hnum:1.0,(HmelRef:1.0,Htim:1.0):0.4323171602384371);", "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/butterfly/3_2/gt_probs_exp.csv");
        run_all(args[0], Integer.parseInt(args[1]), Integer.parseInt(args[2]));
//        String path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/data/simulation/scale_ingroup_half_popsize/net0/";
//        String scale = "short";
//        compute_all_gt_prob(path, scale,false, 5, false, "500", 1, 1, 0);

    }
}

