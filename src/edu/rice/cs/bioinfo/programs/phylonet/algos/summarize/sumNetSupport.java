package edu.rice.cs.bioinfo.programs.phylonet.algos.summarize;
/*
 *@ClassName: sumNetSupport
 *@Description
 *@Author: Zhen Cao
 *@Date:  2019-08-05 14:18
 *@Version: 1.0
 */
//import edu.rice.cs.bioinfo.programs.phylonet.algos.trinetInference.summarizeTrinets;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.library.programming.Tuple3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.treeAugment.treeRoot;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;
import java.util.List;

public class sumNetSupport {
    private int _chainlen = 1100000;
    private int _burnin =  10000;
    private int _sample_freq = 1000;
    private static String _outDir = "";

    public sumNetSupport(int chainlen, int burnin, int sample_freq, String outDir){
        _chainlen = chainlen;
        _burnin = burnin;
        _sample_freq = sample_freq;
        _outDir = outDir;
    }

    public sumNetSupport(String outDir){

        _outDir = outDir;
    }


    public static class Out {
        public Out(String file){
            try {

                PrintStream print = new PrintStream(new File(_outDir+file));
                System.setOut(print);
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
        }

    }

    public Map feature(int method, Map<Network, Double> networkMap){
        switch (method){
            case 0:
                Map<Tree, Double> maxtreemap = maxTree.summarize(networkMap);
                System.out.println("---------------------max tree---------------------");
                for(Tree t: maxtreemap.keySet()){
                    System.out.println(t.toString()+"\tcount:"+maxtreemap.get(t).toString());
                }
                System.out.println("-------------------------End------------------------");
                return maxtreemap;
            case 1:
                Map<Network, Double> backboneMap = backboneNetwork.summarize(networkMap);
                System.out.println("---------------------backbone network---------------------");
                for(Network backnet: backboneMap.keySet()){
                    System.out.println(backnet.toString()+"\t\t"+backboneMap.get(backnet).toString());
                }
                System.out.println("-------------------------End------------------------");
                return backboneMap;
            case 2:
                Map<Tree, Double> majortreemap = majorTree.summarize(networkMap);
                System.out.println("---------------------major tree---------------------");
                for(Tree t: majortreemap.keySet()){
                    System.out.println(t.toString()+"\t\tcount:"+majortreemap.get(t).toString());
                }
                System.out.println("-------------------------End------------------------");
                return majortreemap;
            case 3:
                Map<Network, Double> decomposedtreemap = decomposedTrees.summarize(networkMap);
                System.out.println("---------------------decomposed tree---------------------");
                for(Network t: decomposedtreemap.keySet()){
                    System.out.println(t.toString()+"\t\tcount:"+decomposedtreemap.get(t).toString());
                }
                System.out.println("-------------------------End------------------------");
                return decomposedtreemap;
            case 4:
                Map<netNodeTuple, Double> tripart = tripartition2.summarize(networkMap);
                System.out.println("---------------------tripartition---------------------");
                for(netNodeTuple t: tripart.keySet()){
                    System.out.println(t.toString()+"\t\tcount:"+tripart.get(t).toString());
                }
                System.out.println("-------------------------End------------------------");
                return tripart;
            case 5:
                nonEquivalent.summarize(networkMap);
                return null;
            default:
                return null;
        }
    }

//    public void produceResultMCMCGT(String path){
//        List<String> filenames = treeRoot.getFiles(path);
//        List<String> filteredFiles = new ArrayList<>();
//        for(String filename:filenames){
//            if(filename.contains("MCMC_GT")){
//                filteredFiles.add(filename);
//            }
//        }
//        int start = (int)(_burnin / _sample_freq) + 1;
//        int end = (int)(_chainlen / _sample_freq);
//        Map<String, List<String>> file2samples = new HashMap<>();
//        Map<String, Tuple<String, Double>> file2topsample = new HashMap<>();
//
//        Map<String, Map<Network, Double>> networkDoubleMap = summarizeTrinets.GetInputFromFolder_MCMCGT(filteredFiles, start, end, file2samples, file2topsample);
//        for(String filename : networkDoubleMap.keySet()){
//            Map<Network, Double> candidate = networkDoubleMap.get(filename);
//            HashMap<Network, Double> dummy = new HashMap<>();
//            dummy.putAll(candidate);
//            String[] arr = filename.split("/");
//            String of = arr[arr.length-1];
//            Out out = new Out(of);
//            System.out.println("---------------------CANDIDATE----------------------");
//            for (Network net: dummy.keySet()){
//                System.out.println(candidate.get(net)+"\t"+net.toString());
//                if(candidate.get(net) < 0.05){
//                    candidate.remove(net);
//                }
//            }
//
//            if (candidate.size() < 2){
//                continue;
//
//            }
//            System.out.println("----------------------------------------------------");
//            for(int i = 0; i < 5; i++) {
//                dummy.clear();
//                dummy.putAll(candidate);
//                Map f = feature(i, dummy);
//            }
//        }
//    }
//
//    public static void main(String[] args) {
//        String outdir = "/Users/zhen/Desktop/Zhen/research/phylogenetics/bookChapter/data/Reti3_C/sumOut/";
//        String path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/bookChapter/data/Reti3_C/out2/";
//        sumNetSupport sumnetSupport = new sumNetSupport(outdir);
//        sumnetSupport.produceResultMCMCGT(path);
//    }
}
