package edu.rice.cs.bioinfo.manuscriptsupport.yununnamed2012;

import java.util.*;
import java.io.*;

import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.network.rea.ReaHillClimberSteepestAscent;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.MDCInference_DP;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.*;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.*;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.*;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.MDCOnNetworkYF;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.printing.*;
import edu.rice.cs.bioinfo.library.phylogenetics.*;
import edu.rice.cs.bioinfo.library.phylogenetics.graphadapters.jung.DirectedGraphToGraphAdapter;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.parsers.antlr.ast.RichNewickReaderAST_ANTLR;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.Func4;
import edu.rice.cs.bioinfo.library.phylogenetics.scoring.network.acceptancetesting.Jung.MDCOnNetworkYFFromRichNewickJung;
import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.rea.*;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.HillClimbResult;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.printing.HybridNodeType;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.RichNewickReadResult;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.Networks;
//import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;



/**
 * Created by IntelliJ IDEA.
 * User: yy9
 * Date: 5/9/12
 * Time: 2:04 PM
 * To change this template use File | Settings | File Templates.
 */
public class Simulation201205 extends MDCOnNetworkYFFromRichNewickJung {
    DirectedGraphToGraphAdapter<String,PhyloEdge<String>> _network;
    List<String> _optimalNetworks;
    List<double[]> _optimalProbabilities;
    //List<Integer> _optimalScores;
    int _optimalScore;
    int iteration = 1;

    //String base_path = "/Users/yyu/research/2010-12/experiment/";
    //String result_path = "/Users/yyu/research/2012-03/experiment/";
    //String base_path = "/research/2010-12/experiment/";
    //String result_path = "/research/2012-03/experiment/";
    //String base_path = "/projects/nakhleh/yy9/2010-12/experiment/";
    String base_path = "/users/yy9/project/2010-12/experiment/";
    String result_path = "/users/yy9/project/2012-03/experiment/";
    //String[] cases = {"Laura/","Luay/","Independent1/","Independent2/","Dependent/","Extinction1/","Extinction3/"};
    //String MS_path = "/research/tools/MS/msdir/ms";
    String MS_path = "/Users/yyu/research/tools/MS/msdir/ms";
    //String MrBayes_path = "/Users/yyu/research/tools/MyBayes/";
    int[] loci = {10,25,50,100,500,1000,2000};
    String[] loci_path = {"10loci/","25loci/","50loci/","100loci/","500loci/","1000loci/","2000loci/"};
    double[] interval = {1,2};
    String[] interval_path = {"interval1/","interval2/"};
    String[] allele_path = {"1allele/","2allele/","4allele/","8allele/","16allele/"};
    int[] alleles = {1,2,4,8,16};

    String[] networkStrings = {"((A,((B,C))X#H1),(D,X#H1))R;","(((A,(B)X#H1),(C,X#H1)),((D,(E)Y#H2),(F,Y#H2)))R;","((A,((B,(C)Y#H2),(D,Y#H2))X#H1),(X#H1,E))R;",
            "((A,(B,(C,(D)Y#H2))X#H1),((F,(E,Y#H2)),X#H1))R;","(((A,((E)Z#H3)X#H1),(B,X#H1)),((C,(Z#H3)Y#H2),(D,Y#H2)))R;","(((A,((E)Z#H3)X#H1),(B,(Z#H3)Y#H2)),((C,X#H1),(D,Y#H2)))R;"};
    String[] case_path = {"Luay/","Independent1/","Indepedent2/","Dependent/","Extinction1/","Extinction4/"};
    String[][] probabilities = {{"0.0","0.3","0.5"},{"0.0_0.5","0.3_0.3","0.5_0.5"},{"0.0_0.5","0.3_0.3","0.5_0.5"},
            {"0.0_0.5","0.3_0.3","0.5_0.0","0.5_0.5","0.5_1.0"},{"0.0_0.0_0.5","0.3_0.3_0.3","0.5_0.5_0.5"},{"0.0_0.0_0.5","0.3_0.3_0.3","0.5_0.5_0.5"}};


    public static void main(String[] args) {
        Simulation201205 sim = new Simulation201205();
        //sim.inferSingleAlleleNetworks(Integer.parseInt(args[0]));
        // sim.estimateYeast();
        sim.test();
        //sim.estimateLuay();
        //sim.calculateAccuracy();
        //sim.inferSingleAlleleNetworks(1);
    }

    public Simulation201205()
    {
        super(new RichNewickReaderAST_ANTLR());
    }

    private void test(){
        //_network = makeNetwork("((A:2,((B:1,C:1)K:1)X#1:1::0.3)J:1,(D:2,X#1:1::0.7)L:1)M;");

        long start = System.currentTimeMillis();
        int r=2;
        //_network = makeNetwork(networkStrings[3]);
        ArrayList<Tree> geneTrees = readSTITrees("/Users/yy9/temp/gtrees89");
        Map<String, List<String>> species2alleles = readMapFile1("/Users/yy9/temp/map");
        getMDCStartTree("/Users/yy9/temp/gtrees89","/Users/yy9/temp/map");
        //_network = makeNetwork("((10#H1,((E,(D)10#H1)9,C)79561788-480a-48cf-a494-e6295a046073)c31a4508-0239-4a5c-aeac-f28a6b561109,(B,A)817fc259-6a85-43c2-8d38-78540d84198b)R;");
        _optimalNetworks = new ArrayList<String>();
        _optimalScore = Integer.MAX_VALUE;
        _optimalProbabilities = new ArrayList<double[]>();
        //geneTrees.add(makeNetwork("((((Scer:1,Spar:1),Smik:1):1,Skud:1):1,Sbay:1);"));


        ReticulateEdgeAddition<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>> reaStrategy = new ReticulateEdgeAdditionInPlace<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>>(makeNode, makeEdge);
        ReaHillClimberSteepestAscent<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>,Integer> searcher =
                new ReaHillClimberSteepestAscent<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>,Integer>(reaStrategy, false);

        Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Integer> scorer = getScoreFunction(geneTrees, species2alleles);
        Comparator<Integer> comparator = getScoreComparator();
        HillClimbResult<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,Integer> result = searcher.search(_network, scorer, comparator, r); // search starts here

        int index = _optimalNetworks.size();
        //System.out.println(result.NumIterationsPerformed);
        if(result.GenerationCount < r){
            index = 1;
            int first = 0;
            for(String net: _optimalNetworks){
                int numHybrid = 0;
                for(int i=1; i<=(result.GenerationCount+1); i++){
                    if(net.contains("#H"+i)){
                        numHybrid = i;
                    }
                    else{
                        break;
                    }
                }

                if(index == 1){
                    first = numHybrid;
                }
                else{
                    if(numHybrid > first){

                        break;
                    }
                }
                index ++;
            }
            index--;

        }
        System.out.println((System.currentTimeMillis()-start)/100.0);


        Iterator<String> inferredNetworks = _optimalNetworks.iterator();
        for(int i=0; i<index; i++){
            System.out.println("Inferred network:  "+ inferredNetworks.next().replace(":1.0:100:1.0",""));
        }

/*
         List<Tree> newgts = readSTITrees("/Users/yy9/temp/gtrees0");
        Map<String, List<String>> species2alleles = readMapFile1("/Users/yy9/temp/map");
        try{

        RichNewickReaderAST_ANTLR reader = new RichNewickReaderAST_ANTLR();
        NetworkFactoryFromRNNetwork transformer = new NetworkFactoryFromRNNetwork();
        RichNewickReadResult<Networks> readResult = reader.read(new ByteArrayInputStream(networkStrings[4].getBytes()));
        Network trueNetwork = transformer.makeNetwork(readResult.getNetworks().Networks.iterator().next());

            long start = System.currentTimeMillis();
            edu.rice.cs.bioinfo.programs.phylonet.algos.network.MDCOnNetworkYF mdc = new edu.rice.cs.bioinfo.programs.phylonet.algos.network.MDCOnNetworkYF();

        mdc.countExtraCoal(trueNetwork,newgts,species2alleles);
            double timeUsed = (System.currentTimeMillis()-start)/1000.0;
            System.out.println(timeUsed);
        }catch(Exception e){
            System.err.println(e.getMessage());
            e.getStackTrace();
          }
*/
    }

    private void estimateYeast(){
        String path = "/research/2010-12/experiment/Yeast/";
        readNetwork(path + "tree1");
        System.out.println(networkToString());
        ArrayList<Tree> geneTrees = new ArrayList<Tree>();

        ArrayList<Tree> originalTrees = readSTITrees(path+"Trees_MP/yeast_pruned.trees");
        int[] resolutionNum = new int[originalTrees.size()];
        int index = 0;
        for(Tree gt: originalTrees){
            if(Trees.isBinary(gt)){
                resolutionNum[index] = 1;
                geneTrees.add(gt);
            }
            else{
                int count = 0;
                for(Tree bgt: Trees.getAllBinaryResolution(gt)){
                    geneTrees.add(bgt);
                    count++;
                }
                resolutionNum[index] = count;
            }
            index++;
        }

        _optimalNetworks = new ArrayList<String>();
        _optimalScore = Integer.MAX_VALUE;
        //_optimalScores = new ArrayList<Integer>();
        ReticulateEdgeAddition<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>> reaStrategy = new ReticulateEdgeAdditionInPlace<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>>(makeNode, makeEdge);
        ReaHillClimberSteepestAscent<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>,Integer> searcher =
                new ReaHillClimberSteepestAscent<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>,Integer>(reaStrategy, false);

        Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Integer> scorer = getScoreFunction(geneTrees, null, resolutionNum);
        Comparator<Integer> comparator = getScoreComparator();
        HillClimbResult<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,Integer> result = searcher.search(_network, scorer, comparator, 2); // search starts here

        //System.out.println(networkToString());

        for(String net: _optimalNetworks){
            System.out.println(net);
        }
        System.out.println();
        //System.out.println(_optimalScores);

        System.out.println(result.BestExaminedScore);
    }

    private void inferSingleAlleleNetworks(int id){
        try{
            for(int r=1; r<=3; r++){
                for(int t=0; t<1; t++){
                    Map<String,List<String>> species2alleles = readMapFile1(base_path + case_path[id] + interval_path[t] + "map");
                    for(int p=1; p<2; p++){
                        for(int l=0; l<loci.length; l++){
                            String path = base_path + case_path[id] +  interval_path[t] + probabilities[id][p] + "/" + loci_path[l];
                            File resultDir = new File(result_path + case_path[id] +  interval_path[t] + probabilities[id][p] + "/" + loci_path[l] + "result/");
                            resultDir.mkdirs();
                            System.out.println(path);

                            for(int k=0; k<100; k++){
                                getMDCStartTree(path + "gtrees" + k, base_path + case_path[id] +  interval_path[t] + "map");
                                ArrayList<Tree> geneTrees = readSTITrees(path + "gtrees" + k);
                                long startTime = System.currentTimeMillis();
                                _optimalNetworks = new ArrayList<String>();
                                _optimalScore = Integer.MAX_VALUE;
                                _optimalProbabilities = new ArrayList<double[]>();
                                ReticulateEdgeAddition<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>> reaStrategy = new ReticulateEdgeAdditionInPlace<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>>(makeNode, makeEdge);
                                ReaHillClimberSteepestAscent<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>,Integer> searcher = new ReaHillClimberSteepestAscent<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>,Integer>(reaStrategy, false);
                                Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Integer> scorer = getScoreFunction(geneTrees, species2alleles);
                                Comparator<Integer> comparator = getScoreComparator();
                                HillClimbResult<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,Integer> result = searcher.search(_network, scorer, comparator, r); // search starts here
                                double usedTime = (System.currentTimeMillis() - startTime) / 1e3;

                                int index = _optimalNetworks.size();
                                if(result.GenerationCount < r){
                                    index = 1;
                                    int first = 0;
                                    for(String net: _optimalNetworks){
                                        int numHybrid = 0;
                                        for(int i=1; i<=(result.GenerationCount+1); i++){
                                            if(net.contains("#H"+i)){
                                                numHybrid = i;
                                            }
                                            else{
                                                break;
                                            }
                                        }

                                        if(index == 1){
                                            first = numHybrid;
                                        }
                                        else{
                                            if(numHybrid > first){

                                                break;
                                            }
                                        }
                                        index ++;
                                    }
                                    index--;

                                }

                                File resultFile = new File(result_path + case_path[id] +  interval_path[t] + probabilities[id][p] + "/" + loci_path[l] + "result/result_mdc_"+r);
                                if(k==0){
                                    new File(result_path + case_path[id] +  interval_path[t] + probabilities[id][p] + "/" + loci_path[l] + "result/result_mdc").delete();
                                    resultFile.delete();
                                }
                                BufferedWriter bw = new BufferedWriter(new FileWriter(resultFile, true));
                                bw.append("Rep"+ k + ":");
                                bw.newLine();
                                bw.append("Number of hybrid nodes:  " + result.GenerationCount);
                                bw.newLine();
                                bw.append("Extra lineages:  "+ _optimalScore);
                                bw.newLine();
                                Iterator<String> inferredNetworks = _optimalNetworks.iterator();
                                Iterator<double[]> inferredProbabilities = _optimalProbabilities.iterator();
                                for(int i=0; i<index; i++){
                                    bw.append("Inferred network:  "+ inferredNetworks.next().replace(":1.0:1.0:1.0",""));
                                    bw.newLine();
                                    bw.append("Inferred probability:  "+ inferredProbabilities.next()[0]);
                                    bw.newLine();
                                }
                                bw.append("Running time: " + usedTime);
                                bw.newLine();
                                bw.newLine();
                                bw.close();

                            }

                        }

                    }
                }
            }
        }catch(Exception e){
            System.err.println(e.getMessage());
            e.getStackTrace();
        }
    }


    private void inferMultipleAllelesNetworks(int id){
        try{
            for(int r=1; r<=3; r++){
                for(int t=0; t<2; t++){
                    for(int a=0; a<allele_path.length; a++){
                        Map<String,List<String>> species2alleles = readMapFile1(base_path + case_path[id] + interval_path[t] + allele_path[a] + "map");
                        for(int p=0; p<probabilities[id].length; p++){
                            for(int l=0; l<loci.length; l++){
                                String path = base_path + case_path[id] +  interval_path[t] + allele_path[a] + probabilities[id][p] + "/" + loci_path[l];
                                File resultDir = new File(result_path + case_path[id] +  interval_path[t] + allele_path[a]  + probabilities[id][p] + "/" + loci_path[l] + "result/");
                                resultDir.mkdirs();
                                System.out.println(path);
                                for(int k=0; k<100; k++){
                                    getMDCStartTree(path + "gtrees" + k, base_path + case_path[id] +  interval_path[t] + allele_path[a]  + "map");
                                    ArrayList<Tree> geneTrees = readSTITrees(path + "gtrees" + k);
                                    long startTime = System.currentTimeMillis();
                                    _optimalNetworks = new ArrayList<String>();
                                    _optimalScore = Integer.MAX_VALUE;
                                    _optimalProbabilities = new ArrayList<double[]>();
                                    ReticulateEdgeAddition<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>> reaStrategy = new ReticulateEdgeAdditionInPlace<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>>(makeNode, makeEdge);
                                    ReaHillClimberSteepestAscent<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>,Integer> searcher = new ReaHillClimberSteepestAscent<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>,Integer>(reaStrategy, false);
                                    Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Integer> scorer = getScoreFunction(geneTrees, species2alleles);
                                    Comparator<Integer> comparator = getScoreComparator();
                                    HillClimbResult<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,Integer> result = searcher.search(_network, scorer, comparator, r); // search starts here
                                    double usedTime = (System.currentTimeMillis() - startTime) / 1e3;

                                    int index = _optimalNetworks.size();
                                    if(result.GenerationCount < r){
                                        index = 1;
                                        int first = 0;
                                        for(String net: _optimalNetworks){
                                            int numHybrid = 0;
                                            for(int i=1; i<=(result.GenerationCount+1); i++){
                                                if(net.contains("#H"+i)){
                                                    numHybrid = i;
                                                }
                                                else{
                                                    break;
                                                }
                                            }

                                            if(index == 1){
                                                first = numHybrid;
                                            }
                                            else{
                                                if(numHybrid > first){

                                                    break;
                                                }
                                            }
                                            index ++;
                                        }
                                        index--;

                                    }

                                    File resultFile = new File(result_path + case_path[id] +  interval_path[t] + allele_path[a] + probabilities[id][p] + "/" + loci_path[l] + "result/result_mdc_"+r);
                                    if(k==0){
                                        new File(result_path + case_path[id] +  interval_path[t] + allele_path[a] + probabilities[id][p] + "/" + loci_path[l] + "result/result_mdc").delete();
                                        resultFile.delete();
                                    }
                                    BufferedWriter bw = new BufferedWriter(new FileWriter(resultFile, true));
                                    bw.append("Rep"+ k + ":");
                                    bw.newLine();
                                    bw.append("Number of hybrid nodes:  " + result.GenerationCount);
                                    bw.newLine();
                                    bw.append("Extra lineages:  "+ _optimalScore);
                                    bw.newLine();
                                    Iterator<String> inferredNetworks = _optimalNetworks.iterator();
                                    Iterator<double[]> inferredProbabilities = _optimalProbabilities.iterator();
                                    for(int i=0; i<index; i++){
                                        bw.append("Inferred network:  "+ inferredNetworks.next());
                                        bw.newLine();
                                        bw.append("Inferred probability:  "+ inferredProbabilities.next()[0]);
                                        bw.newLine();
                                    }
                                    bw.append("Running time: " + usedTime);
                                    bw.newLine();
                                    bw.newLine();
                                    bw.close();

                                }

                            }
                        }

                    }
                }
            }
        }catch(Exception e){
            System.err.println(e.getMessage());
            e.getStackTrace();
        }
    }


    private void calculateAccuracy(int id){

        try{
            for(int t=0; t<2; t++){
                for(int p=1; p<probabilities[id].length; p++){
                    for(int l=0; l<loci.length; l++){
                        File resultFile = new File(result_path + case_path +  interval_path[t] + probabilities[id][p] + "/" + loci_path[l] + "result/result_mdc");
                        BufferedReader br = new BufferedReader(new FileReader(resultFile));
                        String line;
                        double[] treeDistance = new double[3];
                        double[] clusterDistance = new double[3];
                        int xl = 0;
                        while((line=br.readLine())!=null){
                            if(line.startsWith("Rep")){
                                double dist1[] = new double[3];
                                double dist2[] = new double[3];
                                int count = 0;
                                while(!(line=br.readLine().trim()).equals("")){
                                    if(line.startsWith("Extra lineages:")){
                                        xl += Integer.parseInt(line.substring(line.indexOf(":")+1).trim());
                                    }
                                    if(line.startsWith("Inferred network:")){
                                        String netExp = line.substring(line.indexOf(":")+1).trim();
                                        RichNewickReaderAST_ANTLR reader = new RichNewickReaderAST_ANTLR();
                                        NetworkFactoryFromRNNetwork transformer = new NetworkFactoryFromRNNetwork();

                                        RichNewickReadResult<Networks> readResult = reader.read(new ByteArrayInputStream(netExp.replaceAll(":1.0","").getBytes()));
                                        //System.out.println("here");
                                        Network inferredNetwork = transformer.makeNetwork(readResult.getNetworks().Networks.iterator().next());
                                        readResult = reader.read(new ByteArrayInputStream(networkStrings[id].getBytes()));
                                        Network trueNetwork = transformer.makeNetwork(readResult.getNetworks().Networks.iterator().next());
                                        int index = 0;
                                        for(double d: edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.computeTreeDistance(inferredNetwork, trueNetwork)){
                                            dist1[index++] += d;
                                        }
                                        index = 0;
                                        for(double d: edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.computeClusterDistance(inferredNetwork, trueNetwork)){
                                            dist2[index++] += d;
                                        }
                                        count++;

                                    }
                                }
                                for(int i=0; i<3; i++){
                                    treeDistance[i] += dist1[i]/count;
                                }
                                for(int i=0; i<3; i++){
                                    clusterDistance[i] += dist2[i]/count;
                                }
                            }
                        }

                        br.close();

                        System.out.println(result_path + case_path +  interval_path[t] + probabilities[id][p] + "/" + loci_path[l]);
                        System.out.println("Tree distance: " + treeDistance[0]/100 + " " + treeDistance[1]/100 + " " + treeDistance[2]/100);
                        System.out.println("Cluster distance: " + clusterDistance[0]/100 + " " + clusterDistance[1]/100 + " " + clusterDistance[2]/100);
                        System.out.println("XL: " + (xl/(loci[l]*100.0)));
                        System.out.println();
                    }
                }
            }
        }catch(Exception e){
            System.err.println(e.getMessage());
            e.getStackTrace();
        }
    }


    private ArrayList<Graph<String,PhyloEdge<String>>> readTrees(String fileName){
        ArrayList<Graph<String,PhyloEdge<String>>> geneTrees = new ArrayList<Graph<String,PhyloEdge<String>>>();
        try{
            BufferedReader br = new BufferedReader(new FileReader(fileName));
            String line;
            while((line = br.readLine()) != null){
                geneTrees.add(makeNetwork(line));
            }
            br.close();
        }catch(Exception e){
            System.err.println(e.getMessage());
            e.getStackTrace();
        }
        return geneTrees;
    }

    private ArrayList<Tree> readSTITrees(String fileName){
        ArrayList<Tree> geneTrees = new ArrayList<Tree>();
        try{
            BufferedReader br = new BufferedReader(new FileReader(fileName));
            String line;
            while((line = br.readLine()) != null){
                NewickReader nr = new NewickReader(new StringReader(line));
                geneTrees.add(nr.readTree());
            }
            br.close();
        }catch(Exception e){
            System.err.println(e.getMessage());
            e.getStackTrace();
        }
        return geneTrees;
    }


    private void readNetwork(String fileName){
        try{
            BufferedReader br = new BufferedReader(new FileReader(fileName));
            String line = br.readLine();
            br.close();
            _network = makeNetwork(line);
        }catch(Exception e){
            System.err.println(e.getMessage());
            e.getStackTrace();
        }
    }

    private HashMap<String,List<String>> readMapFile1(String fileName){
        HashMap<String,List<String>> taxonMap = null;
        try{
            BufferedReader br = new BufferedReader(new FileReader(fileName));
            String line;
            taxonMap = new HashMap<String,List<String>>();
            while ((line = br.readLine()) != null) {
                String[] mapString = line.trim().split(";");
                for(String s : mapString){
                    String species = s.substring(0,s.indexOf(":")).trim();
                    s = s.substring(s.indexOf(":")+1);
                    String[] alleles = s.split(",");
                    List<String> allelelist = Arrays.asList(alleles);
                    taxonMap.put(species, allelelist);
                }
            }
            br.close();

        }catch(Exception e){
            System.err.println(e.getMessage());
            e.getStackTrace();
        }
        return taxonMap;
    }

    private HashMap<String,String> readMapFile2(String fileName) throws Exception{
        BufferedReader br = new BufferedReader(new FileReader(fileName));
        String line;
        HashMap<String,String> taxonMap = new HashMap<String,String>();
        while ((line = br.readLine()) != null) {
            String[] mapString = line.trim().split(";");
            for(String s : mapString){
                String species = s.substring(0,s.indexOf(":")).trim();
                s = s.substring(s.indexOf(":")+1);
                String[] alleles = s.split(",");
                for(String allele: alleles){
                    taxonMap.put(allele, species);
                }
            }
        }
        br.close();
        return taxonMap;
    }

    private void getMDCStartTree(String gtFile, String mapFile){
        try{
            List<Tree> gts = readSTITrees(gtFile);
            Map<String,String> allele2species = readMapFile2(mapFile);
            MDCInference_DP mdc = new MDCInference_DP();
            Tree st = mdc.inferSpeciesTree(gts, allele2species, false, 1, true, 100, true, -1).get(0)._st;
            ((STINode)st.getRoot()).setName("R");
            System.out.println(st.toString());
            _network = makeNetwork(st.toString());

        }catch(Exception e){
            System.err.println(e.getMessage());
            e.getStackTrace();
        }

    }


    private String networkToString(){
        try{
            StringWriter sw = new StringWriter();
            new RichNewickPrinterCompact<String>().print(true, "R", _getNetworkNodeLabel, _getDestinationNodes, _getNetworkDistanceForPrint, _getSupportForPrint, _getProbabilityForPrint, _getHybridNodeIndexForPrint, _getHybridTypeForPrint, sw);
            sw.flush();
            sw.close();
            return sw.toString();
        }catch (Exception e){
            System.err.println(e.getMessage());
            e.getStackTrace();
        }
        return null;
    }

    private Func1<String,String> _getNetworkNodeLabel  = new Func1<String,String>()
    {
        public String execute(String node) {
            return node;
        }
    };

    private Func1<String, Iterable<String>> _getDestinationNodes = new Func1<String, Iterable<String>>() {
        public Iterable<String> execute(String node) {
            return new GetDirectSuccessors<String,PhyloEdge<String>>().execute(_network, node);
        }
    };

    private Func2<GraphReadOnly<String,PhyloEdge<String>>, PhyloEdge<String>, Double> _getNetworkDistance = new Func2<GraphReadOnly<String, PhyloEdge<String>>, PhyloEdge<String>, Double>() {
        public Double execute(GraphReadOnly<String, PhyloEdge<String>> network, PhyloEdge<String> edge) {
            return edge.getBranchLength();
        }
    };

    private Func2<GraphReadOnly<String,PhyloEdge<String>>, PhyloEdge<String>, Double> _getProbability = new Func2<GraphReadOnly<String, PhyloEdge<String>>, PhyloEdge<String>, Double>() {
        public Double execute(GraphReadOnly<String, PhyloEdge<String>> network, PhyloEdge<String> edge) {
            return edge.getProbabilty();
        }
    };

    private Func2<String, String, String> _getNetworkDistanceForPrint = new Func2<String, String, String>() {
        public String execute(String parent, String child) {
            return _network.getEdge(parent, child).getBranchLength()+"";
        }
    };

    private Func2<String, String, String> _getProbabilityForPrint = new Func2<String, String, String>() {
        public String execute(String parent, String child) {
            return _network.getEdge(parent, child).getProbabilty()+"";
        }
    };

    private Func2<String, String, String> _getSupportForPrint = new Func2<String, String, String>() {
        public String execute(String parent, String child) {
            return _network.getEdge(parent, child).getSupport()+"";
        }
    };

    Func1<String, HybridNodeType> _getHybridTypeForPrint = new Func1<String, HybridNodeType>()
    {
        public HybridNodeType execute(String node)
        {
            int inDegree = new GetInDegree<String,PhyloEdge<String>>().execute(_network, node);
            return inDegree == 2 ? HybridNodeType.Hybridization : null;
        }
    };

    Func1<String,String> _getHybridNodeIndexForPrint = new Func1<String, String>() {
        List<String> hybridNodes = new ArrayList<String>();

        public String execute(String node) {
            int inDegree = new GetInDegree<String,PhyloEdge<String>>().execute(_network, node);
            if(inDegree == 2){
                int index = hybridNodes.indexOf(node) + 1;
                if(index == 0){
                    hybridNodes.add(node);
                    return hybridNodes.size()+"";
                }
                else{
                    return index + "";
                }
            }
            else{
                return null;
            }
        }
    };

    private Func4<String,String,Double,Double,PhyloEdge<String>> _makeNetworkEdge = new Func4<String, String, Double, Double, PhyloEdge<String>>() {
        public PhyloEdge<String> execute(String source, String destination, Double distance, Double probability) {
            PhyloEdge<String> edge = new PhyloEdge<String>(source,destination);
            edge.setBranchLength(distance);
            edge.setProbability(probability);
            return edge;
        }
    };

    private Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Integer> getScoreFunction(final ArrayList<Tree> geneTrees, final Map<String, List<String>> species2alleles){
        return new Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Integer>() {
            public Integer execute(DirectedGraphToGraphAdapter<String,PhyloEdge<String>> network) {
                Network ofNetwork = null;
                System.out.println(networkToString());
                try{
                    RichNewickReaderAST_ANTLR reader = new RichNewickReaderAST_ANTLR();
                    NetworkFactoryFromRNNetwork transformer = new NetworkFactoryFromRNNetwork();
                    RichNewickReadResult<Networks> readResult = reader.read(new ByteArrayInputStream(networkToString().replace(":1.0:100.0:1.0","").getBytes()));
                    ofNetwork = transformer.makeNetwork(readResult.getNetworks().Networks.iterator().next());
                }catch (Exception e){

                }
                MDCOnNetworkYF scorer = new MDCOnNetworkYF();                //System.out.println(networkToString());
                List<Integer> scores = scorer.countExtraCoal(ofNetwork, geneTrees, species2alleles);

                int total = 0;

                for(int score: scores){
                    total += score;
                }

                if(total < _optimalScore){
                    _optimalScore = total;
                    _optimalNetworks.clear();
                    _optimalNetworks.add(networkToString());
                    _optimalProbabilities.clear();
                //    _optimalProbabilities.add(scorer.getHybridProbabilities());
                }
                else if(total == _optimalScore){
                    _optimalNetworks.add(networkToString());
                  //  _optimalProbabilities.add(scorer.getHybridProbabilities());
                }
                return total;
            }
        };
    }


    private Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Integer> getScoreFunction(final ArrayList<Tree> geneTrees, final Map<String, List<String>> species2alleles, final int[] resolutionNum){
        return new Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Integer>() {
            public Integer execute(DirectedGraphToGraphAdapter<String,PhyloEdge<String>> network) {
                Network ofNetwork = null;
                try{
                    RichNewickReaderAST_ANTLR reader = new RichNewickReaderAST_ANTLR();
                    NetworkFactoryFromRNNetwork transformer = new NetworkFactoryFromRNNetwork();
                    RichNewickReadResult<Networks> readResult = reader.read(new ByteArrayInputStream(networkToString().replace(":1.0:100.0:1.0","").getBytes()));
                    ofNetwork = transformer.makeNetwork(readResult.getNetworks().Networks.iterator().next());
                }catch (Exception e){

                }
                MDCOnNetworkYF scorer = new MDCOnNetworkYF();
                //System.out.println(networkToString());
                List<Integer> scores = scorer.countExtraCoal(ofNetwork, geneTrees, species2alleles);
                int total = 0;
                /*
                for(int score: scores){
                    total += score;
                }
                */

                Iterator<Integer> it = scores.iterator();
                for(int num: resolutionNum){
                    int min = Integer.MAX_VALUE;
                    for(int i=0; i<num; i++){
                        min = Math.min(min, it.next());
                    }
                    total += min;
                }

                //System.out.println(total);
                if(total < _optimalScore){
                    _optimalScore = total;
                    _optimalNetworks.clear();
                    _optimalNetworks.add(networkToString());
                }
                else if(total == _optimalScore){
                    _optimalNetworks.add(networkToString());
                }
                //System.out.println(scores);
                return total;
            }
        };
    }


    private Comparator<Integer> getScoreComparator(){
        return new Comparator<Integer>() {
            public int compare(Integer o1, Integer o2)
            {
                System.out.println(iteration++);
                //return 0;
                return Double.compare(o2, o1);
            }
        };
    }



}
