package edu.rice.cs.bioinfo.manuscriptsupport.yununnamed2012;

import java.util.*;
import java.io.*;

import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.MDCInference_DP;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.*;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.*;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.*;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.printing.*;
import edu.rice.cs.bioinfo.library.phylogenetics.*;
import edu.rice.cs.bioinfo.library.phylogenetics.graphadapters.jung.DirectedGraphToGraphAdapter;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.parsers.antlr.ast.RichNewickReaderAST_ANTLR;
import edu.rice.cs.bioinfo.library.phylogenetics.scoring.network.MDCOnNetworkYF;
import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.Func4;
import edu.rice.cs.bioinfo.library.phylogenetics.scoring.network.acceptancetesting.Jung.MDCOnNetworkYFFromRichNewickJung;
import edu.rice.cs.bioinfo.library.phylogenetics.rearrangement.network.rea.*;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.HillClimbResult;
import edu.rice.cs.bioinfo.library.phylogenetics.search.hillclimbing.network.rea.ReaHillClimber;
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
public class Simulation201205 extends MDCOnNetworkYFFromRichNewickJung{
    DirectedGraphToGraphAdapter<String,PhyloEdge<String>> _network;
    List<String> _optimalNetworks;
    List<double[]> _optimalProbabilities;
    //List<Integer> _optimalScores;
    int _optimalScore;
    //int iteration = 1;

    String base_path = "/Users/yyu/research/2010-12/experiment/";
    String result_path = "/Users/yyu/research/2012-03/experiment/";
    //String base_path = "/research/2010-12/experiment/";
    //String result_path = "/research/2012-03/experiment/";
    String[] cases = {"Laura/","Luay/","Independent1/","Independent2/","Dependent/","Extinction2/","Extinction3/"};
    //String MS_path = "/research/tools/MS/msdir/ms";
    String MS_path = "/Users/yyu/research/tools/MS/msdir/ms";
    //String MrBayes_path = "/Users/yyu/research/tools/MyBayes/";
    int[] loci = {10,25,50,100,500,1000,2000};
    String[] loci_path = {"10loci/","25loci/","50loci/","100loci/","500loci/","1000loci/","2000loci/"};
    double[] interval = {1,2};
    String[] interval_path = {"interval1/","interval2/"};
    String[] allele_path = {"1allele/","2allele/","4allele/","8allele/","16allele/"};
    int[] alleles = {1,2,4,8,16};
    
    String LuayNetwork = "((A,((B,C))X#1),(D,X#1));";


    
    public static void main(String[] args) {
        Simulation201205 sim = new Simulation201205();
        //sim.estimateYeast();
        sim.test();
        //sim.estimateLuay();
        //sim.calculateAccuracy();
    }
    
    public Simulation201205()
    {
        super(new RichNewickReaderAST_ANTLR());
    }
    
    private void test(){    
        //_network = makeNetwork("((A:2,((B:1,C:1)K:1)X#1:1::0.3)J:1,(D:2,X#1:1::0.7)L:1)M;");
        _network = makeNetwork("((((Scer:1,Spar:1),Smik:1):1,Skud:1):1,Sbay:1)R;");
        ArrayList<Graph<String,PhyloEdge<String>>> geneTrees = new ArrayList<Graph<String,PhyloEdge<String>>>();
        Map<String, List<String>> species2alleles = null;

        _optimalNetworks = new ArrayList<String>();
        _optimalScore = Integer.MAX_VALUE;

        _optimalProbabilities = new ArrayList<double[]>();
        geneTrees.add(makeNetwork("((Scer,Spar),(Smik,(Skud,Sbay)));"));

        ReticulateEdgeAddition<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>> reaStrategy = new ReticulateEdgeAdditionInPlace<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>>(makeNode, makeEdge);
        ReaHillClimber<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>> searcher = new ReaHillClimber<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>>(reaStrategy);

        Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Integer> scorer = getScoreFunction(geneTrees, species2alleles);
        Comparator<Integer> comparator = getScoreComparator();
        HillClimbResult<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,Integer> result = searcher.search(_network, scorer, comparator, 1); // search starts here
        //System.out.println(result.LocalOptimum.toString());
        //System.out.println(result.LocalOptimumScore);
        for(String net: _optimalNetworks){
            System.out.println(net);
        }
        System.out.println(_optimalScore);
        System.out.println(result.LocalOptimumScore);
        
    }

    private void estimateYeast(){
        String path = "/research/2010-12/experiment/Yeast/";
        readNetwork(path + "tree2");
        System.out.println(networkToString());
        ArrayList<Graph<String,PhyloEdge<String>>> geneTrees = new ArrayList<Graph<String,PhyloEdge<String>>>();

        ArrayList<Tree> originalTrees = readSTITrees(path+"Trees_MP/yeast_pruned.trees");
        int[] resolutionNum = new int[originalTrees.size()];
        int index = 0;
        for(Tree gt: originalTrees){
            if(Trees.isBinary(gt)){
                resolutionNum[index] = 1;
                geneTrees.add(makeNetwork(gt.toString()));
            }
            else{
                int count = 0;
                for(Tree bgt: Trees.getAllBinaryResolution(gt)){
                    geneTrees.add(makeNetwork(bgt.toString()));
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
        ReaHillClimber<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>> searcher = new ReaHillClimber<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>>(reaStrategy);
        
        Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Integer> scorer = getScoreFunction(geneTrees, null, resolutionNum);
        Comparator<Integer> comparator = getScoreComparator();
        HillClimbResult<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,Integer> result = searcher.search(_network, scorer, comparator, 1); // search starts here
        
        //System.out.println(networkToString());
        
        for(String net: _optimalNetworks){
            System.out.println(net);
        }
        System.out.println(_optimalScore);
        //System.out.println(_optimalScores);

        System.out.println(result.LocalOptimumScore);
    }

    private void estimateLuay(){
        String case_path = "Luay/";
        String[] probability = {"0.0","0.3","0.5"};

        try{
            for(int t=0; t<2; t++){
                Map<String,List<String>> species2alleles = readMapFile1(base_path + case_path + interval_path[t] + "map");
                for(int p=0; p<probability.length; p++){
                    for(int l=0; l<loci.length; l++){
                        String path = base_path + case_path +  interval_path[t] + probability[p] + "/" + loci_path[l];
                        File resultDir = new File(result_path + case_path +  interval_path[t] + probability[p] + "/" + loci_path[l] + "result/");
                        resultDir.mkdirs();
                        System.out.println(path);
                        
                        for(int k=0; k<100; k++){
                            getMDCStartTree(path + "gtrees" + k, base_path + case_path +  interval_path[t] + "map");
                            ArrayList<Graph<String,PhyloEdge<String>>> geneTrees = readTrees(path + "gtrees" + k);                            
                            long startTime = System.currentTimeMillis();
                            _optimalNetworks = new ArrayList<String>();
                            _optimalScore = Integer.MAX_VALUE;
                            _optimalProbabilities = new ArrayList<double[]>();
                            ReticulateEdgeAddition<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>> reaStrategy = new ReticulateEdgeAdditionInPlace<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>>(makeNode, makeEdge);
                            ReaHillClimber<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>,String,PhyloEdge<String>> searcher = new ReaHillClimber<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, String, PhyloEdge<String>>(reaStrategy);
                            Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Integer> scorer = getScoreFunction(geneTrees, species2alleles);
                            Comparator<Integer> comparator = getScoreComparator();
                            searcher.search(_network, scorer, comparator, 1); // search starts here

                            double usedTime = (System.currentTimeMillis() - startTime) / 1e3;
                           
                            File resultFile = new File(result_path + case_path +  interval_path[t] + probability[p] + "/" + loci_path[l] + "result/result_mdc");
                            if(k==0){
                                resultFile.delete();
                            }
                            BufferedWriter bw = new BufferedWriter(new FileWriter(resultFile, true));
                            bw.append("Rep"+ k + ":");
                            bw.newLine();
                            bw.append("Extra lineages:  "+ _optimalScore);
                            bw.newLine();
                            Iterator<String> inferredNetworks = _optimalNetworks.iterator();
                            Iterator<double[]> inferredProbabilities = _optimalProbabilities.iterator();
                            while(inferredNetworks.hasNext()){
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
        }catch(Exception e){
            System.err.println(e.getMessage());
            e.getStackTrace();
        }
    }

    private void calculateAccuracy(){
        String case_path = "Luay/";
        String[] probability = {"0.0","0.3","0.5"};

        try{
            for(int t=0; t<2; t++){
                for(int p=1; p<probability.length; p++){
                    for(int l=0; l<loci.length; l++){
                        File resultFile = new File(result_path + case_path +  interval_path[t] + probability[p] + "/" + loci_path[l] + "result/result_mdc");
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
                                        readResult = reader.read(new ByteArrayInputStream(LuayNetwork.getBytes()));
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

                        System.out.println(result_path + case_path +  interval_path[t] + probability[p] + "/" + loci_path[l]);
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

    private HashMap<String,List<String>> readMapFile1(String fileName) throws Exception{
        BufferedReader br = new BufferedReader(new FileReader(fileName));
        String line;
        HashMap<String,List<String>> taxonMap = new HashMap<String,List<String>>();
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
            Tree st = mdc.inferSpeciesTree(gts, allele2species, false, 1, true, 1, true, -1).get(0)._st;
            ((STINode)st.getRoot()).setName("R");
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

    private Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Integer> getScoreFunction(final ArrayList<Graph<String,PhyloEdge<String>>> geneTrees, final Map<String, List<String>> species2alleles){
        return new Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Integer>() {
            public Integer execute(DirectedGraphToGraphAdapter<String,PhyloEdge<String>> network) {
                MDCOnNetworkYF scorer = new MDCOnNetworkYF();
                //System.out.println(networkToString());
                List<Integer> scores = scorer.countExtraCoal(network, geneTrees, species2alleles, _getNetworkNodeLabel, _getNetworkNodeLabel, _getNetworkDistance, _getProbability, _getNetworkDistance, _getProbability,
                        _makeNetworkEdge, _makeNetworkEdge);
                int total = 0;

                for(int score: scores){
                    total += score;
                }

                if(total < _optimalScore){
                    _optimalScore = total;
                    _optimalNetworks.clear();
                    _optimalNetworks.add(networkToString());
                    _optimalProbabilities.clear();
                    _optimalProbabilities.add(scorer.getHybridProbabilities());
                }
                else if(total == _optimalScore){
                    _optimalNetworks.add(networkToString());
                    _optimalProbabilities.add(scorer.getHybridProbabilities());
                }
                return total;
            }
        };
    }


    private Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Integer> getScoreFunction(final ArrayList<Graph<String,PhyloEdge<String>>> geneTrees, final Map<String, List<String>> species2alleles, final int[] resolutionNum){
        return new Func1<DirectedGraphToGraphAdapter<String,PhyloEdge<String>>, Integer>() {
            public Integer execute(DirectedGraphToGraphAdapter<String,PhyloEdge<String>> network) {
                MDCOnNetworkYF scorer = new MDCOnNetworkYF();
                //System.out.println(networkToString());
                List<Integer> scores = scorer.countExtraCoal(network, geneTrees, species2alleles, _getNetworkNodeLabel, _getNetworkNodeLabel, _getNetworkDistance, _getProbability, _getNetworkDistance, _getProbability,
                        _makeNetworkEdge, _makeNetworkEdge);
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
                System.out.println(o1 + " vs. " + o2);
          //     return 0;
                return Double.compare(o2, o1);
            }
        };
    }



}
