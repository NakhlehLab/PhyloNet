package edu.rice.cs.bioinfo.programs.phylonet.algos.chi2;
/*
 * @ClassName:   Enumerator
 * @Description:
 * @Author:      Zhen Cao
 * @Date:        12/18/21 2:04 PM
 */

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.counting.NetworkTopologiesCounting;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.GeneTreeProbabilityYF;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.NetworkLikelihoodFromGTT;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.NetworkLikelihoodFromGTT_SingleTreePerLocus;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.File;
import java.io.PrintWriter;
import java.util.*;


public class Enumerator {
    private int _maxRounds = 100;
    private int _maxTryPerBranch = 100;
    private double _maxBranchLength = 6;
    private double _improvementThreshold = 0.001;
    private double _Brent1 = 0.01;
    private double _Brent2 = 0.001;
    private int _parallel = 1;
    private boolean _inferBL = false;
    private List<Tree> _treelist = null;
    private HashMap<String, List<String>> _taxonMap = null;
    private HashMap<String, String> _allelemap = null;

    /* Constructor */
    public Enumerator(int ntaxa) {
        _treelist = NetworkTopologiesCounting.enumerateTreeTopologies(ntaxa);
    }


    public Enumerator(int ntaxa, HashMap taxonmap) {
        this(ntaxa);
        _taxonMap = taxonmap;
        _allelemap = new HashMap<>();
        for (String k: _taxonMap.keySet()){
            for (String v: _taxonMap.get(k)){
//                _allelemap.put(v, k+"_0");
                _allelemap.put(v, k);
//
            }
        }
    }

    public StringBuilder computeLikelihood(Network speciesNetwork){
        NetworkLikelihoodFromGTT scoring = new NetworkLikelihoodFromGTT_SingleTreePerLocus();
        scoring.setSearchParameter(_maxRounds, _maxTryPerBranch, _improvementThreshold, _maxBranchLength, _Brent1, _Brent2, _parallel);
        StringBuilder sb = new StringBuilder();
        sb.append("gt,");
        sb.append("prob_exp\n");
        for(Tree t : _treelist) {
            List<List<MutableTuple<Tree,Double>>> gts = new ArrayList<List<MutableTuple<Tree,Double>>>();
            List<MutableTuple<Tree, Double>> gtsForOneLocus = new ArrayList<>();
            gtsForOneLocus.add(new MutableTuple<Tree, Double>(t, 1.0));
            gts.add(gtsForOneLocus);
            double score = scoring.computeLikelihood(speciesNetwork, gts, _taxonMap, _inferBL);
//            System.out.println(t.toString());
//            System.out.format("%f %f %n", score, Math.exp(score));

            for(TNode tn: t.getRoot().getLeaves()){
                STINode node = (STINode) tn;
//                System.out.println(node.getName());
                String newname  = _allelemap.get(node.getName());
                node.setName(newname);
            }

            sb.append("\""+t.toString()+"\"");
            sb.append(",");
            sb.append(Math.exp(score));
            sb.append("\n");
        }
//        System.out.println(sb.toString());
        return sb;
    }

    public static void outputResult(String output_path, StringBuilder sb){
        try (PrintWriter writer = new PrintWriter(new File(output_path))){
            writer.write(sb.toString());
        }catch (Exception e){
            System.err.println(e.getMessage());
        }
    }

    public static void test_5taxa(){
        String snet = "(((Q:1.5,R:1.5):2.5,L:4):1,(G:2,C:2):3);";
        HashMap<String, List<String>> taxonmap = new HashMap<>();

        Network net = Networks.readNetwork(snet);
        int i = 1;
        for(Object o: net.getLeaves()){
            NetNode node = (NetNode) o;
            System.out.println(node.getName());
            List<String> list = new ArrayList<>();
            list.add(String.valueOf(i));
            taxonmap.put(node.getName(), list);
            i += 1;
        }
        System.out.println(taxonmap);
        Enumerator enu = new Enumerator(5, taxonmap);
        StringBuilder sb = enu.computeLikelihood(net);

        String output_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/experiment/5taxa/gtprob.csv";
        enu.outputResult(output_path, sb);
    }

    public static void test_3taxa(){
        String snet = "((A:1,B:1):1,C:1);";
        HashMap<String, List<String>> taxonmap = new HashMap<>();

        Network net = Networks.readNetwork(snet);
        int i = 1;
        for(Object o: net.getLeaves()){
            NetNode node = (NetNode) o;
            System.out.println(node.getName());
            List<String> list = new ArrayList<>();
            list.add(String.valueOf(i));
            taxonmap.put(node.getName(), list);
            i += 1;
        }
        System.out.println(taxonmap);
        Enumerator enu = new Enumerator(3, taxonmap);
        StringBuilder sb = enu.computeLikelihood(net);

        String output_path = "/Users/zhen/Desktop/Zhen/research/phylogenetics/X2/experiment/3taxa_gtprob.csv";
        enu.outputResult(output_path, sb);
    }


//    public static void test(){
//        Network speciesTree = Networks.readNetwork("((A:1,B:1):1,C:1);");
//        List<String> taxa = new ArrayList<>();
//        Map<String,List<String>> species2alleles = new HashMap<>();
//        for (Object nd : speciesTree.getLeaves()) {
//            String species = ((NetNode)nd).getName();
//            taxa.add(species);
//            species2alleles.put(species, Arrays.asList(species));
//        }
//
//        List<Tree> allTopologies = Trees.generateAllBinaryTrees(taxa.toArray(new String[0]));
//        GeneTreeProbabilityYF calculator = new GeneTreeProbabilityYF();
//        calculator.preProcess(speciesTree, allTopologies, true);
//        double[] probArray = new double[allTopologies.size()];
//        calculator.calculateGTDistribution(speciesTree, allTopologies, species2alleles, probArray);
//        for (int i = 0; i < allTopologies.size(); i++) {
//            System.out.println("[Topology]" + allTopologies.get(i) + "[Probability]" + probArray[i]);
//        }
//    }


    public static void main(String[] args) {
        test_3taxa();
//        test();
    }
}
