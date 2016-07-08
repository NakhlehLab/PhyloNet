package edu.rice.cs.bioinfo.programs.phylonet.algos.clustering;

import edu.rice.cs.bioinfo.library.programming.MutableTuple;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.GLASSInference;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.MDCInference_Rooted;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.MajorityConsensusInference;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.Solution;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.InferNetworkMLFromGTT_SingleTreePerLocus;
import edu.rice.cs.bioinfo.programs.phylonet.algos.network.NetworkPseudoLikelihoodFromGTT_SingleTreePerLocus;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 7/7/16
 * Time: 2:07 PM
 * To change this template use File | Settings | File Templates.
 */
public class InferTreeWrapper {
    public void summarizeData(List originalGTs, Map<String,String> allele2species, List distinctGTs){
        Map<String, MutableTuple<Tree,Double>> exp2tree = new HashMap<String, MutableTuple<Tree, Double>>();
        for(Object list: originalGTs) {
            for (MutableTuple<Tree, Double> gtTuple : (List<MutableTuple<Tree, Double>>)list) {
                for (TNode node : gtTuple.Item1.getNodes()) {
                    node.setParentDistance(TNode.NO_DISTANCE);
                }
                String exp = Trees.getLexicographicNewickString(gtTuple.Item1, allele2species);
                MutableTuple<Tree, Double> existingTuple = exp2tree.get(exp);
                if (existingTuple == null) {
                    existingTuple = gtTuple;
                    exp2tree.put(exp, existingTuple);

                } else {
                    existingTuple.Item2 += gtTuple.Item2;
                }

            }
        }
        distinctGTs.addAll(exp2tree.values());
    }

    protected Tree getMajorityConsensusTree(List originalGTs, Map<String,List<String>> species2alleles){

        List dataForStartingNetwork = new ArrayList();
        summarizeData(originalGTs, null, dataForStartingNetwork);

        MajorityConsensusInference majorityConsensusInference = new MajorityConsensusInference();
        Tree MajorityConsensusTree;
        MajorityConsensusTree = Trees.generateRandomBinaryResolution(majorityConsensusInference.inferSpeciesTree(dataForStartingNetwork, true, 0));
        return MajorityConsensusTree;
    }

    protected Tree getMDCTree(List originalGTs, Map<String,List<String>> species2alleles){

        List dataForStartingNetwork = new ArrayList();
        summarizeData(originalGTs, null, dataForStartingNetwork);

        Map<String,String> allele2species = null;
        if(species2alleles!=null){
            allele2species = new HashMap<String, String>();
            for(Map.Entry<String,List<String>> entry: species2alleles.entrySet()){
                String species = entry.getKey();
                for(String allele: entry.getValue()){
                    allele2species.put(allele,species);
                }
            }
        }

        MDCInference_Rooted mdc = new MDCInference_Rooted();
        Solution sol;
        if(allele2species==null){
            sol = (Solution)(mdc.inferSpeciesTree(dataForStartingNetwork, false, 1, false, true, -1).get(0));
        }
        else{
            sol = (Solution)(mdc.inferSpeciesTree(dataForStartingNetwork, allele2species, false, 1, false, true, -1).get(0));
        }
        Tree MDCTree= Trees.generateRandomBinaryResolution(sol._st);
        return MDCTree;
    }

    protected Tree getDemocraticVoteTree(List originalGTs, Map<String,List<String>> species2alleles){


        List dataForStartingNetwork = new ArrayList();
        summarizeData(originalGTs, null, dataForStartingNetwork);

        List<MutableTuple<Tree,Double>> gts = dataForStartingNetwork;
        double maxCount = 0;
        Tree maxTree = (Tree)(gts.get(0).Item1);

        for(MutableTuple tuple : gts){
            if(maxCount < (double)(tuple.Item2)){
                maxTree = (Tree)(tuple.Item1);
                maxCount = (double)(tuple.Item2);
            }
        }
        return maxTree;
    }

    protected Tree getMLTree(List<List<MutableTuple<Tree,Double>>> originalGTs, Map<String,List<String>> species2alleles, Tree startTree){
        if(startTree == null) {
            String[] leaves = originalGTs.get(0).get(0).Item1.getLeaves();
            int taxa = originalGTs.get(0).get(0).Item1.getLeaves().length;
            startTree = Trees.generateRandomTree(Arrays.copyOfRange(leaves, 0, taxa));
        }

        InferNetworkMLFromGTT_SingleTreePerLocus calculator = new InferNetworkMLFromGTT_SingleTreePerLocus();

        LinkedList<Tuple<Network, Double>> resultList = new LinkedList<Tuple<Network, Double>>();

        //calculator.setParallel(8);
        //calculator.setNumRuns(1);
        calculator.setStartNetwork(Networks.readNetwork(startTree.toNewick()));
        calculator.inferNetwork(originalGTs, species2alleles, 0, 3, false, resultList);

        return Trees.readTree(resultList.get(0).Item1.toString());
    }

    protected Tree getGLASSTree(List<Tree> gtList, Map<String,List<String>> species2alleles){
        GLASSInference glassInference = new GLASSInference();
        Tree inferredST = glassInference.inferSpeciesTree(gtList);
        return inferredST;
    }

    protected String getASTRALTree(List<Tree> gtList) {
        String path = "/Users/zhujiafan/Documents/Luay/ASTRAL/astral.4.10.2.jar";
        String result = "";
        int n = gtList.size();
        try {
            PrintWriter out = new PrintWriter("intree_ASTRAL");
            for(int i = 0 ; i < n ; i++)
                out.println(gtList.get(i).toString());
            out.close();
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }

        try {
            File f = new File("outtree_ASTRAL");
            f.delete();
        } catch (Exception e) {
        }

        try {
            Process proc = Runtime.getRuntime().exec("java -jar " + path + " -i intree_ASTRAL -o outtree_ASTRAL", null, null);
            proc.waitFor();
            File file = new File("outtree_ASTRAL");
            Scanner scanner = new Scanner(file);

            result = scanner.nextLine();

            scanner.close();

        } catch (Exception e) {
            e.printStackTrace();
        }
        result = result.replaceAll("1:", ":");
        return result;

    }

    public String inferTreeByMethod(List<Tree> gtList, Map<String,List<String>> species2alleles, String method){

        //System.out.println("Method: " + method);

        List<List<MutableTuple<Tree, Double>>> originalGTs = new ArrayList<>();

        for(Tree tree : gtList){
            originalGTs.add(Arrays.asList(new MutableTuple(Trees.readTree(tree.toNewick()), 1.0)));
        }


        Network<Object> startingNetwork = null;

        if(startingNetwork == null){


            Tree startingTree = null;
            if(method.equals("MDC"))
                startingTree = getMDCTree(originalGTs, species2alleles);
            else if(method.equals("MajorityConsensus"))
                startingTree = getMajorityConsensusTree(originalGTs, species2alleles);
            else if(method.equals("DemocraticVote"))
                startingTree = getDemocraticVoteTree(originalGTs, species2alleles);
            else if(method.equals("GLASS"))
                startingTree = getGLASSTree(gtList, species2alleles);
            else if(method.equals("ML_Random"))
                startingTree = getMLTree(originalGTs, species2alleles, null);
            else if(method.equals("ML_MDC"))
                startingTree = getMLTree(originalGTs, species2alleles, getMDCTree(originalGTs, species2alleles));
            else if(method.equals("ML_MajorityConsensus"))
                startingTree = getMLTree(originalGTs, species2alleles, getMajorityConsensusTree(originalGTs, species2alleles));
            else if(method.equals("ML_DemocraticVote"))
                startingTree = getMLTree(originalGTs, species2alleles, getDemocraticVoteTree(originalGTs, species2alleles));
            else if(method.equals("ML_GLASS"))
                startingTree = getMLTree(originalGTs, species2alleles, getGLASSTree(gtList, species2alleles));
            else if(method.equals("ASTRAL")) {
                return getASTRALTree(gtList);
            }
            else
                throw new RuntimeException("Unknown inferring method: " + method);

            startingNetwork = edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks.readNetwork(startingTree.toString());

        }


        for(NetNode<Object> node: startingNetwork.dfs()){
            for(NetNode<Object> parent: node.getParents()){
                if(node.getParentDistance(parent) == NetNode.NO_DISTANCE)
                    node.setParentDistance(parent,1.0);
                if(node.isNetworkNode()){
                    if(node.getParentProbability(parent) == NetNode.NO_PROBABILITY)
                        node.setParentProbability(parent, 0.5);
                }
            }
        }

        //System.out.println(startingNetwork.toString());
        return startingNetwork.toString();
    }
}
