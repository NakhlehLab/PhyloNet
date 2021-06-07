package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.LocalGenealogy;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.TreeNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.ParseException;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.*;
import java.util.*;

/**
 * Parse an INDELible control file to build sequence of true genealogies.
 * Only applicable to HCG case with four genealogy classes. Can be extended to general case.
 *
 * Created by Xinhao Liu on 7/8/20.
 */
public class ControlFileParser {
    private static void buildGTNodeHeight(STITree<TreeNodeInfo> gt) {
        for (TNode node:gt.postTraverse()) {
            if (node.isLeaf()) {
                node.setNodeHeight(0);
            } else {
                double height = 0;
                for (TNode child:node.getChildren()) {
                    height = Math.max(height, child.getNodeHeight() + child.getParentDistance());
                }
                node.setNodeHeight(height);
            }
        }
    }

    public static String treeString2GenealogyClass(String newickString) throws IOException, ParseException {
        STITree<TreeNodeInfo> gt = new STITree<>(newickString);
        Trees.convertToLexicographicTree(gt);
        buildGTNodeHeight(gt);

        try {
            if (Trees.haveSameRootedTopology(gt, new STITree<>("((1,3), 2);"))) {
                return "HG";
            } else if (Trees.haveSameRootedTopology(gt, new STITree<>("((2,3), 1);"))) {
                return "CG";
            } else {
                if (Trees.getInternalNodes(gt).get(0).getNodeHeight() >= 5.5) {
                    return "HC2";
                } else if (Trees.getInternalNodes(gt).get(0).getNodeHeight() < 5.5) {
                    return "HC1";
                }
//                if (Trees.getInternalNodes(gt).get(0).getNodeHeight() >= 5.5) {
//                    return "HC";
//                } else if (Trees.getInternalNodes(gt).get(0).getNodeHeight() < 5.5) {
//                    return "HC";
//                }
                System.out.println("ERROR IN ControlFileParser.treeString2GenealogyClass!!!!!!!!");
                System.exit(1);
            }
        } catch (IOException | ParseException e) {
            e.printStackTrace();
        }
        return null;
    }

    public static List<String> parse(String filePath) throws IOException, ParseException {
        File file = new File(filePath);
        BufferedReader br = new BufferedReader(new FileReader(file));
        Map<String, String> treeName2GenealogyClass = new HashMap<>();
        List<String> trueGenealogySequence = new ArrayList<>();

        String line;
        while ((line = br.readLine()) != null) {
            if (line.startsWith("[TREE]")) {
                String[] splited = line.split("\\s+");
                String treeName = splited[1];
                String treeString = splited[2];
                treeName2GenealogyClass.put(treeName, treeString2GenealogyClass(treeString));
            }
            if (line.startsWith("[PARTITIONS]")) {
                String[] splited = line.split("\\s+");
                String treeName = splited[2].substring(1);
                int segmentLength = Integer.parseInt(splited[4].substring(0, splited[4].length() - 1));
                for (int i = 0; i < segmentLength; i++) {
                    trueGenealogySequence.add(treeName2GenealogyClass.get(treeName));
                }
            }
        }
        return trueGenealogySequence;
    }


    public static void main(String[] args) throws IOException, ParseException {
        List<String> sequence = parse("/Users/xinhaoliu/Desktop/Research/Scripts/HMMtestdata/TrueHCG_500000/2/control.txt");
        System.out.println(sequence.size());
        for (String genealogy:sequence) {
            System.out.println(genealogy);
        }
        //System.out.println(treeString2GenealogyClass("(3:5.592,(1:4.296,2:4.296):1.296);"));
    }
}
