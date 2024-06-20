package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq2.move.tree;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq2.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq2.move.Operator;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq2.start.UPGMATree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq2.start.distance.JCDistance;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq2.structs.UltrametricTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq2.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq2.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;

import java.util.*;


/*
 * @ClassName:   CoordindatedWideNNI
 * @Description: adapted from starbeast2/CoordinatedExchange.java pickwide
 * @Author:      Zhen Cao
 */


public class CoordindatedWideNNI extends TreeOperator {

    private double _logHR;

    private TNode[] treeNodes;

    public TNode aNode; // naming follows Rannala & Yang 2015
    public TNode bNode;
    public TNode yNode;
    public TNode cNode;
    public TNode zNode;

    private int nLeafNodes;
    private int nInternalNodes; // excludes the root node
    private int nNodes;

    public CoordindatedWideNNI(UltrametricTree tree) {
        super(tree);
        nLeafNodes = _tree.getTaxa().length;
        nInternalNodes = _tree.getInternalNodes().size();
        nNodes = _tree.getNodes().size();
    }

    @Override
    public double propose() {
        _violate = false;
        _logHR = Utils.INVALID_MOVE;
        treeNodes = _tree.getNodeArray();
        for (int i = 0; i < treeNodes.length; i++){
            STINode n = (STINode) treeNodes[i];
            n.setLabel(i);
        }

        if (!pickWide()) return _logHR;

        exchangeParents(cNode, bNode);
        _logHR = 0.0;

        _violate = true;
        return _logHR;
    }

    @Override
    public void undo() {
        if(_logHR == Utils.INVALID_MOVE) return;
        exchangeParents(cNode, bNode);
    }

    @Override
    public String getName() {
        return "CoordindatedWideNNI";
    }

    private boolean pickWide() {
        final int nNodesExceptRoot = nNodes - 1;
        final TNode rootNode = treeNodes[nNodesExceptRoot];

        // pick an internal node at random (excluding the root)
        final int yNodeNumber = nLeafNodes + Randomizer.getRandomInt(nInternalNodes - 1);
        yNode =  treeNodes[yNodeNumber];
        final double yNodeHeight = yNode.getNodeHeight();


        if (Randomizer.getRandomDouble() < 0.5) {
            aNode = getChild(yNode,0);
            bNode = getChild(yNode,1);
        } else {
            aNode = getChild(yNode,1);
            bNode = getChild(yNode,0);
        }

        // for all internal nodes (excluding the root)
        final STINode[] zNodes = new STINode[nNodesExceptRoot];

        czNodeFinder(yNode, rootNode, yNodeHeight, zNodes);

        // pick a cousin from the available candidates
        int cousinNodeNumber = Randomizer.getRandomInt(nNodesExceptRoot);
        zNode = zNodes[cousinNodeNumber];
        while (zNode == null) {
            cousinNodeNumber = Randomizer.getRandomInt(nNodesExceptRoot);
            //System.out.println(String.format("%d/%d", cousinNodeNumber, nNodesExceptRoot));
            zNode = zNodes[cousinNodeNumber];
        }

        cNode = treeNodes[cousinNodeNumber];

        return true;
    }

    private List<Integer> czNodeFinder(final TNode parentNode, final TNode currentNode, final double parentNodeHeight, final TNode[] zNodes) {
        // valid graft nodes (nodes defining branches which include the height of the parent node)
        final List<Integer> candidateList = new ArrayList<>();
        final double currentNodeHeight = currentNode.getNodeHeight();

        if (parentNode == currentNode) {
            return null;
        } else if (parentNodeHeight >= currentNodeHeight) {
            // this is a candidate node (would be a valid choice to graft parentNode to)
            candidateList.add(((STINode) currentNode).getLabel());
            return candidateList;
        } else {
            List<TNode> children = new ArrayList<TNode>();
            for (TNode n : currentNode.getChildren()) {
                children.add(n);
            }
            final List<Integer> leftCandidateNodeNumbers = czNodeFinder(parentNode, children.get(0), parentNodeHeight, zNodes);
            final List<Integer> rightCandidateNodeNumbers = czNodeFinder(parentNode, children.get(1), parentNodeHeight, zNodes);

            if (leftCandidateNodeNumbers == null) { // parent is the left child or descendant of the left child
                // therefore the current node is the most recent common ancestor connecting the parent and right candidates
                for (final Integer candidateNodeNumber: rightCandidateNodeNumbers) {
                    zNodes[candidateNodeNumber] = currentNode;
                }
                return null;
            } else if (rightCandidateNodeNumbers == null) { // parent is the right child or descendant of the right child
                // therefore the current node is the most recent common ancestor connecting the parent and left candidates
                for (final Integer candidateNodeNumber: leftCandidateNodeNumbers) {
                    zNodes[candidateNodeNumber] = currentNode;
                }
                return null;
            } else {
                candidateList.addAll(leftCandidateNodeNumbers);
                candidateList.addAll(rightCandidateNodeNumbers);
                return candidateList;
            }
        }
    }


    public static void main(String[] args) {
        Utils._SEED = 1;
        Map<String, String> locus = new HashMap<>();
        {
            locus.put("A", "CTTCGTGACGGGCTCGGCTCGTACGGTCAAGGGCACCTGAGCTAGGCAACTCAAGACGGGCGAGAGTCCCGCTACTACGGAAAAGGGTTTGCAGCTTCGAACATAGAACCACGGGACCTCCGGGACTCGCCCAATTGCGGTCGATGGCACTACGATTGGACAGGCGCTTACGCCAATTTACAGGTAGCGAGAGGTCTGCGCAGCAAAAGACTCCGCTTCACCCGGTGGTACCTGACTCGCGGCGCTGACCGTTCGCCGTAAGCGTTGACAATTTTCGGAAGCATCGCCAAAGTGCTGGAGTACGGGCTTCACACCCGCAACCACCGACAGCCCAACGAATCACTTCGCCGGTTGTCCTTCACCCCTGTTGGCAGAGGATGCTTGCGTGACTTTCATCCTTCTGTCCTTTCGTGCGGACTCGCACAGATCCTCCAAACAAGCGAGATCCGACCGATACTCTGCCCCCAGCAAGGCGGGGTTTCAGCGTCCCGTAGCTAGTAGGTATCCAGTCCAGAGGCACGTAGAGCACTCACCGTCCCCAGCCCCCCTCCACTCATCGCTGTGTTGAGGTCAAGTTCCCCGAGATCGTACTAACCGGTATGCAACCCACATCTGCCCATAGTCTACCCATCATACTACAGTCAGTAAACCCGGAGTGTATGGCTCGACTAATGTCCGTACAGCCAGCGCTCGATAAGTGCGTGCTCGAACGCTTACACCCCTGCTGGCACCAGGTAGCGCATGATCCACCAATGCAGCCCGAGCGCTTGCTGTGACTTCCGTCCGAAAGTATCGAGTTAAGCCTCAGCAATATGGGCCACCTTGACTGAGCTACACTCCCTCCCGTACGGGTGAAGTCCCCGCCGGCACAGGGGGGACACTTCTATTGAACATTTCGCCTCACCGAAGTCGCCCCCGAGTACGCCAGTGACGTACAGTCGCACGCGGGTTTGGTAAGGAGGAGCCAGATCGCGTCTACATGTTGAGTAGGTCCCGGCGC");
            locus.put("C", "CTACGTGACGGGCCCGGCTCGCACGGTCAATGGCACCTGCGCTAGGCAACTCGAGACGGGCGAGGGTCCCGCTACTACGGACAAGGGTGTGCAGCGTCGAACATAGAACCACGGGACCTCCGGTACTCGCTCAACCGCGGTCGGTGGTACTACGATTGGACGGGCACTTTCGGCAATTTACAGTTAGCGAGAGGTCTGCGCGGCAAAAGACTCCAATTCACCCGGCGTTACCTGACTCGCGGCACTGACCGTTCGCCGTAAGCGTTGACAATTTTCGGTAGCATCGCCAAAGTGCTGGAGTACGGGTTTCACACTTGTAGCCACCGACAGCCCAACGGATCAATTCGCCAGTTGCCCTTCACTCCTGTTAGCAGAGGATGCTAGCGTGAGTTTCATCCATCTGTCCTGTCGAGCGGACTCGCACAGATACTCCAAACAAGTGAGATCCGTCCGGTACTCTACCCCTCACAAAGCGGGATTTCAGCGTCTCGTAGTGAGTAGGCATCCGGTCCAGGGGCACGTAGAGCACTCACCGTCCTCAGCCCCCTTCTATTCATTGCTGTGTTGAAGCGAAGTTCCCCGAAATCGTACTAACCGGTATGCAACCCAGATCTGCCCATAGTTTCCACATCATACTACAGCCAGAAAACTCGGAGTGTATGGCTTGACTAATGTCCGTACTGCCAGCGCTAGGTAAGTACGTGCTCGAACGCTTACACCCCTGCTGGCACCTGGTTGCGCATGATCCACCAATGCAGCCCGAGCGCTTGCTGTGACTTCCGTCCAAAAGTATTGAATTAAGCCTCAGCAATATGGGCCACCCAGACTGAGCCACGCTCCCCTCCGTACGGATGAAGTCCCCACCGGGGCAGGGGGGACGCTTCTATTGAACACTTCGCCTCACCGAAGTCGCCCCCGAGTACGCCAGTGACGTACGGTCGCACGCGGGTTCGGTAAGGAGGGGTCAGATCGCGTCTGCATGTTGAGTAGGTCCCTGCGC");
            locus.put("G", "CTACGTGACGGGCCTGGCTCGCACGGTCAATGGCACTTGAGCTAGGCAACCCAAGACGGGCGAGGGTCCCGCTACTACGGAAAAGGGTGTGCAGCTTCGAACATAGAACCACGGGACCTCCGGGACTCGCTCAATTGCGGTCGGTGGTACTACGATTGGACGGGCACTTTCGGCAATTTACAATTAGCGAGAGGTCTGCGCAGCAAAAGACTCCAATTCACCCGGCGCTACCTGATTCGCGGCACTGACCGTTCGCCGTAAGCGTTGACAAATTTCGGTAGCATCGCCAAAGTGCTGGAGTACGGGCTTCACACTTGTAGCCACCGACAGCCCAACGAATCAATTCGCCAGTTGCCCTTCACTCCTGTTAGCAGAGGATGCCAGCATGAGATTCATCCATCTGTCCTTTCGAGCGGACTCGCACAGATACTCCAAACAAGCGAGATCCGACCGGTACTCTACCCCTCACAAAGCGGGATTTCAGCGTCTCGTAGTGAGTAGGTATCCGGTCCAGAGGCACGTAGAGCCCTCACTGTCCTCAGCCCCCCTCTACTCATTGCTGTGTTGAGGTGAAGTTCCCCGAAATCGTACTAACCGGTATGCAACCCAGATCTGCCCATAGGCTACACATCATACTACAGCCAGTAAACCCGGAGTGTATGGCTTGACTAAGGTCCGTACAGCCAGCGCTCGATAAGTACGTGCTCGAACGCTTACACCCCTGCTGGCACCAGGTAGCGTATGATCCACCGATGCAGCCCGAGCGCTTGCTGTGACTTCCGTCCGAAAGTATTGAGTTAAGCCTCAGCAATATGGGCCACCCTGACTGAGCTACGCTCCCTTCCGTACGGATGAAGTCCCCACCGGGACAGGGGGGACGCTTCTATTGAACACTTCACCTCACCGAAGTCGCCCCCGAGTACACCAGTGACGTACAGTCGCACGCGGGTTCGATAAGGAGGGGCCAGATCGCGTCTACATGTTGAGTAGGTCCCTGCGC");
//            locus.put("R", "TTGTGCGACGGGCGCGGCCCGGGCGATCAAGGGTACCAGAGTCAGGCAACTTAAGATGGGCGAGAGCCCCGCTAATACGGAAAGGAGTATGCAGCTCCGGACATGGGATCACTGGATCTCCAGGACTCGGTCGACTGCGGTCGGCGTTACTACAATCGGAGAGGCACATTCGCTAACTCATAGTTTGCGAGAGATCTGCGGAGCAAAGGATTCCCCTCTACTCGGCGCTACCTGACTCGCTACGCAGACCGTTCGCCGTAAGTGTTGCCAATCCCCGGAGGCATCGCCAAAGTACTGGAGTACGGGCTTAACACCTACAGCCACCGACAGCGCAACGAATCATTTCACCAGTTGCCTTTTACCCTTGTTAGCTGAGGATGCTAGCTTAACTTTCATCCTATGTTCCTCTCGGGCGGACTCTAACTGATCCCCCGAATGAGCGAGGTCTGACCGATACTCTGCCCCAGGCAAGGGGGGAGCCCATCCCCTTATAGTAAGCAAATCCCCAGTTCAGAGGCACATAGAGCACCCACCGTCCACAGCCCCTTTCCACTCACTGGTGCGCTGAGGTGAAAGTGCCCGAAATCCTACGAATTGGTATGCAACCCAGATCTGTAGGCAGGCTACATGTCATACTACAGCTCGTAAACTTGGAGTGTATGGCTAGACTGATATCCGAACAAACAACGCTCGACAAGCGCGACCTCGACCGCTCACACCCTTGCTGACACCAAGCAACACATGATCCATCAGTGCAGCCCCAACGTTTTTTGTGACCTCCGTCCGAAAGTATGGATTTGAGCCTCAGCAATGTGGCCCACCATGGCCGAGCTACGCTCCCCTACGTACGGATGATTTCCCCGCCGGGACAGGCGGGACGGTTCTATTAAACATCTCACCTTACTGATGTCGCCCCCGGGTACGGCAGCGACGTACAGCCGCACGCGAGCTTGGTAAGGAGGAGCCAGATCGTGCCTACATGTTGAGTAGGTCCCTACTC");
//            locus.put("Q", "CCATACGATGGGCTCGGCTCGTATGATTGAGGGCACCGAAGCTAGGCGACTCAAGATGGGCGAGGGCTCCGCGAATACGGAAAAGGGTATGCAGCTTCGGCCATAGGACCACGTGATCTCCGGGACTCGCCCAATTGAGATCGGCGTTACTACAATTAGACGAGCACATTCGCCAATTTATAGTTAACAAGAGATCTGCGTAGCAAAAGATTCAACGTTACCCGGCGCTACCAGACCCGCGGCACAGGCCGTTTGCCCTAAGCGTTGACAATCTTCGGATTCATCGCCAAAGTGCTGGAGTACAGGCTTCACACCTACAGCCACTGACAGCCTAACGAATCACTTCACCAATTGCCTTTCACCCCTGTTAGCGGAGGGCGCTAGCATAACTTTCGTCCTACGTTCCTCTCGTACGGATTCGGACAGATCCTCCGAGCAAGCGAAGTCCGGCCGATACTCTGCCCCTAGCAAGGCGGGATCCCGTCGCCTTGTAGTGAGCAGATATCCAGTTTGGGAGCACATAGAGCACCGACCGTCCACAGTCCCCTTCTCTTCATTGGTGCGTTGAGGTGAAAGTTCCCGAAATCCTACAAACTGGTATGCAAACCGGATCTGCAGGTTGGCTACGTATCATACTACAGCCCGTAAACTCGGAGTGCATGGCTTGACTAACATCCGTACAAACAGCGCTCGATGAGTGCGACCTCGGCAGCTTACACCCTTGCTGACACCAAATAGCGCATGATCCACCAGTACAGCCCAAGCGCCTTTCGCGTCCTCCGCCCGAATGTATGGATGTAAGCCTGAGTCACGTGGACCACCGTGCCCGAGCTACGCTCCCTTACGTGCGGATGATGTCCCCGCCGGGACAGGCGGAACGCTTCTATTGAACATTTCACCTCGCTGAAGTCGCCCCTGAGTACGCCAGCAACGTACAGTCGCACGCGAGCTCAGTAAGGAGGAGCCATATCGCGTCTACATGATGCGTAGGTCCCTGCGC");
//            locus.put("L", "CTACGCGACGGGCTCGGCTCGTACAATAAAGGGCACCAGAGCTAGGTAACCCAAGATAGGCGAGGGCTTCGCTAATACGGAGAAGGGTATTCAGCTTCGGTCACAGGGCCACGTGATCTCCGGGACTCGCCCAGTTGCGATCGGCGTTACTACAATTGGAGGGGCACATTCGCCAACTTACAGTTAGTGAGGGATCTGCGTATCAAAAAACTCCACACTACTCGGCGCGACCTGGCTCGTGGCGCAGGCCGTTTGCCGTAAGTGTTGAGAATTTTCGGAAGAATGGCCGGAGTGTTAGAGTCCAGGCTTCACACCTACAGCCATTGACAGCTCAACGAGTCATTTCACCAGTTGCCTTTCACCCCTGTTAGCAGAGGGAGCTAGCGCAGCTTTCATCCTATGTCCATCTCGGGCGGACTCGGACAGATCCTCCGAACAAGCGAGGTCCGGTCGATACTCTGCCCCTAACAAGGCAGGATCCCATCGCCTTGTAGTGAGCAGATATCCAGTTTAGGGGCGCATAGAGCACCCACCGTCCACAGCCCCCTTGTACTTATTGGTGTGTTGAGGTGAAAGTCCCCGAAATCCCAGGAACTTGTATGCAACCCAGATCTGCAGGTAGGCTACGTATCATACTGTACCCCGCAAACTCGGAGTGTGTGACTGGACTAACGTCCGTACAAACAGCGCTCGATAGGTGCGACCTCGACAGCTTACACCCTTGCTGGCACCAAATAGCGCGTGATCCACCAGTGCAGACCAAACGTTTTCTGCGCCCTCCGTCCGAAAGTACGGGTTTAAGTCTCAGTAACATGGTCCACTATGACCGAACTATACCCCCTTACGTACGGGTGATGTCCCCGCCGGGACGGGCGGGACGCTTCTACTGAACGCTTCACTTCACTGAAGTCGCCCCTGAGTATGCCAGCAACGTACAGTCGCACGCGAGCTCGGTAAAGAGGAACCAGATCCCGTCGACATGTTGAGTTGGTCCCTATGC");
        }
        Alignment seq = new Alignment(locus);
        UPGMATree upgma = new UPGMATree(new JCDistance( seq.getAlignment() ));
        UltrametricTree template = upgma.getUltrametricTree();
        System.out.println(template.getTree().toNewick());
        {
            UltrametricTree tree = new UltrametricTree(template);
            System.out.println(tree.getTree().toNewick());

            int runs = 10000;
            int counter = 0;
            int test = 0;

            for(int i = 0; i < runs; i++) {
                Operator op = new CoordindatedWideNNI(tree);
                double logHR = op.propose();

                if(logHR != Utils.INVALID_MOVE) counter++;

                else{
                    System.out.println(tree.toString());
                }
                if(tree.checkUltrametric()) test++;
            }
            System.out.println(test == runs);
            System.out.println(tree.getTree().toNewick());

            System.out.printf("%d out of %d\n", counter, runs);
        }
        {
            UltrametricTree tree = new UltrametricTree(template);
            System.out.println(tree.getTree().toNewick());

            int runs = 10000;
            int counter = 0;
            int test = 0;
            int diff_topo = 0;
            for(int i = 0; i < runs; i++) {
                Operator op = new CoordindatedWideNNI(tree);
                double logHR = op.propose();
                if(logHR != Utils.INVALID_MOVE) counter++;
                if(tree.compareTo(template) != 0) {
                    diff_topo++;
                }
                op.undo();
                if(tree.compareTo(template) == 0) {
                    test++;
                } else {
                    System.out.println(tree.getTree().toNewick());
                }
            }
            System.out.println(test == runs);
            System.out.printf("different topologies %d out of %d\n", diff_topo, runs);
            System.out.printf("%d out of %d\n", counter, runs);
            System.out.println(tree.getTree().toNewick());
        }
    }

}
