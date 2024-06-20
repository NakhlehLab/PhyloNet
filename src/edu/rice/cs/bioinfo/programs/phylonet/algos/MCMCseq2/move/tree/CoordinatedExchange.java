package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq2.move.tree;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq2.felsenstein.alignment.Alignment;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq2.move.Operator;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq2.start.UPGMATree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq2.start.distance.JCDistance;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq2.structs.UltrametricTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq2.util.Randomizer;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCseq2.util.Utils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import com.google.common.collect.SetMultimap;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.util.*;

/*
 * @ClassName:   CoordindatedWideNNI
 * @Description: adapted from starbeast2/CoordinatedExchange.java pickwide
 * @Author:      Zhen Cao
 */


public class CoordinatedExchange extends TreeOperator {

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
    private int czBranchCount;

    private List<List<SortedMap<TNode, TNode>>> movedNodes;
    private List<SetMultimap<Integer, TNode>> graftNodes;

    public CoordinatedExchange(UltrametricTree tree) {
        super(tree);
        nLeafNodes = _tree.getTaxa().length;
        nInternalNodes = _tree.getInternalNodes().size();
        nNodes = _tree.getNodes().size();

    }

    @Override
    public double propose() {
        _logHR = Utils.INVALID_MOVE;
        treeNodes = _tree.getNodeArray();
        for (int i = 0; i < treeNodes.length; i++){
            STINode n = (STINode) treeNodes[i];
            n.setLabel(i);
        }

        if (!pickWide()) return _logHR;


//        fillNodes();
        pruneAndRegraft(yNode, cNode, bNode);

        _violate  = true;
        _logHR = 0.0;
        return _logHR;
    }



    // fills forward nodes by destination branch (c through z)
    public void fillNodes() {
        // must be done before any changes are made to the gene or species trees
        TNode childNode = cNode;
        movedNodes = new ArrayList<>();
        graftNodes = new ArrayList<>();
        czBranchCount = 0;
        while (childNode != zNode) {
            czBranchCount++;
            final TNode parentNode = childNode.getParent();
            childNode = parentNode;
        }
    }

    public UltrametricTree getTree(){
        return _tree;
    }

    private void pruneAndRegraft(final TNode nodeToMove, final TNode newChild, final TNode disownedChild) {
        final TNode sourceParent = nodeToMove.getParent();
        final TNode destinationParent = newChild.getParent();

        ((STINode) newChild).setParent2((STINode) nodeToMove);
        ((STINode) disownedChild).setParent2((STINode) sourceParent);
        ((STINode) nodeToMove).setParent2((STINode) destinationParent);

        newChild.setParentDistance(_tree.getNodeHeight(nodeToMove)-_tree.getNodeHeight(newChild));
        disownedChild.setParentDistance(_tree.getNodeHeight(sourceParent) - _tree.getNodeHeight(disownedChild));
        nodeToMove.setParentDistance(_tree.getNodeHeight(destinationParent) - _tree.getNodeHeight(nodeToMove));

    }


    @Override
    public void undo() {
        if(_logHR == Utils.INVALID_MOVE) return;
        final TNode bNodeTmp = bNode;
        final TNode cNodeTmp = cNode;

        bNode = cNodeTmp;
        cNode = bNodeTmp;

        fillNodes();
//        System.out.println(czBranchCount);
        pruneAndRegraft(yNode, cNode, bNode);

    }

    @Override
    public String getName() {
        return "CoordindatedExchange";
    }

    private boolean pickWide() {
        final int nNodesExceptRoot = nNodes - 1;
        final TNode rootNode = treeNodes[nNodesExceptRoot];

        // pick an internal node at random (excluding the root)
        final int yNodeNumber = nLeafNodes + Randomizer.getRandomInt(nInternalNodes - 1);
        yNode = treeNodes[yNodeNumber];
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
//            List<TNode> children = new ArrayList<TNode>();
//            for (TNode n : currentNode.getChildren()) {
//                children.add(n);
//            }
            final List<Integer> leftCandidateNodeNumbers = czNodeFinder(parentNode, getChild(currentNode,0), parentNodeHeight, zNodes);
            final List<Integer> rightCandidateNodeNumbers = czNodeFinder(parentNode, getChild(currentNode,1), parentNodeHeight, zNodes);

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

    public static void testUltrametric() {

        Utils._SEED = 2;
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
        UltrametricTree tree = new UltrametricTree(template);
        System.out.println(tree.getTree().toNewick());

        int runs = 10000;
        int counter = 0;
        int test = 0;

        for(int i = 0; i < runs; i++) {
            Operator op = new CoordinatedExchange(tree);
            double logHR = op.propose();

            if(logHR != Utils.INVALID_MOVE) counter++;
            if(tree.checkUltrametric()) test++;

        }
        System.out.println(test == runs);
        System.out.println(tree.getTree().toNewick());

        System.out.printf("%d out of %d\n", counter, runs);
    }

    public static void testSPR_big(){
        Utils._SEED = 2;
        for (int i = 0; i < 10; i++){
            System.out.println(Randomizer.getRandomDouble());
        }

        String newickSpeciesTree = "((((s1:0.33109175037666511,s3:0.33109175037666511):0.19728320951827943,(s4:0.35745288551058663,s5:0.35745288551058663):0.17092207438435791):0.13952939110386009,(s0:0.23635691859481153,s2:0.23635691859481153):0.43154743240399307):0.24425886940385155,s6:0.91216322040265618)";
        boolean bIsParent = true;
        String bTipLabel = "s5";
        String cTipLabel = null;
        boolean cIsParent = false;

        Tree st = null;
        try{
            st = new STITree<>(newickSpeciesTree);
        }catch (Exception e){
            e.printStackTrace();
        }

        UltrametricTree speciesTree = new UltrametricTree(st);
        TNode cNode = null;
        TNode bNode = null;
        for (String id: speciesTree.getTree().getLeaves()) {
            TNode n = st.getNode(id);
            if (id.equals(bTipLabel)) {
                if (bIsParent) bNode = n.getParent();
                else bNode = n;
            } else if (id.equals(cTipLabel)) {
                if (cIsParent) cNode = n.getParent();
                else cNode = n;
            }
        }

        CoordinatedExchange coex = new CoordinatedExchange(speciesTree);

        TNode yNode = bNode.getParent();
        TNode zNode = yNode.getParent();
        TNode aNode = (bNode == coex.getChild(yNode, 1)) ? coex.getChild(yNode, 0) : coex.getChild(yNode, 1);

        if (cNode == null) {
            cNode = (yNode == coex.getChild(zNode, 1)) ? coex.getChild(zNode, 0) : coex.getChild(zNode, 1);
        }

        coex.aNode = aNode;
        coex.bNode = bNode;
        coex.cNode = cNode;
        coex.yNode = yNode;
        coex.zNode = zNode;
        System.out.println(coex.getTree().toString());
        final double calculatedLogHR = coex.propose();
        System.out.println(coex.getTree().toString());
        System.out.println(calculatedLogHR);
        coex.undo();
        System.out.println(coex.getTree().toString());
    }

    public static void testSPR_small(){
        Utils._SEED = 2;
        for (int i = 0; i < 10; i++){
            System.out.println(Randomizer.getRandomDouble());
        }

        String newickSpeciesTree = "((A:1.0,B:1.0):2.0,C:3.0)";;
        boolean bIsParent = false;
        String bTipLabel = "B";
        String cTipLabel = null;
        boolean cIsParent = false;

        Tree st = null;
        try{
            st = new STITree<>(newickSpeciesTree);
        }catch (Exception e){
            e.printStackTrace();
        }

        UltrametricTree speciesTree = new UltrametricTree(st);
        TNode cNode = null;
        TNode bNode = null;
        for (String id: speciesTree.getTree().getLeaves()) {
            TNode n = st.getNode(id);
            if (id.equals(bTipLabel)) {
                if (bIsParent) bNode = n.getParent();
                else bNode = n;
            } else if (id.equals(cTipLabel)) {
                if (cIsParent) cNode = n.getParent();
                else cNode = n;
            }
        }

        CoordinatedExchange coex = new CoordinatedExchange(speciesTree);

        TNode yNode = bNode.getParent();
        TNode zNode = yNode.getParent();
        TNode aNode = (bNode == coex.getChild(yNode, 1)) ? coex.getChild(yNode, 0) : coex.getChild(yNode, 1);

        if (cNode == null) {
            cNode = (yNode == coex.getChild(zNode, 1)) ? coex.getChild(zNode, 0) : coex.getChild(zNode, 1);
        }

        coex.aNode = aNode;
        coex.bNode = bNode;
        coex.cNode = cNode;
        coex.yNode = yNode;
        coex.zNode = zNode;
        System.out.println(coex.getTree().toString());
        final double calculatedLogHR = coex.propose();
        System.out.println(coex.getTree().toString());
        System.out.println(calculatedLogHR);
        coex.undo();
        System.out.println(coex.getTree().toString());
    }

    public static void testUndo(){
        Map<String, String> locus = new HashMap<>();
        {
            locus.put("A", "CTTCGTGACGGGCTCGGCTCGTACGGTCAAGGGCACCTGAGCTAGGCAACTCAAGACGGGCGAGAGTCCCGCTACTACGGAAAAGGGTTTGCAGCTTCGAACATAGAACCACGGGACCTCCGGGACTCGCCCAATTGCGGTCGATGGCACTACGATTGGACAGGCGCTTACGCCAATTTACAGGTAGCGAGAGGTCTGCGCAGCAAAAGACTCCGCTTCACCCGGTGGTACCTGACTCGCGGCGCTGACCGTTCGCCGTAAGCGTTGACAATTTTCGGAAGCATCGCCAAAGTGCTGGAGTACGGGCTTCACACCCGCAACCACCGACAGCCCAACGAATCACTTCGCCGGTTGTCCTTCACCCCTGTTGGCAGAGGATGCTTGCGTGACTTTCATCCTTCTGTCCTTTCGTGCGGACTCGCACAGATCCTCCAAACAAGCGAGATCCGACCGATACTCTGCCCCCAGCAAGGCGGGGTTTCAGCGTCCCGTAGCTAGTAGGTATCCAGTCCAGAGGCACGTAGAGCACTCACCGTCCCCAGCCCCCCTCCACTCATCGCTGTGTTGAGGTCAAGTTCCCCGAGATCGTACTAACCGGTATGCAACCCACATCTGCCCATAGTCTACCCATCATACTACAGTCAGTAAACCCGGAGTGTATGGCTCGACTAATGTCCGTACAGCCAGCGCTCGATAAGTGCGTGCTCGAACGCTTACACCCCTGCTGGCACCAGGTAGCGCATGATCCACCAATGCAGCCCGAGCGCTTGCTGTGACTTCCGTCCGAAAGTATCGAGTTAAGCCTCAGCAATATGGGCCACCTTGACTGAGCTACACTCCCTCCCGTACGGGTGAAGTCCCCGCCGGCACAGGGGGGACACTTCTATTGAACATTTCGCCTCACCGAAGTCGCCCCCGAGTACGCCAGTGACGTACAGTCGCACGCGGGTTTGGTAAGGAGGAGCCAGATCGCGTCTACATGTTGAGTAGGTCCCGGCGC");
            locus.put("C", "CTACGTGACGGGCCCGGCTCGCACGGTCAATGGCACCTGCGCTAGGCAACTCGAGACGGGCGAGGGTCCCGCTACTACGGACAAGGGTGTGCAGCGTCGAACATAGAACCACGGGACCTCCGGTACTCGCTCAACCGCGGTCGGTGGTACTACGATTGGACGGGCACTTTCGGCAATTTACAGTTAGCGAGAGGTCTGCGCGGCAAAAGACTCCAATTCACCCGGCGTTACCTGACTCGCGGCACTGACCGTTCGCCGTAAGCGTTGACAATTTTCGGTAGCATCGCCAAAGTGCTGGAGTACGGGTTTCACACTTGTAGCCACCGACAGCCCAACGGATCAATTCGCCAGTTGCCCTTCACTCCTGTTAGCAGAGGATGCTAGCGTGAGTTTCATCCATCTGTCCTGTCGAGCGGACTCGCACAGATACTCCAAACAAGTGAGATCCGTCCGGTACTCTACCCCTCACAAAGCGGGATTTCAGCGTCTCGTAGTGAGTAGGCATCCGGTCCAGGGGCACGTAGAGCACTCACCGTCCTCAGCCCCCTTCTATTCATTGCTGTGTTGAAGCGAAGTTCCCCGAAATCGTACTAACCGGTATGCAACCCAGATCTGCCCATAGTTTCCACATCATACTACAGCCAGAAAACTCGGAGTGTATGGCTTGACTAATGTCCGTACTGCCAGCGCTAGGTAAGTACGTGCTCGAACGCTTACACCCCTGCTGGCACCTGGTTGCGCATGATCCACCAATGCAGCCCGAGCGCTTGCTGTGACTTCCGTCCAAAAGTATTGAATTAAGCCTCAGCAATATGGGCCACCCAGACTGAGCCACGCTCCCCTCCGTACGGATGAAGTCCCCACCGGGGCAGGGGGGACGCTTCTATTGAACACTTCGCCTCACCGAAGTCGCCCCCGAGTACGCCAGTGACGTACGGTCGCACGCGGGTTCGGTAAGGAGGGGTCAGATCGCGTCTGCATGTTGAGTAGGTCCCTGCGC");
            locus.put("G", "CTACGTGACGGGCCTGGCTCGCACGGTCAATGGCACTTGAGCTAGGCAACCCAAGACGGGCGAGGGTCCCGCTACTACGGAAAAGGGTGTGCAGCTTCGAACATAGAACCACGGGACCTCCGGGACTCGCTCAATTGCGGTCGGTGGTACTACGATTGGACGGGCACTTTCGGCAATTTACAATTAGCGAGAGGTCTGCGCAGCAAAAGACTCCAATTCACCCGGCGCTACCTGATTCGCGGCACTGACCGTTCGCCGTAAGCGTTGACAAATTTCGGTAGCATCGCCAAAGTGCTGGAGTACGGGCTTCACACTTGTAGCCACCGACAGCCCAACGAATCAATTCGCCAGTTGCCCTTCACTCCTGTTAGCAGAGGATGCCAGCATGAGATTCATCCATCTGTCCTTTCGAGCGGACTCGCACAGATACTCCAAACAAGCGAGATCCGACCGGTACTCTACCCCTCACAAAGCGGGATTTCAGCGTCTCGTAGTGAGTAGGTATCCGGTCCAGAGGCACGTAGAGCCCTCACTGTCCTCAGCCCCCCTCTACTCATTGCTGTGTTGAGGTGAAGTTCCCCGAAATCGTACTAACCGGTATGCAACCCAGATCTGCCCATAGGCTACACATCATACTACAGCCAGTAAACCCGGAGTGTATGGCTTGACTAAGGTCCGTACAGCCAGCGCTCGATAAGTACGTGCTCGAACGCTTACACCCCTGCTGGCACCAGGTAGCGTATGATCCACCGATGCAGCCCGAGCGCTTGCTGTGACTTCCGTCCGAAAGTATTGAGTTAAGCCTCAGCAATATGGGCCACCCTGACTGAGCTACGCTCCCTTCCGTACGGATGAAGTCCCCACCGGGACAGGGGGGACGCTTCTATTGAACACTTCACCTCACCGAAGTCGCCCCCGAGTACACCAGTGACGTACAGTCGCACGCGGGTTCGATAAGGAGGGGCCAGATCGCGTCTACATGTTGAGTAGGTCCCTGCGC");
            locus.put("R", "TTGTGCGACGGGCGCGGCCCGGGCGATCAAGGGTACCAGAGTCAGGCAACTTAAGATGGGCGAGAGCCCCGCTAATACGGAAAGGAGTATGCAGCTCCGGACATGGGATCACTGGATCTCCAGGACTCGGTCGACTGCGGTCGGCGTTACTACAATCGGAGAGGCACATTCGCTAACTCATAGTTTGCGAGAGATCTGCGGAGCAAAGGATTCCCCTCTACTCGGCGCTACCTGACTCGCTACGCAGACCGTTCGCCGTAAGTGTTGCCAATCCCCGGAGGCATCGCCAAAGTACTGGAGTACGGGCTTAACACCTACAGCCACCGACAGCGCAACGAATCATTTCACCAGTTGCCTTTTACCCTTGTTAGCTGAGGATGCTAGCTTAACTTTCATCCTATGTTCCTCTCGGGCGGACTCTAACTGATCCCCCGAATGAGCGAGGTCTGACCGATACTCTGCCCCAGGCAAGGGGGGAGCCCATCCCCTTATAGTAAGCAAATCCCCAGTTCAGAGGCACATAGAGCACCCACCGTCCACAGCCCCTTTCCACTCACTGGTGCGCTGAGGTGAAAGTGCCCGAAATCCTACGAATTGGTATGCAACCCAGATCTGTAGGCAGGCTACATGTCATACTACAGCTCGTAAACTTGGAGTGTATGGCTAGACTGATATCCGAACAAACAACGCTCGACAAGCGCGACCTCGACCGCTCACACCCTTGCTGACACCAAGCAACACATGATCCATCAGTGCAGCCCCAACGTTTTTTGTGACCTCCGTCCGAAAGTATGGATTTGAGCCTCAGCAATGTGGCCCACCATGGCCGAGCTACGCTCCCCTACGTACGGATGATTTCCCCGCCGGGACAGGCGGGACGGTTCTATTAAACATCTCACCTTACTGATGTCGCCCCCGGGTACGGCAGCGACGTACAGCCGCACGCGAGCTTGGTAAGGAGGAGCCAGATCGTGCCTACATGTTGAGTAGGTCCCTACTC");
            locus.put("Q", "CCATACGATGGGCTCGGCTCGTATGATTGAGGGCACCGAAGCTAGGCGACTCAAGATGGGCGAGGGCTCCGCGAATACGGAAAAGGGTATGCAGCTTCGGCCATAGGACCACGTGATCTCCGGGACTCGCCCAATTGAGATCGGCGTTACTACAATTAGACGAGCACATTCGCCAATTTATAGTTAACAAGAGATCTGCGTAGCAAAAGATTCAACGTTACCCGGCGCTACCAGACCCGCGGCACAGGCCGTTTGCCCTAAGCGTTGACAATCTTCGGATTCATCGCCAAAGTGCTGGAGTACAGGCTTCACACCTACAGCCACTGACAGCCTAACGAATCACTTCACCAATTGCCTTTCACCCCTGTTAGCGGAGGGCGCTAGCATAACTTTCGTCCTACGTTCCTCTCGTACGGATTCGGACAGATCCTCCGAGCAAGCGAAGTCCGGCCGATACTCTGCCCCTAGCAAGGCGGGATCCCGTCGCCTTGTAGTGAGCAGATATCCAGTTTGGGAGCACATAGAGCACCGACCGTCCACAGTCCCCTTCTCTTCATTGGTGCGTTGAGGTGAAAGTTCCCGAAATCCTACAAACTGGTATGCAAACCGGATCTGCAGGTTGGCTACGTATCATACTACAGCCCGTAAACTCGGAGTGCATGGCTTGACTAACATCCGTACAAACAGCGCTCGATGAGTGCGACCTCGGCAGCTTACACCCTTGCTGACACCAAATAGCGCATGATCCACCAGTACAGCCCAAGCGCCTTTCGCGTCCTCCGCCCGAATGTATGGATGTAAGCCTGAGTCACGTGGACCACCGTGCCCGAGCTACGCTCCCTTACGTGCGGATGATGTCCCCGCCGGGACAGGCGGAACGCTTCTATTGAACATTTCACCTCGCTGAAGTCGCCCCTGAGTACGCCAGCAACGTACAGTCGCACGCGAGCTCAGTAAGGAGGAGCCATATCGCGTCTACATGATGCGTAGGTCCCTGCGC");
            locus.put("L", "CTACGCGACGGGCTCGGCTCGTACAATAAAGGGCACCAGAGCTAGGTAACCCAAGATAGGCGAGGGCTTCGCTAATACGGAGAAGGGTATTCAGCTTCGGTCACAGGGCCACGTGATCTCCGGGACTCGCCCAGTTGCGATCGGCGTTACTACAATTGGAGGGGCACATTCGCCAACTTACAGTTAGTGAGGGATCTGCGTATCAAAAAACTCCACACTACTCGGCGCGACCTGGCTCGTGGCGCAGGCCGTTTGCCGTAAGTGTTGAGAATTTTCGGAAGAATGGCCGGAGTGTTAGAGTCCAGGCTTCACACCTACAGCCATTGACAGCTCAACGAGTCATTTCACCAGTTGCCTTTCACCCCTGTTAGCAGAGGGAGCTAGCGCAGCTTTCATCCTATGTCCATCTCGGGCGGACTCGGACAGATCCTCCGAACAAGCGAGGTCCGGTCGATACTCTGCCCCTAACAAGGCAGGATCCCATCGCCTTGTAGTGAGCAGATATCCAGTTTAGGGGCGCATAGAGCACCCACCGTCCACAGCCCCCTTGTACTTATTGGTGTGTTGAGGTGAAAGTCCCCGAAATCCCAGGAACTTGTATGCAACCCAGATCTGCAGGTAGGCTACGTATCATACTGTACCCCGCAAACTCGGAGTGTGTGACTGGACTAACGTCCGTACAAACAGCGCTCGATAGGTGCGACCTCGACAGCTTACACCCTTGCTGGCACCAAATAGCGCGTGATCCACCAGTGCAGACCAAACGTTTTCTGCGCCCTCCGTCCGAAAGTACGGGTTTAAGTCTCAGTAACATGGTCCACTATGACCGAACTATACCCCCTTACGTACGGGTGATGTCCCCGCCGGGACGGGCGGGACGCTTCTACTGAACGCTTCACTTCACTGAAGTCGCCCCTGAGTATGCCAGCAACGTACAGTCGCACGCGAGCTCGGTAAAGAGGAACCAGATCCCGTCGACATGTTGAGTTGGTCCCTATGC");
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
            int diff_topo = 0;
            for(int i = 0; i < runs; i++) {
                Operator op = new CoordinatedExchange(tree);
                double logHR = op.propose();
                if(logHR != Utils.INVALID_MOVE) counter++;
                if(tree.compareTo(template) != 0) {
                    diff_topo++;
                }

                System.out.println(tree.toString());
                op.undo();
                System.out.println(tree.toString());
                System.out.println("---------");
                if(tree.compareTo(template) == 0) {
                    test++;
                }

            }
            System.out.println(test == runs);
            System.out.printf("different topologies %d out of %d\n", diff_topo, runs);
            System.out.printf("%d out of %d\n", counter, runs);
            System.out.println(tree.getTree().toNewick());
        }
    }

    public static void main(String[] args) {

//        testUndo();
        testUndo();


    }

}
