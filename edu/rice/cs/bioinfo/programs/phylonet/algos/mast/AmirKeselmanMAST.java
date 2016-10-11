package edu.rice.cs.bioinfo.programs.phylonet.algos.mast;

import edu.rice.cs.bioinfo.programs.phylonet.algos.SymmetricDifference;
import edu.rice.cs.bioinfo.programs.phylonet.algos.lca.SchieberVishkinLCA;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.MutableTree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 9/25/16
 * Time: 2:59 PM
 * To change this template use File | Settings | File Templates.
 */
public class AmirKeselmanMAST {

    private boolean _verbose = false;
    private Iterable<? extends Tree> _trees;

    public Tree computeRMAST(Iterable<? extends Tree> trees) {

        int degree = 2;
        String [] leaves = null;
        _trees = trees;

        for(Tree t : trees) {
            Trees.removeBinaryNodes((MutableTree) t);
            if(leaves == null) {
                leaves = t.getLeaves();
                Arrays.sort(leaves);
            } else {
                String [] curLeaves = t.getLeaves();
                Arrays.sort(curLeaves);
                if(!Arrays.equals(curLeaves, leaves)) {
                    throw new RuntimeException("Leaves inconsistant.");
                }
            }
            for(TNode node : t.postTraverse()) {
                if(node.getChildCount() > degree) {
                    throw new RuntimeException("Only implemented for binary trees.");
                }
            }
        }

        //Construct all Maximal Decomposable Subsets
        Set<String> maximalDecomposableSet = new HashSet<>();
        Map<Set<String>, List<Set<String>>> decompose = new HashMap<>();
        Map<Set<String>, Integer> decomposeHeightCount = new HashMap<>();
        for(String a1 : leaves) {
            for(String a2 : leaves) {
                if(a1.compareTo(a2) >= 0) continue;
                Set<String> s1 = new HashSet<>();
                Set<String> s2 = new HashSet<>();
                for(Tree tree : trees) {
                    List<TNode> children = new ArrayList<>();
                    getLCA(tree.getRoot(), a1, a2, children);
                    Set<String> st1 = new HashSet<>();
                    Set<String> st2 = new HashSet<>();
                    getSubtreeLeaves(children.get(0), st1);
                    getSubtreeLeaves(children.get(1), st2);

                    if(s1.isEmpty() && s2.isEmpty()) {
                        s1 = st1;
                        s2 = st2;
                    } else {
                        s1.retainAll(st1);
                        s2.retainAll(st2);
                    }
                }

                Set<String> s = new HashSet<>();
                s.addAll(s1);
                s.addAll(s2);
                if(!decompose.containsKey(s)) {
                    List<Set<String>> setList = new ArrayList<>();
                    setList.add(s1);
                    setList.add(s2);
                    decompose.put(s, setList);
                    decomposeHeightCount.put(s, 0);
                    for(Tree tree : trees) {
                        int count = getDecomposeLinkCount(tree, s);
                        if(count > decomposeHeightCount.get(s)) {
                            decomposeHeightCount.put(s, count);
                        }
                    }
                }

                if(s.size() > maximalDecomposableSet.size()
                        || (s.size() == maximalDecomposableSet.size()
                                && decomposeHeightCount.get(s) < decomposeHeightCount.get(maximalDecomposableSet))) {
                    maximalDecomposableSet.clear();
                    maximalDecomposableSet.addAll(s);
                }

                if(_verbose) {
                    System.out.println("Tuple: " + a1 + " " + a2);
                    System.out.println("Decompose: " + setToString(s1) + ", " + setToString(s2));
                    //System.out.println("Count: " + decomposeHeightCount.get(s));
                }
            }
        }
        if(_verbose)
            System.out.println("Maximal Decomposable Set: " + setToString(maximalDecomposableSet));

        //construct the tree
        Tree mast = constructTreeFromDecomposition(decompose, decomposeHeightCount, maximalDecomposableSet);

        return mast;
    }

    private String setToString(Set<String> set) {
        StringBuilder sb = new StringBuilder();
        for(String str : set) {
            sb.append(str);
            sb.append(" ");
        }
        return sb.toString();
    }

    private int getDecomposeHeight(Tree tree, Set<String> leaves) {
        int count = 0;
        for(String leafName : leaves) {
            TNode node = tree.getNode(leafName);
            while(!node.isRoot()) {
                node = node.getParent();
                count++;
            }
        }
        return count;
    }

    private int getDecomposeLinkCount(Tree tree, Set<String> leaves) {
        SchieberVishkinLCA lcaClass = new SchieberVishkinLCA(tree);
        Set<TNode> leafNodes = new HashSet<>();
        for(String leafName : leaves) {
            leafNodes.add(tree.getNode(leafName));
        }
        TNode lcaNode = lcaClass.getLCA(leafNodes);
        int OUT = 0;
        int IN = 1;
        Map<TNode, String> marks = new HashMap<>();

        Set<TNode> otherLeafNodes = new HashSet<>();

        for(TNode node : lcaNode.postTraverse()) {
            marks.put(node, "");
            if(node.isLeaf() && !leafNodes.contains(node))
                otherLeafNodes.add(node);
        }

        for(TNode leaf : leafNodes) {
            TNode node = leaf;
            marks.put(node, leaf.getName());
            while(node != lcaNode){
                node = node.getParent();
                if(marks.get(node).equals(""))
                    marks.put(node, leaf.getName());
                else marks.put(node, marks.get(node) + leaf.getName());
            };
        }

        Map<String, Integer> linkCount = new HashMap<>();
        Set<TNode> counted = new HashSet<>();

        for(TNode node : otherLeafNodes) {
            do {
                if(marks.get(node).length() > 0){
                    if(counted.contains(node))
                        break;
                    counted.add(node);
                    if(!linkCount.containsKey(marks.get(node))){
                        linkCount.put(marks.get(node), 0);
                    }
                    linkCount.put(marks.get(node), linkCount.get(marks.get(node)) + 1);
                    break;
                }
                node = node.getParent();
            }while(node != lcaNode);
        }

        if(linkCount.size() == 0)
            linkCount.put("", 0);

        return Collections.max(linkCount.values());
    }

    private int computeDecomposeLinkCount(Set<String> greaterSet, Set<String> innerSet) {
        int ret = 0;
        Set<String> tmp = new HashSet<>();
        tmp.add("4");
        tmp.add("6");
        tmp.add("7");
        if(greaterSet.equals(tmp))
            System.out.println("!");
        for(Tree tree : _trees) {
            SchieberVishkinLCA lcaClass = new SchieberVishkinLCA(tree);
            Set<TNode> greaterLeafNodes = new HashSet<>();
            for(String leafName : greaterSet) {
                greaterLeafNodes.add(tree.getNode(leafName));
            }
            TNode lcaNode = lcaClass.getLCA(greaterLeafNodes);

            Set<TNode> innerLeafNodes = new HashSet<>();
            for(String leafName : innerSet) {
                innerLeafNodes.add(tree.getNode(leafName));
            }

            Map<TNode, String> marks = new HashMap<>();

            for(TNode node : lcaNode.postTraverse()) {
                marks.put(node, "");
            }

            for(TNode leaf : innerLeafNodes) {
                TNode node = leaf;
                marks.put(node, leaf.getName());
                while(node != lcaNode){
                    node = node.getParent();
                    if(marks.get(node).equals(""))
                        marks.put(node, leaf.getName());
                    else marks.put(node, marks.get(node) + "+" + leaf.getName());
                };
            }

            Map<String, Integer> linkCount = new HashMap<>();
            Set<TNode> counted = new HashSet<>();

            for(TNode node : greaterLeafNodes) {
                if(!innerLeafNodes.contains(node)) {
                    while (node != lcaNode) {
                        node = node.getParent();
                        if (marks.get(node).length() > 0) {
                            if (counted.contains(node))
                                break;
                            counted.add(node);
                            if (!linkCount.containsKey(marks.get(node))) {
                                linkCount.put(marks.get(node), 0);
                            }
                            linkCount.put(marks.get(node), linkCount.get(marks.get(node)) + 1);
                            break;
                        }
                    }
                }
            }

            if(linkCount.size() == 0)
                linkCount.put("", 0);

            ret = Math.max(ret, Collections.max(linkCount.values()));
        }


        return ret;
    }

    private STITree<Double> constructTreeFromDecomposition(Map<Set<String>, List<Set<String>>> decompose, Map<Set<String>, Integer> decomposeHeightCount, Set<String> cur) {
        if(cur.size() == 1) {
            STITree<Double> tree = new STITree<>();
            tree.getRoot().setName(cur.iterator().next());
            return tree;
        }

        if(cur.size() == 2) {
            STITree<Double> tree = new STITree<>();
            for(String leafName : cur) {
                tree.getRoot().createChild(leafName);
            }
            return tree;
        }

        Set<String> leftLeaves = null;
        if(decompose.get(cur).get(0).size() <= 2) {
            leftLeaves = decompose.get(cur).get(0);
        } else {
            for (Set<String> subset : decompose.keySet()) {
                if (decompose.get(cur).get(0).containsAll(subset)) {
                    if (leftLeaves == null || subset.size() > leftLeaves.size()
                            || (subset.size() == leftLeaves.size()
                                && computeDecomposeLinkCount(decompose.get(cur).get(0), subset) < computeDecomposeLinkCount(decompose.get(cur).get(0), leftLeaves))) {
                        leftLeaves = subset;
                    }
                }
            }
        }

        Set<String> rightLeaves = null;
        if(decompose.get(cur).get(1).size() <= 2) {
            rightLeaves = decompose.get(cur).get(1);
        } else {
            for (Set<String> subset : decompose.keySet()) {
                if (decompose.get(cur).get(1).containsAll(subset)) {
                    if (rightLeaves == null || subset.size() > rightLeaves.size()
                            || (subset.size() == rightLeaves.size()
                                && computeDecomposeLinkCount(decompose.get(cur).get(1), subset) < computeDecomposeLinkCount(decompose.get(cur).get(1), rightLeaves))){
                        rightLeaves = subset;
                    }
                }
            }
        }

        STITree<Double> leftTree = constructTreeFromDecomposition(decompose, decomposeHeightCount, leftLeaves);
        STITree<Double> rightTree = constructTreeFromDecomposition(decompose, decomposeHeightCount, rightLeaves);

        STITree<Double> tree = new STITree<>();
        tree.getRoot().adoptChild(leftTree.getRoot());
        tree.getRoot().adoptChild(rightTree.getRoot());
        return tree;
    }

    private void getSubtreeLeaves(TNode node, Set<String> leaves) {
        for(TNode n : node.postTraverse()) {
            if(n.isLeaf())
                leaves.add(n.getName());
        }
    }

    private TNode getLCA(TNode node, String a1, String a2) {
        List<TNode> children = new ArrayList<>();
        return getLCA(node, a1, a2, children);
    }

    private TNode getLCA(TNode node, String a1, String a2, List<TNode> children) {
        if(node.getName().equals(a1)) return node;
        if(node.getName().equals(a2)) return node;
        if(node.isLeaf()) return null;

        List<TNode> subtree = new ArrayList<TNode>();
        for(TNode child : node.getChildren()) {
            subtree.add(child);
        }
        TNode left = getLCA(subtree.get(0), a1, a2, children);
        TNode right = getLCA(subtree.get(1), a1, a2, children);
        if(left != null && right != null) {
            if(left.getName().equals(a1)) {
                children.add(subtree.get(0));
                children.add(subtree.get(1));
            } else {
                children.add(subtree.get(1));
                children.add(subtree.get(0));
            }
            return node;
        }

        if(left != null) return left;
        else return right;
    }

    /*private int getSubtreeRetain(Tree tree, List<String> leaves) {

        Map<TNode, Integer> marks = new HashMap<>();

        int OUT = 1;
        int IN = 2;

        for(TNode node : tree.postTraverse()) {
            marks.put(node, IN);
        }

        marks.put(tree.getRoot(), OUT);

        for(String leafName : tree.getLeaves()) {
            TNode leaf = tree.getNode(leafName);
            if(!leaves.contains(leaf.getName())) {
                TNode node = leaf;
                while(!node.isRoot()) {
                    marks.put(node, OUT);
                    node = node.getParent();
                }
            }
        }

        int count = 0;
        for(TNode node : tree.postTraverse()) {
            if(!node.isRoot()) {
                if(!marks.get(node).equals(marks.get(node.getParent()))) {
                    count++;
                }
            }
        }

        boolean done = false;
        while(!done) {
            done = true;
            for (TNode node : tree.postTraverse()) {
                if (!leaves.contains(node.getName()) && node.getChildCount() == 0) {
                    TNode parent = node.getParent();
                    parent.re
                    for (NetNode<Object> parent : node.getParents()) {
                        parent.removeChild(node);
                        done = false;
                    }
                }
            }
        }

        return count;
    }*/
}
