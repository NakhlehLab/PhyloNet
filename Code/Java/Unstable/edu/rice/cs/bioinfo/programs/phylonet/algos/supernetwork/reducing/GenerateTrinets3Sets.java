package edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork.reducing;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.util.*;

/**
 * Created by Xinhao Liu.
 * Date: 11/29/18
 *
 * Generate set of trinets that cover all nodes in a set of input gene trees.
 * To cover a node, instead of choosing one leaf from one child and two leaves from the other child,
 * we choose one leaf from each child, and one leaf from parent's sibling. This way we also covers
 * all edges in the original network.
 */
public class GenerateTrinets3Sets {
  private ArrayList<Tree> gts = new ArrayList<>();

  /**
   * Constructor.
   */
  public GenerateTrinets3Sets() {

  }

  /**
   * Add a tree to the set of trees for generating trinets.
   *
   * @param gt string representation of gene tree
   */
  public void addTree(String gt) {
    gts.add(Trees.readTree(gt));
  }

  /**
   * Getter for all added gene trees.
   *
   * @return all added gene trees
   */
  public ArrayList<Tree> getGts() {
    return gts;
  }


  /**
   * Generate the set of trinets that cover all gene tree nodes.
   * @return set of trinets
   */
  public Set<String[]> generate() {
    HashSet<String[]> trinets = new HashSet<>();
    for (Tree gt:gts) {
      for (TNode node:gt.getNodes()) {
        if (!node.isLeaf() && !node.isRoot() && !isCovered(trinets, node)) {
          trinets.add(generateTrinetOfANode(node));
        }
      }
    }
    return trinets;
  }

  /**
   * Test if a node has already been covered by the given trinets.
   * @param trinets trinets already established
   * @param node new node
   * @return whether the node has already been covered
   */
  public static boolean isCovered(HashSet<String[]> trinets, TNode node) {
    ArrayList<TNode> allchildren = new ArrayList<>();
    HashSet<String> child1Set = new HashSet<>();
    HashSet<String> child2Set = new HashSet<>();

    for (TNode child:node.getChildren()) {
      allchildren.add(child);
    }
    // all leaves in two children
    for (TNode leaf:allchildren.get(0).getLeaves()) {
      child1Set.add(leaf.getName());
    }
    for (TNode leaf:allchildren.get(1).getLeaves()) {
      child2Set.add(leaf.getName());
    }

    // all leaves from sibling
    TNode otherSibling = null;
    for (TNode sibling:node.getParent().getChildren()) {
      if (!sibling.equals(node)) {
        otherSibling = sibling;
      }
    }
    HashSet<String> siblingSet = new HashSet<>();
    for (TNode leaf:otherSibling.getLeaves()) {
      siblingSet.add(leaf.getName());
    }

//    for (String[] trinet:trinets) {
//      Set<String> threeTaxonSet = new HashSet<>(Arrays.asList(trinet));
//      Set<String> intersectionWithChild1 = new HashSet<>(threeTaxonSet);
//      intersectionWithChild1.retainAll(child1Set);
//      Set<String> intersectionWithChild2 = new HashSet<>(threeTaxonSet);
//      intersectionWithChild2.retainAll(child2Set);
//      if (intersectionWithChild1.size() == 1 && intersectionWithChild2.size() == 1) {
//        return true;
//      }
//    }

    for (String[] trinet:trinets) {
      Set<String> threeTaxonSet = new HashSet<>(Arrays.asList(trinet));
      Set<String> intersectionWithChild1 = new HashSet<>(threeTaxonSet);
      intersectionWithChild1.retainAll(child1Set);
      Set<String> intersectionWithChild2 = new HashSet<>(threeTaxonSet);
      intersectionWithChild2.retainAll(child2Set);
      Set<String> intersectionWithSibling = new HashSet<>(threeTaxonSet);
      intersectionWithSibling.retainAll(siblingSet);
      if (intersectionWithChild1.size() == 1 && intersectionWithChild2.size() == 1 && intersectionWithSibling.size() == 1) {
        return true;
      }
    }
    return false;
  }


  public static String[] generateTrinetOfANode(TNode node) {
    ArrayList<TNode> allChildren = new ArrayList<>();
    for (TNode child:node.getChildren()) {
      allChildren.add(child);
    }

    TNode child1 = allChildren.get(0);
    TNode child2 = allChildren.get(1);

    ArrayList<String> child1Leaves = new ArrayList<>();
    ArrayList<String> child2Leaves = new ArrayList<>();

    // all leaves in two children
    for (TNode leaf:child1.getLeaves()) {
      child1Leaves.add(leaf.getName());
    }
    for (TNode leaf:child2.getLeaves()) {
      child2Leaves.add(leaf.getName());
    }

    ArrayList<String> trinet = new ArrayList<>();
    Random rand = new Random();
    trinet.add(child1Leaves.get(rand.nextInt(child1Leaves.size())));
    trinet.add(child2Leaves.get(rand.nextInt(child2Leaves.size())));

    TNode otherSibling = null;
    for (TNode sibling:node.getParent().getChildren()) {
      if (!sibling.equals(node)) {
        otherSibling = sibling;
      }
    }
    ArrayList<String> siblingLeaves = new ArrayList<>();
    for (TNode leaf:otherSibling.getLeaves()) {
      siblingLeaves.add(leaf.getName());
    }
    trinet.add(siblingLeaves.get(rand.nextInt(siblingLeaves.size())));

    return trinet.toArray(new String[0]);
  }

  /**
   * Generate a trinet that covers a node given its two children.
   * Randomly select a leaf from one child and two leaves from the other.
   *
   * @param node1 child of node to cover
   * @param node2 child of node to cover
   * @return a trinet that covers the node whose children has been given
   */
  public static String[] generateTrinetGivenTwoSets(TNode node1, TNode node2) {
    ArrayList<String> node1Leaves = new ArrayList<>();
    ArrayList<String> node2Leaves = new ArrayList<>();

    // all leaves in two children
    for (TNode leaf:node1.getLeaves()) {
      node1Leaves.add(leaf.getName());
    }
    for (TNode leaf:node2.getLeaves()) {
      node2Leaves.add(leaf.getName());
    }

    ArrayList<String> trinet = new ArrayList<>();
    if (node1Leaves.size() == 1) {    // if one of the children is a leaf
      trinet.add(node1Leaves.get(0));
      // Randomly select two distinct leaves from the other child
      ArrayList<Integer> indices = new ArrayList<>();
      for (int i = 0; i < node2Leaves.size(); i++) {
        indices.add(i);
      }
      Collections.shuffle(indices);
      trinet.add(node2Leaves.get(indices.get(0)));
      trinet.add(node2Leaves.get(indices.get(1)));
    } else if (node2Leaves.size() == 1) {    // if one of the children is a leaf
      trinet.add(node2Leaves.get(0));
      // Randomly select two distinct leaves from the other child
      ArrayList<Integer> indices = new ArrayList<>();
      for (int i = 0; i < node1Leaves.size(); i++) {
        indices.add(i);
      }
      Collections.shuffle(indices);
      trinet.add(node1Leaves.get(indices.get(0)));
      trinet.add(node1Leaves.get(indices.get(1)));
    } else {
      // Select a leaf from one child.
      Random rand = new Random();
      trinet.add(node1Leaves.get(rand.nextInt(node1Leaves.size())));

      // Randomly select two distinct leaves from the other child
      ArrayList<Integer> indices = new ArrayList<>();
      for (int i = 0; i < node2Leaves.size(); i++) {
        indices.add(i);
      }
      Collections.shuffle(indices);
      trinet.add(node2Leaves.get(indices.get(0)));
      trinet.add(node2Leaves.get(indices.get(1)));
    }

    return trinet.toArray(new String[0]);
  }


  public static void main(String[] args) {
    Tree t1 = Trees.readTree("((A,(B,C)),D);");

    Tree t2 = Trees.readTree("(((D,(F,E)),(C,(B,A))),(J,(K,(I,G))));");

//    HashSet<String[]> trinets = new HashSet<>();
//    trinets.add(generateTrinetOfANode(t2.getRoot().getChildren().iterator().next()));
//    trinets.add(generateTrinetOfANode(t2.getRoot()));
//    //System.out.println(Arrays.toString(generateTrinetOfANode(t2.getRoot().getChildren().iterator().next())));
//    //System.out.println(Arrays.toString(generateTrinetOfANode(t2.getRoot())));
//
//    System.out.println(isCovered(trinets, t2.getRoot()));
//    System.out.println(isCovered(trinets, t2.getRoot().getChildren().iterator().next()));

    GenerateTrinets3Sets gen = new GenerateTrinets3Sets();
    //gen.addTree("((A,(B,C)),D);");
    gen.addTree("(((D,(F,E)),(C,(B,A))),(J,(K,(I,G))));");

    for (String[] trinet:gen.generate()) {
      System.out.println(Arrays.toString(trinet));
    }
  }
}
