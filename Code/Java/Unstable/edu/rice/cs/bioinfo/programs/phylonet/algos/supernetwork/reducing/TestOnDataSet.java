package edu.rice.cs.bioinfo.programs.phylonet.algos.supernetwork.reducing;

import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;

public class TestOnDataSet {

  /**
   * Root an unrooted tree.
   */
  public static String rootTree(String unrootedTree) {
    Tree gt_ = Trees.readTree(unrootedTree);

    STITree gt = (STITree) gt_;

    int minLeafCount = Integer.MAX_VALUE;
    STINode childToReroot = null;
    for (Object child:gt.getRoot().getChildren()) {
      STINode stiChild = (STINode) child;
      if (stiChild.getLeafCount() < minLeafCount) {
        minLeafCount = stiChild.getLeafCount();
        childToReroot = stiChild;
      }
    }

    childToReroot.removeNode();

    STITree newTree = new STITree();
    newTree.getRoot().adoptChild(gt.getRoot());
    newTree.getRoot().adoptChild(childToReroot);

    return (newTree.toString());
  }

//  public static void testTrinet() {
//    GenerateTrinets gen = new GenerateTrinets();
//
//    File data = new File("/Users/liu/Desktop/testdata/trinets/Reti5_B_gts.tre");
//
//    try (BufferedReader br = new BufferedReader(new FileReader(data))) {
//      String line;
//      while ((line = br.readLine()) != null) {
//        gen.addTree(rootTree(line));
//      }
//    } catch (IOException e) {
//      e.printStackTrace();
//    }
//
//    int count = 0;
//    for (String[] trinet:gen.generate()) {
//      System.out.println(Arrays.toString(trinet));
//      count++;
//    }
//
//    System.out.println(count);
//  }

  public static void test3SetTrinet() {
    GenerateTrinets3Sets gen = new GenerateTrinets3Sets();

    File data = new File("/Users/liu/Desktop/testdata/trinets/all/Reti5_D/gts.tre");

    try (BufferedReader br = new BufferedReader(new FileReader(data))) {
      String line;
      while ((line = br.readLine()) != null) {
        gen.addTree(rootTree(line));
      }
    } catch (IOException e) {
      e.printStackTrace();
    }

    int count = 0;
    for (String[] trinet:gen.generate()) {
      System.out.println(Arrays.toString(trinet));
      count++;
    }

    System.out.println(count);
  }

  public static void main(String[] args) {
    test3SetTrinet();
  }
}
