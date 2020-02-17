package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.util;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.TreeNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.ParseException;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.IOException;
import java.math.BigDecimal;
import java.util.*;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

public class trash {
    public static void main(String[] args) throws IOException, ParseException {
//        STITree<TreeNodeInfo> tree = model.getTree();
//        for (TNode node:tree.postTraverse()) {
//            STINode<TreeNodeInfo> stiNode = (STINode<TreeNodeInfo>) node;
//            System.out.println(stiNode.getData().getPopSize());
//        }
        double generationHeight = 45000;
        int N0 = 10000;
        double coalescentHeight = generationHeight/N0;
        System.out.println(coalescentHeight);

        int popSize = 25000;
        double relativePopsize = (double)popSize / N0;
        System.out.println(relativePopsize);


        STITree<TreeNodeInfo> t1 = new STITree<>("((1:37.422,2:37.422):36.297,(3:37.819,4:37.819):35.899);");
        STITree<TreeNodeInfo> t2 = new STITree<>("((3:36.195,4:36.195):37.524,(1:41.137,2:41.137):32.582);");
//        System.out.println(Trees.getLexicographicNewickString(t1, null));
//        System.out.println(Trees.getLexicographicNewickString(t2, null));
        Trees.convertToLexicographicTree(t1);
        Trees.convertToLexicographicTree(t2);
        System.out.println(t1.toNewick());
        System.out.println(t2.toNewick());

        Set<Integer> s1 = new HashSet<>();
        s1.add(1);
        s1.add(2);
        Set<Integer> s2 = new HashSet<>();
        s2.add(1);
        s2.add(2);
        s2.add(3);
        System.out.println(s2.containsAll(s1));


        List<Integer> l1 = new LinkedList<>();
        l1.add(1);
        l1.add(2);
        l1.add(10);
        List<Integer> l2 = new LinkedList<>();
        l2.add(1);
        l2.add(2);
        l2.add(10);
        System.out.println(l1.equals(l2));


        double[][] mat = new double[][]{{1.1,2.2,3.3}, {4.4,5.5,6.6}, {7.7,8.8,9.9}};
        System.out.println(Arrays.deepToString(mat));
        double[] ar = Arrays.stream(mat).mapToDouble(arr -> DoubleStream.of(arr).sum()).toArray();
        System.out.println(Arrays.toString(ar));
        System.out.println(DoubleStream.of(mat[0]).sum());
        System.out.println(DoubleStream.of(mat[1]).sum());
        System.out.println(DoubleStream.of(mat[2]).sum());



        System.out.println(System.getProperty("java.only"));



        STITree<TreeNodeInfo> tt1 = new STITree<>("[&U] (3.3:283041.51180, (2.2:203562.74288, 1.1:203562.74288):79478.76892);");
        System.out.println(tt1.isRooted());
        System.out.println(tt1.toNewick());
        Trees.convertToLexicographicTree(tt1);
        System.out.println(tt1.toNewick());
        STITree tt2 = new STITree("[&U] (3.3:283041.51180, (2.2:203562.74288, 1.1:203562.74288):79478.76892);");
        System.out.println(Trees.haveSameRootedTopology(tt1, tt2));
        for (TNode node:tt1.postTraverse()) {
            if (node.isRoot()) {
                System.out.println("!!!!!!");
            }
        }

        double num = 1.5E-8;
        System.out.println(num);
        System.out.println(BigDecimal.valueOf(num).toPlainString());

    }
}
