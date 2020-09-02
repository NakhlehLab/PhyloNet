package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.util;

import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.ModelTree;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.structs.TreeNodeInfo;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.test.HCGModelBuilder;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variational.distribution.Prior;
import edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCcoal.variationalMulti.MultiVariateGaussian;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.ParseException;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;
import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

import java.io.IOException;
import java.math.BigDecimal;
import java.util.*;
import java.util.stream.DoubleStream;

public class trash {
    public static void main(String[] args) throws IOException, ParseException {
//        testPopsizePrior();
//        NormalDistribution dist =  new NormalDistribution(0, 2);
//        System.out.println(dist.density(0.6));
        trash5();
    }

    public static void trash1() throws IOException, ParseException {
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

    public static void trash2() {
        NormalDistribution dist = new NormalDistribution(1, 0.25);

        System.out.println(dist.sample());
        System.out.println(dist.probability(1));
        System.out.println(dist.logDensity(1));

        double scaleFactor = 1E-8;



        System.out.println(8 * scaleFactor);

        System.out.println(1E-8/1E-8);
    }

    public static void testPopsizePrior() {
        ModelTree model = HCGModelBuilder.getHCGModelInit();
        Prior prior = new Prior();
        System.out.println(prior.logPrior(model));
    }

    public static void testLogTime() {
        long startTime = System.currentTimeMillis();
        System.out.println(startTime);
        trash2();
        long endTime = System.currentTimeMillis();
        System.out.println((double)(endTime - startTime) / 1000);
        
    }

    public static void trash3() {
        int b = 2;
        System.out.println(1.0 / b);
    }

    public static void trash4() {
        double[] mean = new double[]{0, 0, 0, 0};
        double[][] matrixArray = new double[][]{
                {2.6708760497203845E8, 6334072.357230105, -8.628748560974622E7, -4500290.9125669},
                {0.0, 3.3590023193667424E8, 1.749093196237955E7, -8.292803416448978E7},
                {0.0, 0.0, 4.801235808567999E7, -2.027825400179901E7},
                {0.0, 0.0, 0.0, 3.43132036513448E7}
        };
//        RealMatrix matrix = MatrixUtils.createRealMatrix(MultiVariateGaussian.full(matrixArray));
//        MatrixUtils.inverse(matrix);
        double[][] fullMatrix = MultiVariateGaussian.full(matrixArray);
        for (double[] row:fullMatrix) {
            System.out.println(Arrays.toString(row));
        }

        EigenDecomposition eigend = new EigenDecomposition(MatrixUtils.createRealMatrix(fullMatrix));
        System.out.println("Eigen values are: " + Arrays.toString(eigend.getRealEigenvalues()));

        MultivariateNormalDistribution dist = new MultivariateNormalDistribution(mean, fullMatrix);
    }

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

    public static void trash5() throws IOException, ParseException {
        STITree<TreeNodeInfo> tree = new STITree<>("((H:150000,C:150000):150000, G:300000);");
//        STITree<TreeNodeInfo> tree = new STITree<>("((A: 5, B: 5):5, (C: 5, D: 5):5);");
        System.out.println(tree.toNewick());
//        buildGTNodeHeight(tree);
        System.out.println(tree.getRoot().getChildren().iterator().next().getParentDistance());
    }

}
