package edu.rice.cs.bioinfo.programs.phylonet.algos.SNAPP;

import convolutionlib.JNIConvolution;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.algorithm.FelsensteinAlgorithm;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations.NucleotideObservation;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import jeigen.DenseMatrix;
import org.apache.commons.math3.util.ArithmeticUtils;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 * This class holds all the various algorithms from http://mbe.oxfordjournals.org/content/29/8/1917.
 */
public class Algorithms
{
    private static final boolean PRINT_DETAILS = false;
    private static final boolean USING_CONVOLUTIONLIB = false;
    private static final int MAXIMUM_NUMBER_OF_TAXA_FOR_CONVOLUTIONLIB = 6;

    /**
     * This class represents the data SNAPP associates with each branch.
     * Note that Tree nodes correspond to the branches above them.
     * IE: (A,B)C.
     * Node A would hold the data for branch A-C.
     * Node B would hold the data for branch B-C.
     * Node C would hold the data for branch C-infinity.
     */
    static class SNAPPData
    {
        /**
         * The maximum number of lineages (aka taxa) at this node.
         */
        int mx;

        /**
         * The "nucleotideIndex" of a node.
         * This index is what is used in R.getNum(index).
         */
        int nucleotideIndex;

        /**
         * The F Bottom matrix as described in the article.
         */
        FMatrix FBottom;

        /**
         * The F Top matrix as described in the article.
         */
        FMatrix FTop;
    }



    /**
     * Calculates an element of the F bottom matrix for a leaf.
     *
     * @param n    The number of lineages used to index the F matrix.
     * @param r    The R vector used to index the F matrix.
     * @param node The node to test.
     * @return The F matrix value for the given n and r for that node.
     */
    private static double getFBottomLeaf(int n, R r, STINode<SNAPPData> node)
    {
        if ((n == 1) && (r.getType() == node.getData().nucleotideIndex))
        {
            return 1.0;
        } else return 0.0;
    }


    /**
     * Calculates an element of the F bottom matrix for a non-leaf node
     *
     * @param nBot The number of lineages used to index the F matrix.
     * @param rBot The R vector used to index the F matrix.
     * @param node The node to test.
     * @return The F matrix value for the given n and r for that node.
     */
    private static double getFBottom(int nBot, R rBot, STINode<SNAPPData> node)
    {
        Iterator<STINode<SNAPPData>> children = node.getChildren().iterator();
        STINode<SNAPPData> leftChild = children.next();
        STINode<SNAPPData> rightChild = children.next();

        double sum = 0;
        for (int nTop = 1; nTop <= nBot - 1; nTop++)
        {

            if (nTop > leftChild.getData().mx || (nBot - nTop) > rightChild.getData().mx)
                continue;

            for (R rTop : R.loopOver(nTop))
            {

                R otherRTop = rBot.subtract(rTop);

                if (otherRTop == null)
                    continue;

                sum += leftChild.getData().FTop.get(rTop) *
                        rightChild.getData().FTop.get(otherRTop) *
                        rBot.getProbabilityOfSelecting(rTop);
                //System.out.println(leftChild.getData().FTop.get(rTop) + "*" + rightChild.getData().FTop.get(otherRTop) + "*" +rBot.getProbabilityOfSelecting(rTop));
            }
        }
        return sum;
    }

    /**
     * This transforms a nucleotide observation map of strings to characters to a "nucleotideIndex" map of strings to ints.
     *
     * @param obs The nucleotide observations.
     * @return The nucleotideIndex observations.
     */
    private static Map<String, Integer> getNucleotideIndexMap(NucleotideObservation obs)
    {
        Map<String, Integer> colorMap = new HashMap<String, Integer>();
        for (String allele : obs.getAlleles())
        {
            colorMap.put(allele, alleleToColor(obs.getObservationForAllele(allele)));
        }
        return colorMap;
    }

    /**
     * Converts a nucleotide to a "nucleotideIndex".
     *
     * @param character The nucleotide to convert.
     * @return The "nucleotideIndex".
     */
    private static int alleleToColor(char character)
    {
        if (R.getNumberOfTypes() == 2)
        {
            switch (character)
            {
                case '0':
                    return 0;
                case '1':
                    return 1;
                default:
                    throw new RuntimeException("Invalid nucleotide only 0 or 1 allowed. Color was: " + character);
            }
        } else
        {
            int number = FelsensteinAlgorithm.nucleotides.indexOf(character);

            if (number == -1)
                throw new RuntimeException("Bad nucleotide was " + character);
            return number;
        }
    }

    /**
     * This is the only public method of this file.
     * It calculates the probabilty of observing a species tree given some observations and a MatrixQ.
     * In other words, it calculates P(obs | speciesTree, Q)
     *
     * @param speciesTree The species tree.
     * @param obs         The nucleotide observations.
     * @param Q           The transition matrix for SNAPP.
     * @return The probability.
     */
    @SuppressWarnings("unchecked")
    public static double getProbabilityObservationGivenTree(Tree speciesTree, NucleotideObservation obs, MatrixQ Q)
    {
        Map<String, Integer> nucleotideIndexMap = getNucleotideIndexMap(obs);

        STITree<SNAPPData> castedSpeciesTree = (STITree<SNAPPData>) speciesTree;

        for (TNode node : castedSpeciesTree.postTraverse())
        {
            STINode<SNAPPData> actualNode = (STINode<SNAPPData>) node;
            processNode(Q, nucleotideIndexMap, actualNode);
        }

        return getProbabilityOfTree(castedSpeciesTree, Q);
    }

    /**
     * This method calculates the probabilty of data given a processed species tree.
     *
     * @param speciesTree An already processed species tree.
     * @param Q           A SNAPP transition matrix.
     * @return The probability
     */
    private static double getProbabilityOfTree(STITree<SNAPPData> speciesTree, MatrixQ Q)
    {
        STINode<SNAPPData> root = speciesTree.getRoot();

        DenseMatrix eq = Q.getEquilibrium();
        double sum = 0;
        for (int n = 1; n <= root.getData().mx; n++)
        {
            for (R r : R.loopOver(n))
            {
                sum += root.getData().FBottom.get(r) * eq.get(r.getIndex(), 0);
            }
        }
        return sum;
    }

    /**
     * Process a node for the snapp algorithm.
     *
     * @param Q                  The transition matrix for SNAPP.
     * @param nucleotideIndexMap The observations encoded as indices.
     * @param node               The node to process.
     */
    private static void processNode(MatrixQ Q, Map<String, Integer> nucleotideIndexMap, STINode<SNAPPData> node)
    {
        if(PRINT_DETAILS){
            System.out.println("\nNode " + node.getName());
        }

        SNAPPData data = new SNAPPData();
        node.setData(data);

        if (node.isLeaf())
        {
            processLeaf(nucleotideIndexMap, node, data);
        } else
        {
            processInternalNode(node, data);
        }

        if(PRINT_DETAILS){
            System.out.println("FBottom:");
            System.out.println(data.FBottom);
        }
        processTopBranchOfNode(Q, node, data);
        if(PRINT_DETAILS){
            System.out.println("FTop:");
            System.out.println(data.FTop);
        }
    }

    /**
     * Process the top branch of a node.
     *
     * @param Q    The transition matrix for SNAPP.
     * @param node The node to process.
     * @param data The data for said node.
     */
    private static void processTopBranchOfNode(MatrixQ Q, STINode<SNAPPData> node, SNAPPData data)
    {
        data.FTop = new FMatrix(data.mx);

        if (Double.isNaN(node.getParentDistance()) || Double.isInfinite(node.getParentDistance()))
            throw new RuntimeException("Snapp only works with finite branch distances: " + node.getParentDistance());

        DenseMatrix result = Q.getProbabilityForColumn(node.getParentDistance(), data.FBottom.arr);

        data.FTop.setMatrix(result);
    }

    /**
     * Process an internal non-leaf node.
     *
     * @param node The node to process.
     * @param data The data for that node.
     */
    private static void processInternalNode(STINode<SNAPPData> node, SNAPPData data)
    {
        Iterator<STINode<SNAPPData>> children = node.getChildren().iterator();
        if (node.getChildCount() != 2)
            throw new RuntimeException("SNAPP does not work on trees with 2 or more children per node");

        STINode<SNAPPData> leftChild = children.next();
        STINode<SNAPPData> rightChild = children.next();

        data.mx = leftChild.getData().mx + rightChild.getData().mx;


        if (USING_CONVOLUTIONLIB)
            data.FBottom = getFBottomConvolution(data, leftChild, rightChild);
        else
            data.FBottom = getFBottomNormal(node, data);
    }

    /**
     * The slow normal calcluation for F bottom.
     *
     * @param node The node to process.
     * @param data The data for that node.
     * @return The F bottom matrix.
     */
    private static FMatrix getFBottomNormal(STINode<SNAPPData> node, SNAPPData data)
    {
        FMatrix FBottom = new FMatrix(data.mx);
        for (int n = 1; n <= data.mx; n++)
        {
            for (R r : R.loopOver(n))
                FBottom.set(r, getFBottom(n, r, node));
        }

        return FBottom;
    }

    /**
     * Perform the initial processing of a leaf.
     *
     * @param nucleotideIndexMap The observations encoded as indices.
     * @param node               The node to process.
     * @param data               The data for that node.
     */
    private static void processLeaf(Map<String, Integer> nucleotideIndexMap, STINode<SNAPPData> node, SNAPPData data)
    {
        data.mx = 1;
        data.nucleotideIndex = nucleotideIndexMap.get(node.getName());

        data.FBottom = new FMatrix(data.mx);
        for (int n = 1; n <= data.mx; n++)
        {
            for (R r : R.loopOver(n))
                data.FBottom.set(r, getFBottomLeaf(n, r, node));
        }
    }

    /**
     * This is the local JNIConvolution available.
     */
    private static final ThreadLocal<JNIConvolution> myConvolution =
            new ThreadLocal<JNIConvolution>()
            {
                @Override
                protected JNIConvolution initialValue()
                {
                    return createConvolution();
                }
            };


    /**
     * Create a new JNIConvolution.
     * Synchronized as the JNIConvlution constructor is not thread safe
     * (because fftw's contstructor is not thread safe.)
     * @return A new JNIConvolution
     */
    private static synchronized JNIConvolution createConvolution()
    {
        JNIConvolution foo = new JNIConvolution(R.getNumberOfTypes());
        for (int n = 1; n <= MAXIMUM_NUMBER_OF_TAXA_FOR_CONVOLUTIONLIB; n++)
            foo.prepareEngine(n + 1);
        return foo;
    }

    /**
     * Calculates the F Bottom matrix for a node using the faster fourier transform/convolution method.
     *
     * @param data       The data for that node.
     * @param leftChild  The left child.
     * @param rightChild The right child.
     * @return The F Bottom matrix.
     */
    private static FMatrix getFBottomConvolution(SNAPPData data, STINode<SNAPPData> leftChild, STINode<SNAPPData> rightChild)
    {
        int size = data.mx + 1;
        double[] arr1 = new double[ArithmeticUtils.pow(size, R.getNumberOfTypes())];
        double[] arr2 = new double[ArithmeticUtils.pow(size, R.getNumberOfTypes())];


        for (int n = 1; n <= data.mx; n++)
        {
            for (R r : R.loopOver(n))
            {
                if (n <= leftChild.getData().mx)
                    arr1[r.getBoxIndex(size)] = leftChild.getData().FTop.get(r) * r.getLikelihoodProduct();

                if (n <= rightChild.getData().mx)
                    arr2[r.getBoxIndex(size)] = rightChild.getData().FTop.get(r) * r.getLikelihoodProduct();
            }
        }


        double[] result = myConvolution.get().convolute(arr1, arr2, size);


        FMatrix otherMethod = new FMatrix(data.mx);

        for (int n = 1; n <= data.mx; n++)
        {
            for (R r : R.loopOver(n))
            {
                otherMethod.set(r, result[r.getBoxIndex(size)] / r.getLikelihoodProduct());
            }
        }
        return otherMethod;
    }

}




