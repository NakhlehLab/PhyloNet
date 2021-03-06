/*
 * Copyright (c) 2012 Rice University.
 *
 * This file is part of PhyloNet.
 *
 * PhyloNet is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PhyloNet is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PhyloNet.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.rice.cs.bioinfo.programs.phylonet.algos.simulator;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.RichNewickReaderAST;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.parsers.antlr.ast.ANTLRRichNewickParser;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReadResult;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.NetworkFactoryFromRNNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.ByteArrayInputStream;
import java.util.*;

/**
 * Created by Yun Yu
 * Date: 11/2/11
 * Time: 10:47 AM
 *
 * This class is for simulating gene trees within branches of a given species networks
 * Note that the gene trees generated do NOT have branch lengths
 */
public class SimGTInNetwork
{
    private Map<NetNode, Map<Integer, double[]>> _node2gij = new HashMap<NetNode, Map<Integer, double[]>>();
    private boolean _printDetails = false;
    private Random _random;
    private Long _seed = null;

    /**
     * This function is to set seed to control the randomness
     */
    public void setSeed(Long seed){
        _seed = seed;
    }



    /**
     * This is the main function for generating gene trees given a species network
     *
     * @param network           the species network that the gene trees are generated from
     * @param species2alleles   mapping from species to alleles sampled from it
     * @param numGTs            the number of gene trees need to be generated
     *
     * @return   the list of simulated gene trees
     */
    public List<Tree> generateGTs(Network network, Map<String, List<String>> species2alleles, int numGTs){
        if(_seed == null){
            _random = new Random();
        }
        else{
            _random = new Random(_seed);
        }
        List<Tree> gts = new ArrayList<Tree>();
        if(species2alleles==null){
            species2alleles = new HashMap<String, List<String>>();
            for(Object leaf: network.getLeaves()){
                String species = ((NetNode)leaf).getName();
                List<String> alleles = new ArrayList<String>();
                alleles.add(species);
                species2alleles.put(species, alleles);
            }
        }
		for(int i=0; i<numGTs; i++){
            if(_printDetails){
                System.out.println("\n\nConstructing gt#" + (i+1) +" ...");
            }

            STITree gt = new STITree();
            STINode root = gt.getRoot();

            Map<Tuple<NetNode,NetNode>, List<TNode>> netEdge2geneLineages = new HashMap<Tuple<NetNode,NetNode>, List<TNode>>();
            for(Object nodeObject: Networks.postTraversal(network)){
                NetNode node = (NetNode)nodeObject;
                if(_printDetails){
                    System.out.println("\nNetNode " + node.getName());
                }
                List<TNode> geneLineages = new ArrayList<TNode>();
                if(node.isLeaf()){
                    for(String allele: species2alleles.get(node.getName())){
                        TNode newNode = root.createChild(allele);
                        geneLineages.add(newNode);
                    }
                }
                else{
                    for(Object childObject: node.getChildren()){
                        NetNode child = (NetNode)childObject;
                        geneLineages.addAll(netEdge2geneLineages.get(new Tuple<NetNode, NetNode>(node, child)));
                    }
                }

                if(_printDetails){
                    System.out.println(geneLineages);
                }

                Map<Integer, double[]> gijMap = _node2gij.get(node);
                if(gijMap == null){
                    gijMap = new HashMap<Integer, double[]>();
                    _node2gij.put(node, gijMap);
                }

                if(node.isRoot()){
                    randomlyCoalGeneLineages(geneLineages, geneLineages.size(), 2, root);
                }
                else if(node.isTreeNode()){
                    NetNode parent = (NetNode)node.getParents().iterator().next();
                    double distance = node.getParentDistance(parent);
                    randomlyCoalGeneLineages(geneLineages, gijMap, distance, root);
                    netEdge2geneLineages.put(new Tuple<NetNode, NetNode>(parent, node), geneLineages);
                }
                else{
                    NetNode parent1 = (NetNode)node.getParents().iterator().next();
                    double inheritanceProb = node.getParentProbability(parent1);
                    if(inheritanceProb == NetNode.NO_PROBABILITY){
                        throw new RuntimeException("The network has network node that doesn't have inheritance probability.");

                    }
                    List<TNode> geneLineages1 = new ArrayList<TNode>();
                    List<TNode> geneLineages2 = new ArrayList<TNode>();
                    for(TNode gl: geneLineages){
                        double random = _random.nextDouble();
                        if(random < inheritanceProb){
                            geneLineages1.add(gl);
                        }
                        else{
                            geneLineages2.add(gl);
                        }
                    }
                    if(_printDetails){
                        System.out.println("Dividing into :" + geneLineages1);
                        System.out.println("Dividing into :" + geneLineages2);
                    }

                    int index = 0;
                    for(Object parentObject: node.getParents()){
                        NetNode parent = (NetNode)parentObject;
                        double distance = node.getParentDistance(parent);
                        if(index == 0){
                            geneLineages = geneLineages1;
                        }
                        else{
                            geneLineages = geneLineages2;
                        }
                        randomlyCoalGeneLineages(geneLineages, gijMap, distance, root);
                        netEdge2geneLineages.put(new Tuple<NetNode, NetNode>(parent, node), geneLineages);
                        index++;
                    }
                }

            }
            Trees.removeBinaryNodes(gt);
            gts.add(gt);
        }

        return gts;
	}


    /**
     * This function is to randomly coalesce lineages on a branch
     *
     * @param geneLineages     gene lineages entering a branch
     * @param gijMap           pre-computed gij values
     * @param length           length of the branch
     * @param root             the root of the gene tree
     */
    private void randomlyCoalGeneLineages(List<TNode> geneLineages, Map<Integer, double[]> gijMap, double length, STINode root){
        if(geneLineages.size()<2) return;
        int glInNum = geneLineages.size();
        double[] gijCache = gijMap.get(glInNum);
        if(gijCache == null){
            gijCache = new double[glInNum];
            calculateGij(gijCache, length);
            gijMap.put(glInNum, gijCache);
        }
        int glOutNum = getRandomNumber(gijCache);
        randomlyCoalGeneLineages(geneLineages, glInNum, glOutNum, root);
    }



    /**
     * This function is to randomly coalesce lineages on a branch with a specified number of coalescent events
     *
     * @param geneLineages     gene lineages entering a branch
     * @param glInNum          the number of gene lineages entering the branch
     * @param glOutNum         the number of gene lineages leaving the branch
     * @param root             the root of the gene tree
     */
    private void randomlyCoalGeneLineages(List<TNode> geneLineages, int glInNum, int glOutNum, STINode root){
        if(glInNum < glOutNum){
            return;
        }
        for(int k=glInNum; k>glOutNum; k--){
            int coal1 = (int)(_random.nextDouble()*k);
            int coal2 = coal1;
            while(coal1==coal2){
                coal2 = (int)(_random.nextDouble()*k);
            }
            STINode newNode = root.createChild();
            if(coal1 > coal2){
                newNode.adoptChild((TMutableNode)geneLineages.remove(coal1));
                newNode.adoptChild((TMutableNode)geneLineages.remove(coal2));
            }
            else{
                newNode.adoptChild((TMutableNode)geneLineages.remove(coal2));
                newNode.adoptChild((TMutableNode)geneLineages.remove(coal1));
            }
            geneLineages.add(newNode);

        }
    }



    /**
     * This function is to get random number of lineages leaving a branch according to their probabilities
     */
    private int getRandomNumber(double[] cache){
        double random = _random.nextDouble();
        int coalNum = -1;
        for(int i=0; i<cache.length; i++){
            if(random<cache[i]){
                coalNum = i+1;
                break;
            }
        }
        return coalNum;
    }



    /**
     * This function is to compute the gij value for a branch length
     *
     * @param cache     the resulting gij values
     * @param length    the branch length
     */
    private void calculateGij(double[] cache, double length){
        if(length == NetNode.NO_DISTANCE){
            throw new RuntimeException("The network has branch that doesn't have branch length.");
        }
        int i = cache.length;
        double total = 0;
        for(int j=1; j<=i; j++){
            total += gij(length, i, j);
            cache[j-1] = total;
        }
    }


    /**
     * The function is to calculate the g_{ij} function.
     *
     * @param	length	the branch length
     * @param 	i	the number of lineages in
     * @param	j	the number of lineages out
     *
     * @return	the resulting probability
     */
    private double gij(double length, int i, int j){
        if(length==0){
            if(i==j)
                return 1;
            else
                return 0;
        }
        if(i==0){
            return 1;
        }

        double result = 0;
        for(int k=j; k<=i; k++){
            double temp = Math.exp(0.5 * k * (1.0 - k) * length) * (2.0 * k - 1.0) * Math.pow(-1, k - j)*(fact(j, j + k - 2))*(fact(i - k + 1, i));
            double denom = fact(1, j)*(fact(1, k - j))*(fact(i, i + k - 1));
            result += temp/denom;
        }
        return result;
    }


    /**
     * The function is to calculate factorial
     * @param	start	the first number
     * @param 	end		the last number
     *
     * @return	the resulting factorial
     */
    private double fact(int start, int end){
        double result = 1;
        for(int i=start; i<=end; i++){
            result = result*i;
        }

        return result;
    }



}
