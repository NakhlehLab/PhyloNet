package edu.rice.cs.bioinfo.programs.phylonet.structs.network.util;

import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.*;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.parsers.antlr.ast.ANTLRRichNewickParser;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReadResult;
import edu.rice.cs.bioinfo.programs.phylonet.commands.NetworkTransformer;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.io.ExNewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TMutableNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.TNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STINode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.util.Trees;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.StringReader;
import java.math.BigDecimal;
import java.util.*;

/**
 * Created by ethan_000 on 7/29/2014.
 */
public class HmmNetworkUtils
{
    public static <E> List<E> makeList(Iterable<E> iter) {
        List<E> list = new ArrayList<E>();
        for (E item : iter) {
            list.add(item);
        }
        return list;
    }

    private static <T> STINode<T> clone(NetNode<T> parentNode, NetNode<T> node, STINode<T> parentTree, Map<String, List<String>> copies)
    {
        STINode<T> child = parentTree.createChildWithUniqueName();
        List<String> list = copies.get(node.getName());
        if (list == null)
        {
            list = new ArrayList<String>();
            copies.put(node.getName(),list);
        }
        list.add(child.getName());
        if (parentNode != null)
            child.setParentDistance(node.getParentDistance(parentNode));

        return child;

    }

    private static <T> STINode<T>  mulTree(NetNode<T> networkParent, NetNode<T> node,STINode<T> treeParent,Map<String,List<String>> copies)
    {
        List<NetNode<T>> children = makeList(node.getChildren());


        if (children.size() == 2)
        {
            STINode<T> current = clone(networkParent, node, treeParent, copies);

            mulTree(node, children.get(0), current, copies);
            mulTree(node, children.get(1), current, copies);

            return current;
        }
        else if (children.size() == 0)
        {
            STINode<T> current = clone(networkParent, node, treeParent,copies);
            return current;
        }
        else if (children.size() == 1)
        {

            STINode<T> chain = mulTree(node,children.get(0),treeParent,copies);

            if (networkParent != null)//happens in the case of this node being the root.
                chain.setParentDistance(chain.getParentDistance()+node.getParentDistance(networkParent));

            return chain;

        }
        else
        {
            throw new RuntimeException("I cannot handle children of other sizes at this point "+node.getName()+","+children.get(0).getName());
        }
    }

    private static <A,B> Map<B,List<A>> inverse(Map<A,B> map)
    {
        Map<B,List<A>> result = new HashMap<B,List<A>>();

        for (Map.Entry<A,B> entry : map.entrySet())
        {
            if (!result.containsKey(entry.getValue()))
                result.put(entry.getValue(),new ArrayList<A>());

            result.get(entry.getValue()).add(entry.getKey());

        }

        return result;
    }

    private static  List<Map<String,String>> possibleAlleleAssignments(String allele, List<String> nodeOptions)
    {
        List<Map<String,String>> result = new ArrayList<Map<String,String>>();
        for (String node : nodeOptions)
        {
            Map<String,String> map = new HashMap<String,String>();
            map.put(allele,node);
            result.add(map);
        }
        return result;
    }

    private static <A,B> List<Map<A,B>> mergeMapCombinations(List<Map<A,B>> firstMaps, List<Map<A,B>> secondMaps)
    {
        List<Map<A,B>> result = new ArrayList<Map<A,B>>();

        for (Map<A,B> firstMap : firstMaps)
            for (Map<A,B> secondMap : secondMaps)
            {
                Map<A, B> combination = new HashMap<A, B>();
                combination.putAll(firstMap);
                combination.putAll(secondMap);
                result.add(combination);
            }

        return result;
    }

    private static List<Map<String,String>> possibleSpeciesAssignments(List<String> alleleOptions, List<String> nodeOptions)
    {
        List<Map<String,String>> start =new ArrayList<Map<String,String>>();
        start.add(new HashMap<String,String>());

        for (String allele : alleleOptions)
        {
            start = mergeMapCombinations(start,possibleAlleleAssignments(allele,nodeOptions));
        }

        return start;
    }

    private static List<Map<String,String>> allPossibleAssignments(Map<String, List<String>> allelesForSpecies, Map<String, List<String>> nodesForAlleles)
    {
        List<Map<String,String>> start =new ArrayList<Map<String,String>>();
        start.add(new HashMap<String,String>());

        for (String species : allelesForSpecies.keySet())
        {
            start = mergeMapCombinations(start,possibleSpeciesAssignments(allelesForSpecies.get(species),nodesForAlleles.get(species)));
        }
        return start;
    }

    private static <T> boolean applyAlleleAssignment(STINode<T> mulTree, Map<String, List<String>> alleleMapping)
    {
        if (mulTree.isLeaf())
        {
            if (alleleMapping.containsKey(mulTree.getName()))
            {
                List<String> names = alleleMapping.get(mulTree.getName());

                for (String name : names)
                {
                    //mulTree.setName(name);
                    STINode<T> child = mulTree.createChild(name);
                    child.setParentDistance(0);
                }


                return true;
            }
            else
                return false;
        }
        else
        {
            List<STINode<T>> children = new ArrayList<STINode<T>>();
            for (STINode<T> child : mulTree.getChildren())
            {
                children.add(child);
            }

            for (STINode<T> child: children)
            {
                if (!applyAlleleAssignment(child,alleleMapping))
                    mulTree.removeChild(child, false);
            }

            if (mulTree.getChildCount() == 1 && mulTree.getParent() != null)
            {
                STINode<T> onlyChild = mulTree.getChildren().iterator().next();
                double totalDist = mulTree.getParentDistance() + onlyChild.getParentDistance();

                mulTree.getParent().removeChild(mulTree,true);
                onlyChild.setParentDistance(totalDist);
                return true;
            }

            return mulTree.getChildCount()>0;
        }

    }

    /**
     * Creates all possible parental trees from a given network.
     * It first creates a MUL Tree, and then applies each possible allele mapping to it.
     * @param network A species network
     * @param speciesOptions A mapping of species names to all the alleles for that species
     * @return A list of all possible parental trees.
     */
    public static <T> List<STITree<T>> generateParentalTrees(Network<T> network, Map<String, List<String>> speciesOptions)
    {

        NetNode<T> root = network.getRoot();

        Map<String,List<String>> alleleNamesToNodeNames = new HashMap<String,List<String>>();
        STITree<T> tree = new STITree<T>(root.getName(),true);
        STINode<T> mulTree = mulTree(null, root, tree.getRoot(), alleleNamesToNodeNames);


        if (!allSpeciesOptionsInNetwork(speciesOptions,alleleNamesToNodeNames))
            throw new RuntimeException("Some of the alleles in the species options are not in the species network.");


        List<Map<String,String>> alleleMappingsOfAlleleNamesToNodeName = allPossibleAssignments(speciesOptions, alleleNamesToNodeNames);

        List<Map<String,List<String>>> alleleAssignmentsOfNodeNamesToAlleleNames = inverseList(alleleMappingsOfAlleleNamesToNodeName);

        List<STITree<T>> results = applyAllAlleleAssignementAndClone(mulTree,alleleAssignmentsOfNodeNamesToAlleleNames);
        //System.out.println(results);
        return results;
    }

    private static <T> STITree<T> applyAlleleAssignementAndClone(STINode<T> mulTree, Map<String,List<String>> alleleAssignment)
    {
        STITree<T> temp = new STITree<T>(mulTree);
        applyAlleleAssignment(temp.getRoot(), alleleAssignment);

        while (temp.getRoot().getChildCount() == 1)
            temp = new STITree<T>(temp.getRoot().getChildren().iterator().next());
        return temp;
    }

    private static <T> List<STITree<T>> applyAllAlleleAssignementAndClone(STINode<T> mulTree, List<Map<String,List<String>>> alleleAssignments)
    {
        List<STITree<T>> result = new ArrayList<STITree<T>>();
        for (Map<String,List<String>> alleleAssignment : alleleAssignments) {
            STITree<T> newTree = applyAlleleAssignementAndClone(mulTree, alleleAssignment);
            Trees.removeBinaryNodes(newTree);
            result.add(newTree);
        }
        return result;
    }

    private static <A,B> List<Map<B, List<A>>> inverseList(List<Map<A, B>> listToInverse)
    {
        List<Map<B, List<A>>> result = new ArrayList<Map<B, List<A>>>();
        for (Map<A,B> item : listToInverse)
            result.add(inverse(item));
        return result;
    }

    private static boolean allSpeciesOptionsInNetwork(Map<String, List<String>> speciesOptions, Map<String, List<String>> alleleNamesToNodeNames)
    {
        for (String species: speciesOptions.keySet())
            if (!alleleNamesToNodeNames.containsKey(species))
                return false;

        return true;
    }

    public static <T> Network<T> fromString(String s) {
        ExNewickReader<T> read = new ExNewickReader<T>(new StringReader(s));

        try {
            return read.readNetwork();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Creates a network from a tree.
     *
     * @param speciesTree The tree to convert.
     * @return The tree in network form. Note the parametrized type (so it can work with the GeneTreeProbability class).
     */
    public static <T> Network<T> toNetwork(STITree<String> speciesTree)
    {
        String newickString = speciesTree.toNewickWD();

        ExNewickReader<T> reader = new ExNewickReader<T>(new StringReader(newickString));

        try
        {
            return reader.readNetwork();
        } catch (IOException e)
        {
            e.printStackTrace();
            throw new RuntimeException(e);
        }
    }

    public static void annotateHeights(Network net)
    {
        Network<Double> clone = (Network<Double>) net;
        for (NetNode<Double> node : clone.bfs())
        {
            node.setData(null);
        }

        annotateHeightsNode(clone.getRoot());
    }

    private static double annotateHeightsNode(NetNode<Double> root)
    {
        if (root.getData() != null)
            return root.getData();

        if (root.isLeaf())
        {
            root.setData(0.0);
            return 0;
        }
        else
        {
            double dist = Double.NEGATIVE_INFINITY;
            for (NetNode<Double> child : root.getChildren())
            {
                double childDist = annotateHeightsNode(child) + child.getParentDistance(root);
                if (dist == Double.NEGATIVE_INFINITY)
                    dist = childDist;

                if (Math.abs(dist - childDist) > .01)
                    throw new RuntimeException("The network passed in is not ultrametric");
            }
            root.setData(dist);
            return dist;
        }
    }

    public static <T> Network<T> fromENewickString(String net)
    {
        return (Network<T>) Networks.readNetwork(net);
    }


    /**
     * The function is to convert a network to a multilabel tree.
     * @param	net 	the given network
     */

    public static STITree generateMulTree(Network net, Map<String,List<String>> species2alleles, Map<NetNode,List<TNode>> netnode2treenodes, List<Map<String, List<String>>> alleleMappings){
        STITree mulTree = new STITree<Double>();
        ((STINode<Double>)(mulTree.getRoot())).setData(1.0);
        ((STINode<Double>)(mulTree.getRoot())).setName("root");
        Queue<NetNode<Double>> source = new LinkedList<NetNode<Double>>();
        Queue<TMutableNode> dest = new LinkedList<TMutableNode>();
        NetNode root = net.getRoot();
        if(root.getName().equals("")){
            root.setName("root");
        }
        source.offer(root);
        dest.offer((TMutableNode) mulTree.getRoot());
        int nameid = 1;
        while(!source.isEmpty()){
            NetNode<Double> parent = source.poll();
            TMutableNode peer = dest.poll();
            for (NetNode<Double> child : parent.getChildren()) {
                if (child.getName().equals(NetNode.NO_NAME)) {
                    child.setName("i" + (nameid++));
                }

                List<TNode> correspondingTreeNodes = netnode2treenodes.get(child);
                if(correspondingTreeNodes == null){
                    correspondingTreeNodes = new ArrayList<>();
                    netnode2treenodes.put(child, correspondingTreeNodes);
                }
                String name = child.getName();
                String newname = name + "#" + (correspondingTreeNodes.size()+1);
                if(child.isNetworkNode()){
                    newname += "~" + parent.getName();
                }
                TMutableNode copy = peer.createChild(newname);
                correspondingTreeNodes.add(copy);

                // Update the distance and data for this child.
                double distance = child.getParentDistance(parent);
                if (distance == NetNode.NO_DISTANCE) {
                    copy.setParentDistance(0);
                }
                else {
                    copy.setParentDistance(distance);
                }

                double gamma = child.getParentProbability(parent);
                gamma = gamma==NetNode.NO_PROBABILITY?1.0:gamma;
                ((STINode<Double>)copy).setData(gamma);

                // Continue to iterate over the children of nn and tn.
                source.offer(child);
                dest.offer(copy);
                //index ++;
            }
        }

        Map<String,List<String>> speciesNamesToNodeNames = new HashMap<String,List<String>>();
        for(Map.Entry<NetNode,List<TNode>> entry: netnode2treenodes.entrySet()){
            String speciesName = "";
            if(entry.getKey().isLeaf()){
                speciesName = entry.getKey().getName();
            }
            List<String> nodeNames = speciesNamesToNodeNames.get(speciesName);
            if(nodeNames == null){
                nodeNames = new ArrayList<>();
                speciesNamesToNodeNames.put(speciesName, nodeNames);
            }
            for(TNode node: entry.getValue()){
                nodeNames.add(node.getName());
            }
        }

        List<Map<String,String>> alleleMappingsOfAlleleNamesToNodeName = allPossibleAssignments(species2alleles, speciesNamesToNodeNames);
        alleleMappings.addAll(inverseList(alleleMappingsOfAlleleNamesToNodeName));

        return mulTree;
    }

}
