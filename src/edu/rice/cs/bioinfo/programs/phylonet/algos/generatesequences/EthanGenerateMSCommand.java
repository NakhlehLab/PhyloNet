package edu.rice.cs.bioinfo.programs.phylonet.algos.generatesequences;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.HmmNetworkUtils;

import java.util.*;

public class EthanGenerateMSCommand
{

    private String command = "";
    private Map<Integer,String> numberToSpeciesMapping = new HashMap<Integer,String>();
    private Map<String,List<String>> speciesToAlleleMap = new HashMap<String,List<String>>();

    private static class MSData
    {
        final double height;

        MSData(double height)
        {
            this.height = height;
        }

        List<Integer> speciesArriving = new ArrayList<Integer>();

        void addSpecies(int species)
        {
            speciesArriving.add(species);
        }

        boolean processed = false;

        void setProcessed()
        {
            processed = true;
        }
    }


    private List<NetNode<MSData>> getSortedNetNodeAndTimes(Network<MSData> net)
    {

        List<NetNode<MSData>> result = new ArrayList<NetNode<MSData>>();
        for (NetNode<MSData> node : net.bfs())
        {
            result.add(node);
        }


        Collections.sort(result, new Comparator<NetNode<MSData>>()
        {
            @Override
            public int compare(NetNode<MSData> o1, NetNode<MSData> o2)
            {
                return Double.compare(o1.getData().height, o2.getData().height);
            }
        });

        return result;
    }

    public EthanGenerateMSCommand(String network, Map<String, Integer> speciesCounts, int numeberOfLoci)
    {

        int totalCount = 0;
        for (String species: speciesCounts.keySet())
        {
            int currentSize = speciesCounts.get(species);
            totalCount += currentSize;
        }

        command = totalCount + " 1 -r 100 " + numeberOfLoci + " -T -I " + speciesCounts.keySet().size() + " ";

        Network<MSData> net = HmmNetworkUtils.fromENewickString(network);

        prepareNetworkForProcessing(net);

        int currentCount = 1;
        for (int i = 1; i <= speciesCounts.size();i++)
        {
            String species = numberToSpeciesMapping.get(i);
            speciesToAlleleMap.put(species, new ArrayList<String>());
            for (int j = 0 ; j < speciesCounts.get(species);j++)
                speciesToAlleleMap.get(species).add(Integer.toString(currentCount++));
            command+= speciesCounts.get(species) + " ";
        }


        List<NetNode<MSData>> nodes = getSortedNetNodeAndTimes(net);


        for (NetNode<MSData> node : nodes)
        {
            processNode(node);
        }

    }

    private void prepareNetworkForProcessing(Network<MSData> net)
    {
        HmmNetworkUtils.annotateHeights(net);

        numberToSpeciesMapping = new HashMap<Integer, String>();

        for (NetNode<MSData> node: net.bfs())
        {
            NetNode fake = (NetNode) node;
            node.setData(new MSData((Double)fake.getData()));
        }

        for (NetNode<MSData> leaf : net.getLeaves())
        {
            processLeaf(leaf);
        }
    }

    private void processLeaf(NetNode<MSData> node)
    {
        if (node.isLeaf())
        {
            if (node.getIndeg() != 1 || node.getOutdeg() != 0)
                throw new RuntimeException("Not a valid child");

            node.getParents().iterator().next().getData().addSpecies(nextAlleleId);
            numberToSpeciesMapping.put(nextAlleleId, node.getName());
            nextAlleleId+=1;
        }
        else
            throw new RuntimeException("Should be a leaf");

        node.getData().setProcessed();
    }

    private int nextAlleleId = 1;
    private void processNode(NetNode<MSData> node)
    {
        if (node.getData().processed)
            return;

        for (NetNode<MSData> child : node.getChildren())
        {
            processNode(child);
        }

        if (node.getData().speciesArriving.size() != node.getOutdeg())
            throw new RuntimeException("Lying node");
        else if (node.isTreeNode())
        {
            if ((node.getIndeg() != 1 && node.getIndeg() != 0))
                throw new RuntimeException("Not a valid treenode");

            int resultingSpecies;

            if (node.getOutdeg() == 1)
                resultingSpecies = node.getData().speciesArriving.get(0);
            else if (node.getOutdeg() == 2)
            {
                resultingSpecies = Collections.max(node.getData().speciesArriving);
                int badSpecies = Collections.min(node.getData().speciesArriving);
                String commandPart = String.format("-ej %f %d %d ", node.getData().height, badSpecies, resultingSpecies);

                command += commandPart;
                System.out.println(commandPart);
            }
            else
            {
                throw new RuntimeException("Does not support more than 2 children");
            }

            if (node.getIndeg() != 0)
                node.getParents().iterator().next().getData().addSpecies(resultingSpecies);
        }
        else if (node.isNetworkNode())
        {
            if (node.getIndeg() != 2 || node.getOutdeg() != 1)
                throw new RuntimeException("Not a valid network node");

            if (node.getData().speciesArriving.size() != 1)
                throw new RuntimeException("Should only be one child of network");

            Iterator<NetNode<MSData>> parents = node.getParents().iterator();

            NetNode<MSData> firstParent = parents.next();
            int currentSpecies = node.getData().speciesArriving.get(0);
            firstParent.getData().addSpecies(currentSpecies);

            NetNode<MSData> secondParent = parents.next();
            int newSpecies = nextAlleleId;
            nextAlleleId+=1;
            secondParent.getData().addSpecies(newSpecies);

            String commandPart = String.format("-es %f %d %f ", node.getData().height, currentSpecies,node.getParentProbability(secondParent));

            command += commandPart;
            System.out.println(commandPart  + "making " + newSpecies);
        }
        else
            throw new RuntimeException("Invalid type of node?");

        node.getData().setProcessed();

    }

    public String getCommand()
    {
        return command;
    }

    public Map<String,List<String>> getSeciesToAlleleMap()
    {
        return speciesToAlleleMap;
    }


    public static void main(String[] args)
    {
        String net = "((A:3.1,(B:3)ANC#H1:.1):.9,ANC#H1:1);";

        Map<String,Integer> speciesToAlleles = new HashMap<String,Integer>();
        speciesToAlleles.put("A", 1);
        speciesToAlleles.put("B", 1);

        EthanGenerateMSCommand a = new EthanGenerateMSCommand(net,speciesToAlleles,10000);
        System.out.println(a.getCommand());
        System.out.println(a.getSeciesToAlleleMap());
    }

}
