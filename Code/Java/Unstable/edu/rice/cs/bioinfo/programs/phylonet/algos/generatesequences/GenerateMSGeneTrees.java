package edu.rice.cs.bioinfo.programs.phylonet.algos.generatesequences;

import java.util.*;

public class GenerateMSGeneTrees
{
    private final EthanGenerateMSCommand generateMSCommand;
    private final List<String> geneTrees;

    public GenerateMSGeneTrees(String networkString, Map<String,Integer> speciesToNumberOfAllelesMap, int numberOfLoci)
    {
        generateMSCommand = new EthanGenerateMSCommand(networkString,speciesToNumberOfAllelesMap,numberOfLoci);
        geneTrees = computeGeneTrees("ms.exe "+generateMSCommand.getCommand());
    }

    public List<String> getGeneTrees()
    {
        return geneTrees;
    }

    private static List<String> computeGeneTrees(String command)
    {
        System.out.println(command);

        String[] args = command.split("\\s+");

        List<String> lines = RunCommandUtilities.runCommandOnArguments(args);

        List<String> actualTrees = lines.subList(4,lines.size());

        return actualTrees;
    }



    public Map<String,List<String>> getSpeciesToAlleleMap()
    {
        return generateMSCommand.getSeciesToAlleleMap();
    }

    public static void main(String[] args)
    {
        //String net = "(((1:.75)ANC1#H1:0)ANC2#H2:7.25,(ANC2#H2:3,(ANC1#H1:0,2:.75)Earliest:3)Medium:4.25)Root;";
        String net = "((1:3,ANC#H1:2.25)AA:2,((2:.75)ANC#H1:0,3:.75)CA:4.25)ROOT;";

        Map<String,Integer> mp = new HashMap<String,Integer>();
        mp.put("1",1);
        mp.put("2",1);
        mp.put("3",1);

        List<String> trees = new GenerateMSGeneTrees(net,mp,10000).getGeneTrees();
        System.out.println(trees);
    }
}
