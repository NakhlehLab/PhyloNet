package edu.rice.cs.bioinfo.programs.phylonet.algos.generatesequences;

import edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.Configuration;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.ConfigurationBuilder;
import org.apache.commons.io.IOUtils;

import java.io.IOException;
import java.io.InputStream;
import java.util.*;

public class GenerateNexusFile
{
    private final Map<String, String> alleleToDnaMap;
    private final Map<String, List<String>> speciesToAlleleMap;
    private final String network;
    private final Configuration configuration;


    public GenerateNexusFile(Map<String, String> alleleToDnaMap, Map<String, List<String>> speciesToAlleleMap, String network)
    {
        this(alleleToDnaMap, speciesToAlleleMap, network, ConfigurationBuilder.getSNAPP());
    }

    public GenerateNexusFile(Map<String, String> alleleToDnaMap, Map<String, List<String>> speciesToAlleleMap, String network, ConfigurationBuilder configuration)
    {
        this.alleleToDnaMap = alleleToDnaMap;
        this.speciesToAlleleMap = speciesToAlleleMap;
        this.network = network;
        this.configuration = configuration.build("");
    }

    public String getDnaMatrix()
    {
        String result = "";
        for (Map.Entry<String, String> alleleEntry : alleleToDnaMap.entrySet())
        {
            result += alleleEntry.getKey() + " " + alleleEntry.getValue() + "\n";
        }
        return result;
    }

    //From http://stackoverflow.com/a/63201/406009
    public static String join(Iterable<? extends CharSequence> s, String delimiter) {
        Iterator<? extends CharSequence> iter = s.iterator();
        if (!iter.hasNext()) return "";
        StringBuilder buffer = new StringBuilder(iter.next());
        while (iter.hasNext()) buffer.append(delimiter).append(iter.next());
        return buffer.toString();
    }

    public String getAlleleMap()
    {
        String result = "<";

        List<String> entries = new ArrayList<String>();
        for (String species : speciesToAlleleMap.keySet())
            entries.add(getAlleleMapEntry(species));

        result += join(entries,"; ");
        result += ">";
        return result;
    }

    public String getAlleleMapEntry(String species)
    {
        String result = species + ":";
        List<String> alleles = speciesToAlleleMap.get(species);
        result += join(alleles,",");
        return result;
    }

    public String getNexusOptionalParameters(Configuration configuration)
    {
        String result = "-threads " + configuration.threads + " -iterations " + configuration.ITERATIONS;
        if (configuration.PPATraining)
            result += " -ppatraining " + configuration.numberOfFolds;
        return result;
    }

    public List<String> getNexusFile()
    {
        try
        {
            InputStream nexusTemplateResource = getClass().getResourceAsStream("/nexusTemplate.txt");
            String nexusTemplate = IOUtils.toString(nexusTemplateResource);

            int numberOfLoci = alleleToDnaMap.values().iterator().next().length();

            String filledInTemplate = String.format(nexusTemplate, network, alleleToDnaMap.size(),
                    numberOfLoci, getDnaMatrix(), getAlleleMap(),getNexusOptionalParameters(configuration));

            return Arrays.asList(filledInTemplate);
        } catch (IOException e)
        {
           throw new RuntimeException(e);
        }
    }

    public static void main(String[] args)
    {
        String net = "((1:3,ANC#H1:2.25)AA:2,((2:.75)ANC#H1:0,3:.75)CA:4.25)ROOT;";

        Map<String,Integer> mp = new HashMap<String,Integer>();
        mp.put("1",1);
        mp.put("2",1);
        mp.put("3",1);

        GenerateMSGeneTrees ms = new GenerateMSGeneTrees(net,mp,1000);
        GenerateSeqGenData seq = new GenerateSeqGenData(ms.getGeneTrees(),1000,.2);

        List<String> result = new GenerateNexusFile(seq.getAlleleToDNAMap(),ms.getSpeciesToAlleleMap(),net, ConfigurationBuilder.getSNAPP().withIterations(100).withPPATraining(4)).getNexusFile();

        System.out.println(result);
    }
}
