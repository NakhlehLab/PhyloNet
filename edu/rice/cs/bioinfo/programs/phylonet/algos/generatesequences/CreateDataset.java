package edu.rice.cs.bioinfo.programs.phylonet.algos.generatesequences;

import com.google.gson.Gson;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.ConfigurationBuilder;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.util.JsonUtilities;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.HmmNetworkUtils;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

public class CreateDataset
{

    private static CreateDatasetInput loadSettings(Path settingsFile)
    {
        return JsonUtilities.readJson(settingsFile,new Gson(),CreateDatasetInput.class);
    }

    public static void main(String[] args)
    {

        Path settingsFile = Paths.get(args[0]);
        Path finalDirectory = Paths.get(args[1]);

        CreateDatasetInput settings = loadSettings(settingsFile);

        System.out.println(settings);

        System.out.println(HmmNetworkUtils.fromENewickString(settings.network));

        createDataset(finalDirectory, settings);


    }

    public static void createDataset(Path finalDirectory, CreateDatasetInput settings)
    {
        GenerateMSGeneTrees geneTrees = new GenerateMSGeneTrees(settings.network, settings.networkAlleleCounts, settings.numberOfLoci);

        GenerateSeqGenData seqGenData = new GenerateSeqGenData(geneTrees.getGeneTrees(), settings.numberOfLoci, settings.mutationRate);

        ConfigurationBuilder result = ConfigurationBuilder.getSNAPP().withIterations(settings.iterations).withThreads(settings.threads);

        GenerateNexusFile nexus = new GenerateNexusFile(seqGenData.getAlleleToDNAMap(), geneTrees.getSpeciesToAlleleMap(), settings.network, result);

        try
        {
            Files.createDirectories(finalDirectory);

            Path geneTreeFile = finalDirectory.resolve("treefile");
            Path nexusFile = finalDirectory.resolve("nexusFile.nxs");

            Files.write(geneTreeFile, geneTrees.getGeneTrees(), Charset.defaultCharset());
            Files.write(nexusFile, nexus.getNexusFile(), Charset.defaultCharset());

        } catch (IOException e)
        {
            throw new RuntimeException(e);
        }
    }
}
