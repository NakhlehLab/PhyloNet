package edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan;

import com.google.gson.Gson;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.model.ProcessOutput;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.model.ProcessedOutputInformation;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.model.RawOutputInformation;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.util.JsonUtilities;

import java.nio.file.Path;
import java.nio.file.Paths;

public class ReprocessOutput
{

    public static void main(String[] args)
    {
        String directory = args[0];
        double threshold = Double.parseDouble(args[1]);

        Path pathToDirectory = Paths.get(directory);

        if (!pathToDirectory.toFile().isDirectory())
            throw new RuntimeException("That path is not a directory");

        Gson g = new Gson();

        Path pathToRawData = pathToDirectory.resolve("rawOutput.json");
        RawOutputInformation info =  JsonUtilities.readJson(pathToRawData, g, RawOutputInformation.class);

        ProcessedOutputInformation process = new ProcessOutput().getProcessedOutput(info,threshold);

        Path pathToProcessed = pathToDirectory.resolve("output.json");
        JsonUtilities.writeJson(pathToProcessed, g, process);

    }

}
