package edu.rice.cs.bioinfo.programs.phylonet.algos.generatesequences;

import org.apache.commons.io.IOUtils;

import java.io.*;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class GenerateSeqGenData
{

    Map<String,String> alleleMap;

    public GenerateSeqGenData(List<String> trees, int loci, double mutationRate)
    {
        alleleMap = generateSequences(trees, loci, mutationRate);
    }

    private static Map<String,String> generateSequences(List<String> trees, int loci, double mutationRate)
    {
        try
        {
            List<String> resultLines = getRawFastaDataFromSeqgen(trees, loci, mutationRate);

            Map<String,String> dnaMap = new HashMap<>();
            for (int i = 1; i < resultLines.size() ;i++)
            {
                String currentLine = resultLines.get(i);

                String[] parts =  currentLine.split("\\s+");
                String dna = parts[1];
                String number = parts[0];

                dnaMap.put(number,dna);
            }

            return dnaMap;

        } catch (Exception e)
        {
            throw new RuntimeException(e);
        }
    }

    private static List<String> getRawFastaDataFromSeqgen(List<String> trees, int loci, double mutationRate) throws IOException, InterruptedException
    {
        Path tempFile = Files.createTempFile("tempSeqGenTrees", "tempFile");
        tempFile.toFile().deleteOnExit();

        try (BufferedWriter writer = Files.newBufferedWriter(tempFile, Charset.defaultCharset()))
        {
            IOUtils.writeLines(trees, "\n", writer);
        }

        Path tempError = Files.createTempFile("tempSeqGenTreesError", "tempFileError");
        tempError.toFile().deleteOnExit();


        ProcessBuilder pb = new ProcessBuilder("seq-gen.exe","-mHKY","-p",Integer.toString(loci),"-l",
                Integer.toString(loci),"-s",Double.toString(mutationRate));

        pb.redirectError(ProcessBuilder.Redirect.INHERIT);
        pb.redirectInput(tempFile.toFile());

        Process p = pb.start();

        List<String> resultLines = IOUtils.readLines(p.getInputStream());
        p.getInputStream().close();


        p.waitFor();
        return resultLines;
    }

    public Map<String,String> getAlleleToDNAMap()
    {
        return alleleMap;
    }

    public static void main(String[] args)
    {
        Map<String,Integer> counts = new HashMap<>();
        counts.put("1",1);
        counts.put("2",1);
        counts.put("3",1);

        String net = "((1:3,ANC#H1:2.25)AA:2,((2:.75)ANC#H1:0,3:.75)CA:4.25)ROOT;";
        int numberOfLoci = 10000;
        GenerateMSGeneTrees foo = new GenerateMSGeneTrees(net,counts,10000);

        GenerateSeqGenData blah = new GenerateSeqGenData(foo.getGeneTrees(),10000,.2);
    }


}
