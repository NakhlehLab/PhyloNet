package edu.rice.cs.bioinfo.programs.phylonet.algos.generatesequences;

import org.apache.commons.io.IOUtils;

import java.util.List;

public class RunCommandUtilities
{
    public static List<String> runCommandOnArguments(String[] args)
    {
        try
        {
            ProcessBuilder builder = new ProcessBuilder(args);
            builder.redirectError(ProcessBuilder.Redirect.INHERIT);
            Process p = builder.start();
            p.getOutputStream().close();


            List<String> lines = IOUtils.readLines(p.getInputStream());
            p.waitFor();
            return lines;
        }
        catch (Exception e)
        {
            throw new RuntimeException(e);
        }
    }
}
