package edu.rice.cs.bioinfo.programs.phylonet;
import org.apache.commons.io.FileUtils;
import org.junit.Test;
import org.junit.Assert;
import sun.net.idn.StringPrep;

import java.awt.geom.Path2D;
import java.io.*;
import java.util.Random;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/23/11
 * Time: 2:31 PM
 * To change this template use File | Settings | File Templates.
 */
public class PhyloNetAAT {


    private static String delim = "\n===\n";

    private static Random _rand = new Random(23);


    @Test
    public void testAATScripts() throws IOException
    {
        File currentDir = new File(".");

        for(File file : currentDir.listFiles(new FilenameFilter() {
            public boolean accept(File dir, String name) {
                return name.endsWith(".txt");
            }
            }))
        {
            System.out.println("Testing " + file.getAbsolutePath());
            String fileContents = FileUtils.readFileToString(file);
            fileContents = fileContents.replace("\r", "");
            String[] parts = fileContents.split(delim);

            if(parts.length == 3)
            {
                String nexus = parts[0];
                String display = parts[1];
                String error = parts[2];
                checkTest(nexus, display, error, file.getName());
            }
            else if(parts.length == 2)
            {
                String nexus = parts[0];
                String display = parts[1];
                String error = "";
                checkTest(nexus, display, error, file.getName());
            }
            else
            {
                throw new RuntimeException("Bad AAT script " + file.toString());
            }


        }
    }

    private void checkTest(String nexus, String expectedStdOut, String expectedStdError, String testFile) throws IOException
    {
        String faultMessage = testFile + " failed.";
        ByteArrayOutputStream display = new ByteArrayOutputStream();
        ByteArrayOutputStream error = new ByteArrayOutputStream();
        Program.run(new ByteArrayInputStream(nexus.getBytes()), new PrintStream(error), new PrintStream(display), _rand);
        Assert.assertEquals(faultMessage, expectedStdError, error.toString().replace("\r", ""));
        Assert.assertEquals(faultMessage,  expectedStdOut, display.toString().replace("\r", ""));
    }
}
