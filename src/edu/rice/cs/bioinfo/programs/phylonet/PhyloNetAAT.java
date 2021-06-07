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

package edu.rice.cs.bioinfo.programs.phylonet;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;
import org.apache.commons.io.FileUtils;
import org.junit.Assert;

import java.io.*;
import java.math.BigDecimal;
import java.util.Random;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/23/11
 * Time: 2:31 PM
 * To change this template use File | Settings | File Templates.
 */

public class PhyloNetAAT  {

    public PhyloNetAAT()
    {

    }

    static class AATTestCase extends TestCase
    {
        private final File _file;

        public AATTestCase(File file)
        {
            _file = file;
            this.setName(_file.getName());
        }

        @Override
        protected void runTest() throws Throwable
        {


            System.out.println("Testing " + _file.getAbsolutePath());
            String fileContents = FileUtils.readFileToString(_file);
            fileContents = fileContents.replace("\r", "");
            String[] parts = fileContents.split(delim);

            if(parts.length == 3)
            {
                String nexus = parts[0];
                String display = parts[1];
                String error = parts[2];
                checkTest(nexus, display, error, _file.getName());
            }
            else if(parts.length == 2)
            {
                String nexus = parts[0];
                String display = parts[1];
                String error = "";
                checkTest(nexus, display, error, _file.getName());
            }
            else
            {
                throw new RuntimeException("Bad AAT script " + _file.toString());
            }
        }
    }


    private static String delim = "\n===\n";

    private static Random _rand = new Random(23);

    public static TestSuite suite() throws Throwable
    {
        TestSuite aatSuite = new TestSuite();

        File currentDir = new File(".");

        for(final File file : currentDir.listFiles(new FilenameFilter() {
            public boolean accept(File dir, String name) {
                return name.endsWith(".txt");
            }
        }))
        {

            new AATTestCase(file).runTest();

        }


        return aatSuite;
    }


    private static void checkTest(String nexus, String expectedStdOut, String expectedStdError, String testFile) throws IOException
    {
        String faultMessage = testFile + " failed.";
        ByteArrayOutputStream display = new ByteArrayOutputStream();
        ByteArrayOutputStream error = new ByteArrayOutputStream();
        Program.run(new ByteArrayInputStream(nexus.getBytes()), new PrintStream(error), new PrintStream(display), _rand, BigDecimal.ZERO);
        Assert.assertEquals(faultMessage, expectedStdError, error.toString().replace("\r", ""));
        Assert.assertEquals(faultMessage,  expectedStdOut, display.toString().replace("\r", ""));
    }
}
