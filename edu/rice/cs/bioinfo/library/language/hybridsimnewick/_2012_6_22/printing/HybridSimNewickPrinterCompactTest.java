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

package edu.rice.cs.bioinfo.library.language.hybridsimnewick._2012_6_22.printing;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func2;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import junit.framework.Assert;
import org.junit.Test;

import java.io.StringWriter;
import java.util.Arrays;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/23/12
 * Time: 3:03 PM
 * To change this template use File | Settings | File Templates.
 */
public class HybridSimNewickPrinterCompactTest
{
    @Test
    public void testPrint1()
    {

        Func1<String, String> getLabel = new Func1<String, String>() {
            public String execute(String input) {
                return input;
            }
        };

        Func2<String,String,String> getBranchLength = new Func2<String,String,String>()
        {
            public String execute(String tail, String tip) {

                if(tail.equals("R") && tip.equals("J"))
                    return "1";
                if(tail.equals("R") && tip.equals("K"))
                    return "2";
                if(tail.equals("J") && tip.equals("X"))
                    return "3";
                if(tail.equals("K") && tip.equals("X"))
                    return "4";
                if(tail.equals("J") && tip.equals("A"))
                    return "5";
                if(tail.equals("K") && tip.equals("B"))
                    return "6";

                throw new RuntimeException();

            }
        };

        Func2<String,String,String> getProbability = new Func2<String, String, String>() {
            public String execute(String tail, String tip) {

                if(tail.equals("J") && tip.equals("X"))
                    return ".7";

                if(tail.equals("K") && tip.equals("X"))
                    return ".3";

                throw new RuntimeException();

            }
        };

        Func1<String, String> getHybridIndex = new Func1<String, String>() {
            public String execute(String node) {
                return node.equals("X") ? "1" : null;
            }
        };

        Func1<String, Iterable<String>> getDestinationNodes = new Func1<String, Iterable<String>>() {
            public Iterable<String> execute(String node) {
                if(node.equals("R"))
                    return Arrays.asList("J", "K");
                if(node.equals("J"))
                    return Arrays.asList("X", "A");
                if(node.equals("K"))
                    return Arrays.asList("X", "B");
                if(node.equals("X"))
                    return Arrays.asList();
                if(node.equals("A"))
                    return Arrays.asList();
                if(node.equals("B"))
                    return Arrays.asList();

                throw new RuntimeException();
            }
        };

        Func1<String, Tuple<String,String>> getHybridParents = new Func1<String, Tuple<String, String>>() {
            public Tuple<String, String> execute(String node) {

                if(node.equals("X"))
                    return new Tuple<String, String>("K", "J");

                throw new RuntimeException();
            }
        };

        StringWriter writer = new StringWriter();
        new HybridSimNewickPrinterCompact(getLabel, getBranchLength, getProbability, getHybridIndex).print("R", getDestinationNodes, getHybridParents, "500", writer);
        String result = writer.toString();

        Assert.assertEquals("((B:6,X#.3:4)K:2,(A:5,X#.3:3)J:1)R:500;", result);
    }

     @Test
    public void testPrint2()
    {

        Func1<String, String> getLabel = new Func1<String, String>() {
            public String execute(String input) {
                return input;
            }
        };

        Func2<String,String,String> getBranchLength = new Func2<String,String,String>()
        {
            public String execute(String tail, String tip) {

                if(tail.equals("R") && tip.equals("J"))
                    return "1";
                if(tail.equals("R") && tip.equals("K"))
                    return "2";
                if(tail.equals("J") && tip.equals("X"))
                    return "3";
                if(tail.equals("K") && tip.equals("X"))
                    return "4";
                if(tail.equals("J") && tip.equals("A"))
                    return "5";
                if(tail.equals("K") && tip.equals("B"))
                    return "6";
                if(tail.equals("X") && tip.equals("M"))
                    return "7";
                if(tail.equals("M") && tip.equals("N"))
                    return "8";
                if(tail.equals("M") && tip.equals("P"))
                    return "9";
                if(tail.equals("N") && tip.equals("Z"))
                    return "10";
                if(tail.equals("P") && tip.equals("Z"))
                    return "11";


                throw new RuntimeException();

            }
        };

        Func2<String,String,String> getProbability = new Func2<String, String, String>() {
            public String execute(String tail, String tip) {

                if(tail.equals("J") && tip.equals("X"))
                    return ".7";

                if(tail.equals("K") && tip.equals("X"))
                    return ".3";

                if(tail.equals("N") && tip.equals("Z"))
                    return ".4";

                if(tail.equals("P") && tip.equals("Z"))
                    return ".6";

                throw new RuntimeException();

            }
        };

        Func1<String, String> getHybridIndex = new Func1<String, String>() {
            public String execute(String node) {
                if(node.equals("X"))
                    return "1";

                if(node.equals("Z"))
                    return "2";

                return null;
            }
        };

        Func1<String, Iterable<String>> getDestinationNodes = new Func1<String, Iterable<String>>() {
            public Iterable<String> execute(String node) {
                if(node.equals("R"))
                    return Arrays.asList("J", "K");
                if(node.equals("J"))
                    return Arrays.asList("X", "A");
                if(node.equals("K"))
                    return Arrays.asList("X", "B");
                if(node.equals("X"))
                    return Arrays.asList("M");
                if(node.equals("A"))
                    return Arrays.asList();
                if(node.equals("B"))
                    return Arrays.asList();
                if(node.equals("M"))
                    return Arrays.asList("N", "P");
                if(node.equals("N"))
                    return Arrays.asList("Z");
                if(node.equals("P"))
                    return Arrays.asList("Z");
                if(node.equals("Z"))
                    return Arrays.asList();

                throw new RuntimeException();
            }
        };

        Func1<String, Tuple<String,String>> getHybridParents = new Func1<String, Tuple<String, String>>() {
            public Tuple<String, String> execute(String node) {

                if(node.equals("X"))
                    return new Tuple<String, String>("K", "J");

                if(node.equals("Z"))
                    return new Tuple<String, String>("N", "P");

                throw new RuntimeException();
            }
        };

        StringWriter writer = new StringWriter();
        new HybridSimNewickPrinterCompact(getLabel, getBranchLength, getProbability, getHybridIndex).print("R", getDestinationNodes, getHybridParents, "500", writer);
        String result = writer.toString();

        Assert.assertEquals("((B:6,(((Z#.6:11)P:9,(Z#.6:10)N:8)M:7)X#.3:4)K:2,(A:5,X#.3:3)J:1)R:500;", result);
    }
}
