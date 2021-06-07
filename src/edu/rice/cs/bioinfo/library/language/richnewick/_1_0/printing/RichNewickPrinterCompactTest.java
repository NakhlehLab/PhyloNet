package edu.rice.cs.bioinfo.library.language.richnewick._1_0.printing;

import edu.rice.cs.bioinfo.library.programming.Func1;
import edu.rice.cs.bioinfo.library.programming.Func2;
import junit.framework.Assert;
import org.junit.Test;

import java.io.IOException;
import java.io.StringWriter;
import java.util.HashMap;
import java.util.LinkedList;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 5/14/12
 * Time: 5:03 PM
 * To change this template use File | Settings | File Templates.
 */
public class RichNewickPrinterCompactTest
{
    @Test
    public void test1() throws IOException
    {
        LinkedList<String> rDesc = new LinkedList<String>();
        rDesc.add("X");
        rDesc.add("I");

        LinkedList<String> xDesc = new LinkedList<String>();
        xDesc.add("A");
        xDesc.add("I");

        LinkedList<String> iDesc = new LinkedList<String>();
        iDesc.add("C");
        iDesc.add("B");

        LinkedList<String> aDesc = new LinkedList<String>();
        LinkedList<String> bDesc = new LinkedList<String>();
        LinkedList<String> cDesc = new LinkedList<String>();

        final HashMap<String, Iterable<String>> destinations = new HashMap<String, Iterable<String>>();
        destinations.put("R", rDesc);
        destinations.put("I", iDesc);
        destinations.put("A", aDesc);
        destinations.put("B", bDesc);
        destinations.put("C", cDesc);
        destinations.put("X", xDesc);

        Func1<String, String> getLabel = new Func1<String, String>() {
            public String execute(String input) {
                return input;
            }
        };

        Func1<String, Iterable<String>> nodeToDest = new Func1<String, Iterable<String>>()
        {
            public Iterable<String> execute(String input)
            {
                return destinations.get(input);
            }
        };

        Func2<String,String,String> getBranchLength = new Func2<String, String, String>() {
            public String execute(String input1, String input2) {
                return input1 == "R" && input2 == "X" ? "2" : "1";
            }
        };

         Func2<String,String,String> getSupport = new Func2<String, String, String>() {
            public String execute(String input1, String input2) {
                return input1 == "R" && input2 == "A" ? ".3" : ".5";
            }
        };

        Func2<String,String,String> getProbability = new Func2<String, String, String>() {
            public String execute(String input1, String input2) {
                return input1 == "I" && input2 == "C" ? ".2" : ".6";
            }
        };

        Func1<String,HybridNodeType> getHybridType = new Func1<String,HybridNodeType>()
        {
            public HybridNodeType execute(String input)
            {
                return input == "I" ? HybridNodeType.Hybridization : null;
            }
        };

        Func1<String,String> getHybridNodeIndex = new Func1<String, String>() {
            public String execute(String input) {
                return input == "I" ? "1" : null;
            }
        };

        StringWriter sw = new StringWriter();
        RichNewickPrinterCompact printer =  new RichNewickPrinterCompact<String>();
        printer.setGetBranchLength(getBranchLength);
        printer.setGetSupport(getSupport);
        printer.setGetProbability(getProbability);
        printer.print(true, "R", getLabel, nodeToDest, getHybridNodeIndex, getHybridType, sw);
        sw.flush();
        sw.close();
        String result = sw.toString();
        Assert.assertEquals("(I#H1:1:.5:.6,((B:1:.5:.6,C:1:.5:.2)I#H1:1:.5:.6,A:1:.5:.6)X:2:.5:.6)R;", result);





    }

    @Test
    public void test2() throws IOException
    {
        LinkedList<String> rDesc = new LinkedList<String>();
        rDesc.add("A");
        rDesc.add("I");

        LinkedList<String> iDesc = new LinkedList<String>();
        iDesc.add("C");
        iDesc.add("B");

        LinkedList<String> aDesc = new LinkedList<String>();
        LinkedList<String> bDesc = new LinkedList<String>();
        LinkedList<String> cDesc = new LinkedList<String>();

        final HashMap<String, Iterable<String>> destinations = new HashMap<String, Iterable<String>>();
        destinations.put("R", rDesc);
        destinations.put("I", iDesc);
        destinations.put("A", aDesc);
        destinations.put("B", bDesc);
        destinations.put("C", cDesc);

        Func1<String, String> getLabel = new Func1<String, String>() {
            public String execute(String input) {
                return input;
            }
        };

        Func1<String, Iterable<String>> nodeToDest = new Func1<String, Iterable<String>>()
        {
            public Iterable<String> execute(String input)
            {
                return destinations.get(input);
            }
        };

        Func2<String,String,String> getBranchLength = new Func2<String, String, String>() {
            public String execute(String input1, String input2) {
                return input1 == "R" && input2 == "A" ? "2" : "1";
            }
        };

         Func2<String,String,String> getSupport = new Func2<String, String, String>() {
            public String execute(String input1, String input2) {
                return input1 == "R" && input2 == "A" ? ".3" : ".5";
            }
        };

        Func2<String,String,String> getProbability = new Func2<String, String, String>() {
            public String execute(String input1, String input2) {
                return input1 == "R" && input2 == "A" ? ".2" : ".6";
            }
        };

        Func1<String,HybridNodeType> getHybridType = new Func1<String,HybridNodeType>()
        {
            public HybridNodeType execute(String input)
            {
                return null;
            }
        };

        Func1<String,String> getHybridNodeIndex = new Func1<String, String>() {
            public String execute(String input) {
                return null;
            }
        };

        StringWriter sw = new StringWriter();
        RichNewickPrinterCompact printer = new RichNewickPrinterCompact<String>();
        printer.setGetBranchLength(getBranchLength);
        printer.setGetProbability(getProbability);
        printer.setGetSupport(getSupport);
        printer.print(false, "R", getLabel, nodeToDest, getHybridNodeIndex, getHybridType, sw);
        sw.flush();
        sw.close();
        String result = sw.toString();
        Assert.assertEquals("[&U]((B:1:.5:.6,C:1:.5:.6)I:1:.5:.6,A:2:.3:.2)R;", result);





    }

    @Test
    public void test3() throws IOException
    {
        LinkedList<String> rDesc = new LinkedList<String>();
        rDesc.add("A");
        rDesc.add("I");

        LinkedList<String> iDesc = new LinkedList<String>();
        iDesc.add("R");
        iDesc.add("C");
        iDesc.add("B");

        LinkedList<String> aDesc = new LinkedList<String>();
        aDesc.add("R");
        LinkedList<String> bDesc = new LinkedList<String>();
        bDesc.add("I");
        LinkedList<String> cDesc = new LinkedList<String>();
        cDesc.add("I");

        final HashMap<String, Iterable<String>> destinations = new HashMap<String, Iterable<String>>();
        destinations.put("R", rDesc);
        destinations.put("I", iDesc);
        destinations.put("A", aDesc);
        destinations.put("B", bDesc);
        destinations.put("C", cDesc);

        Func1<String, String> getLabel = new Func1<String, String>() {
            public String execute(String input) {
                return input;
            }
        };

        Func1<String, Iterable<String>> nodeToDest = new Func1<String, Iterable<String>>()
        {
            public Iterable<String> execute(String input)
            {
                return destinations.get(input);
            }
        };

        Func2<String,String,String> getBranchLength = new Func2<String, String, String>() {
            public String execute(String input1, String input2) {
                return null;
            }
        };

         Func2<String,String,String> getSupport = new Func2<String, String, String>() {
            public String execute(String input1, String input2) {
                 return null;
            }
        };

        Func2<String,String,String> getProbability = new Func2<String, String, String>() {
            public String execute(String input1, String input2) {
                 return null;
            }
        };

        Func1<String,HybridNodeType> getHybridType = new Func1<String,HybridNodeType>()
        {
            public HybridNodeType execute(String input)
            {
                return null;
            }
        };

        Func1<String,String> getHybridNodeIndex = new Func1<String, String>() {
            public String execute(String input) {
                return null;
            }
        };

        StringWriter sw = new StringWriter();
        RichNewickPrinterCompact printer = new RichNewickPrinterCompact<String>();
        printer.setGetBranchLength(getBranchLength);
        printer.setGetProbability(getProbability);
        printer.setGetSupport(getSupport);
        printer.print(false, "A", getLabel, nodeToDest, getHybridNodeIndex, getHybridType, sw);
        sw.flush();
        sw.close();
        String result = sw.toString();
        Assert.assertEquals("[&U](((B,C)I)R,A);", result);





    }
}
