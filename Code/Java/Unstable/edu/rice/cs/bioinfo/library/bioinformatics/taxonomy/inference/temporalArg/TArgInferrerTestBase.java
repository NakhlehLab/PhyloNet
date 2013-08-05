package edu.rice.cs.bioinfo.library.bioinformatics.taxonomy.inference.temporalArg;

import edu.rice.cs.bioinfo.library.language.dot_2013_1.printing.DigraphDotPrinter;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import junit.framework.Assert;
import org.junit.Test;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 7/18/13
 * Time: 1:34 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class TArgInferrerTestBase<T extends TArgInferrer<Tuple<String,Integer>,?>>
{
    public abstract T makeInferrer();

    @Test
    public void testCase1()
    {
        Tuple<String,Integer> sequencing1  = new Tuple<String, Integer>("AAAAA", 0);

        Tuple<String,Integer> sequencing2  = new Tuple<String, Integer>("AATAA", 1);

        Tuple<String,Integer> sequencing3  = new Tuple<String, Integer>("AATAA", 2);

        Tuple<String,Integer> sequencing4  = new Tuple<String, Integer>("AAATA", 3);
        Tuple<String,Integer> sequencing5  = new Tuple<String, Integer>("ATTAA", 3);

        Tuple<String,Integer> sequencing6  = new Tuple<String, Integer>("ATTAA", 4);
        Tuple<String,Integer> sequencing7  = new Tuple<String, Integer>("AATAT", 4);

        Tuple<String,Integer> sequencing8  = new Tuple<String, Integer>("TAATA", 5);

        Tuple<String,Integer> sequencing9  = new Tuple<String, Integer>("AAATA", 6);
        Tuple<String,Integer> sequencing10  = new Tuple<String, Integer>("TAATA", 6);
        Tuple<String,Integer> sequencing11 = new Tuple<String, Integer>("ATATA", 6);
        Tuple<String,Integer> sequencing12 = new Tuple<String, Integer>("ATTAA", 6);
        Tuple<String,Integer> sequencing13 = new Tuple<String, Integer>("ATTAT", 6);
        Tuple<String,Integer> sequencing14 = new Tuple<String, Integer>("ATTAT", 6);
        Tuple<String,Integer> sequencing15 = new Tuple<String, Integer>("AATAT", 6);
        Tuple<String,Integer> sequencing16 = new Tuple<String, Integer>("AATAA", 6);

        T inferrer = makeInferrer();
        inferrer.InferGraph(Arrays.asList(sequencing1, sequencing2, sequencing3,
                                          sequencing4, sequencing5, sequencing6,
                                          sequencing7, sequencing8, sequencing9,
                                          sequencing10, sequencing11, sequencing12,
                                          sequencing11,  sequencing12, sequencing13,
                                          sequencing14, sequencing15, sequencing16));

        // identical
        Assert.assertTrue(isEdge(sequencing4, sequencing9, inferrer));
        Assert.assertTrue(isEdge(sequencing8, sequencing10, inferrer));
        Assert.assertTrue(isEdge(sequencing6, sequencing12, inferrer));
        Assert.assertTrue(isEdge(sequencing7, sequencing15, inferrer));
        Assert.assertTrue(isEdge(sequencing3, sequencing16, inferrer));
        Assert.assertTrue(isEdge(sequencing5, sequencing6, inferrer));
        Assert.assertTrue(isEdge(sequencing2, sequencing3, inferrer));

        // off by 1
        Assert.assertTrue(isEdge(sequencing1, sequencing4, inferrer));
        Assert.assertTrue(isEdge(sequencing1, sequencing2, inferrer));
        Assert.assertTrue(isEdge(sequencing3, sequencing5, inferrer));
        Assert.assertTrue(isEdge(sequencing3, sequencing7, inferrer));
        Assert.assertTrue(isEdge(sequencing4, sequencing8, inferrer));
    }

    @Test
    public void testCaseNIH()
    {
        Tuple<String,Integer> sequencing1_1  = new Tuple<String, Integer>("00000000000000000000000000000000000000000", 615);
        Tuple<String,Integer> sequencing1_2  = new Tuple<String, Integer>("11100000000000000000000000000000000000000", 616);
        Tuple<String,Integer> sequencing1_3  = new Tuple<String, Integer>("00000000000000000000000000000000000000000", 617);
        Tuple<String,Integer> sequencing1_4  = new Tuple<String, Integer>("00000000000000000000000000000000000000000", 619);
        Tuple<String,Integer> sequencing1_5  = new Tuple<String, Integer>("00000000000000000000000000000000000000000", 627);
        Tuple<String,Integer> sequencing1_6  = new Tuple<String, Integer>("11100000000000000000000000000000000000000", 630);
        Tuple<String,Integer> sequencing1_7  = new Tuple<String, Integer>("11100000000000000000000000000000000000000", 630);
        Tuple<String,Integer> sequencing2    = new Tuple<String, Integer>("00000000000000000000000001000000011110000", 805);
        Tuple<String,Integer> sequencing3    = new Tuple<String, Integer>("00000000000000000000000000000000011110000", 815);
        Tuple<String,Integer> sequencing4    = new Tuple<String, Integer>("01100100000001001101000010000000100001011", 823);
        Tuple<String,Integer> sequencing5    = new Tuple<String, Integer>("00000000000000000000000000000000011110000", 829);
        Tuple<String,Integer> sequencing6    = new Tuple<String, Integer>("00010110000001001101000010001000100001011", 915);
        Tuple<String,Integer> sequencing7    = new Tuple<String, Integer>("00010110000001001101000010001000100001011", 919);
        Tuple<String,Integer> sequencing8    = new Tuple<String, Integer>("10000000010000110000100100000000000001011", 922);
        Tuple<String,Integer> sequencing9    = new Tuple<String, Integer>("00000100000101001101000010110100100001011", 926);
        Tuple<String,Integer> sequencing10    = new Tuple<String, Integer>("00010110000001001101000010001000100001011", 1006);
        Tuple<String,Integer> sequencing11   = new Tuple<String, Integer>("00011110000001001101000010001011100001011", 1017);
        Tuple<String,Integer> sequencing12    = new Tuple<String, Integer>("00010111000001001111000010001000100001011", 1017);
        Tuple<String,Integer> sequencing13    = new Tuple<String, Integer>("00010111000001001111000010001000100001011", 1103);
        Tuple<String,Integer> sequencing14    = new Tuple<String, Integer>("00010110000001001101000010001000100001011", 1110);
        Tuple<String,Integer> sequencing15    = new Tuple<String, Integer>("00010110000001001101000010001000100001011", 1117);
        Tuple<String,Integer> sequencing16    = new Tuple<String, Integer>("00010110100011001101000010001000100001011", 1118);
        Tuple<String,Integer> sequencing17    = new Tuple<String, Integer>("00010110000001001101000010001000100001011", 1127);
        Tuple<String,Integer> sequencing18    = new Tuple<String, Integer>("00010111001001001111000010011000100001011", 1214);

        List<Tuple<String,Integer>> sequencings = Arrays.asList(sequencing1_1, sequencing1_2, sequencing1_3,
                                                  sequencing1_4 , sequencing1_5, sequencing1_6,
                                                  sequencing1_7,
                                                  sequencing2, sequencing3, sequencing5,
                                                  sequencing10, sequencing11, sequencing12,
                                                  sequencing11,  sequencing12, sequencing13,
                                                  sequencing14, sequencing15, sequencing16, sequencing17, sequencing18);

        T inferrer = makeInferrer();
        inferrer.InferGraph(sequencings);

        // identical
        Assert.assertTrue(isEdge(sequencing1_1, sequencing1_3, inferrer));
        Assert.assertTrue(isEdge(sequencing1_3, sequencing1_4, inferrer));
        Assert.assertTrue(isEdge(sequencing1_4, sequencing1_5, inferrer));
        Assert.assertTrue(isEdge(sequencing3, sequencing5, inferrer));

        // off by 1
        Assert.assertTrue(isEdge(sequencing2, sequencing3, inferrer));

        Set<Tuple< Tuple<String,Integer>, Tuple<String,Integer>>> edges = new HashSet<Tuple< Tuple<String,Integer>,  Tuple<String,Integer>>>();

        for(Tuple<String,Integer> node1 : sequencings)
        {
            for(Tuple<String,Integer> node2 : sequencings)
            {
                if(isEdge(node1, node2,inferrer))
                    edges.add(new Tuple< Tuple<String,Integer>,  Tuple<String,Integer>>(node1, node2));
            }
        }

        DigraphDotPrinter<Tuple<String,Integer>,Tuple< Tuple<String,Integer>, Tuple<String,Integer>>> printer = new DigraphDotPrinter< Tuple<String,Integer>, Tuple< Tuple<String,Integer>, Tuple<String,Integer>>>()
        {
            @Override
            protected String getNodeLabel(Tuple<String,Integer> node)
            {
               return node.Item1;
            }

            @Override
            protected  Tuple<String,Integer> getSource(Tuple< Tuple<String,Integer>,  Tuple<String,Integer>> edge)
            {
                return edge.Item1;
            }

            @Override
            protected  Tuple<String,Integer> getDestination(Tuple< Tuple<String,Integer>,  Tuple<String,Integer>> edge)
            {
                return edge.Item2;  //To change body of implemented methods use File | Settings | File Templates.
            }
        };
        String dot = printer.toDot(new HashSet<Tuple<String,Integer>>(sequencings), edges);
        int i =0;


    }

    protected abstract boolean isEdge(Tuple<String,Integer> from, Tuple<String,Integer> to, T inferrer);
}
