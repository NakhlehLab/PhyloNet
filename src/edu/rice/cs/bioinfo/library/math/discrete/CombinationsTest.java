package edu.rice.cs.bioinfo.library.math.discrete;

import junit.framework.Assert;
import org.junit.Test;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/15/13
 * Time: 2:33 PM
 * To change this template use File | Settings | File Templates.
 */
public class CombinationsTest
{
    @Test
    public void testIterator()
    {
        ArrayList<String> elements = new ArrayList();
        elements.add("0");
        elements.add("1");
        elements.add("2");
        elements.add("3");
        elements.add("4");

        Set<Set<String>> expected = new HashSet<Set<String>>();
        expected.add(makeSet("0", "1", "2"));
        expected.add(makeSet("0", "1", "3"));
        expected.add(makeSet("0", "1", "4"));
        expected.add(makeSet("0", "2", "3"));
        expected.add(makeSet("0", "2", "4"));
        expected.add(makeSet("0", "3", "4"));
        expected.add(makeSet("1", "2", "3"));
        expected.add(makeSet("1", "2", "4"));
        expected.add(makeSet("1", "3", "4"));
        expected.add(makeSet("2", "3", "4"));

        for(Set<String> combination : new Combinations<String>(elements, 3))
        {
            expected.remove(combination);
        }

        Assert.assertTrue(expected.size() == 0);
    }

    @Test
    public void testIteratorDynamic()
    {
        for(int n = 1; n<10; n++)
        {
            HashSet<Integer> elements = new HashSet<Integer>();
            for(int e = 1; e<=n; e++)
            {
                elements.add(e);
            }

            for(int r = 1; r<=n; r++)
            {
                int numChoices = 0;
                for(Set<Integer> choice : new Combinations<Integer>(elements, r))
                {
                    numChoices++;
                }

                // multiplicative formula value of binomial coefficient
                int expectedNumberOfChoices = 1;
                for(int i = 1; i<=r; i++)
                {
                    expectedNumberOfChoices*=(n-(r-i));
                    expectedNumberOfChoices/=i;
                }

                Assert.assertEquals(expectedNumberOfChoices, numChoices);
            }

        }
    }

    private Set<String> makeSet(String a, String b, String c)
    {
        HashSet<String> set = new HashSet<String>();
        set.add(a);
        set.add(b);
        set.add(c);
        return set;
    }
}
