package edu.rice.cs.bioinfo.programs.soranus.models.analysis;

import junit.framework.Assert;
import org.junit.Test;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/25/13
 * Time: 4:00 PM
 * To change this template use File | Settings | File Templates.
 */
public class RecombUnderInfiniteSitesTest
{
    class RecombUnderInfiniteSitesImpl extends RecombUnderInfiniteSites<String,Integer>
    {
        @Override
        protected boolean isReferenceNucleotide(String sequence, Integer site)
        {
            return sequence.charAt(site) == 'A';
        }
    }

    @Test
    public void test()
    {
        RecombUnderInfiniteSites<String,Integer>.RecombProof proof =
                new RecombUnderInfiniteSitesImpl().tryFindRecombinationProof(Arrays.asList("AA", "AC", "CA", "CC"), Arrays.asList(0, 1));

        Assert.assertNotNull(proof);
        Set<String> unexpected = new HashSet<String>(Arrays.asList(proof.Sequence1, proof.Sequence2, proof.Sequence3, proof.Sequence4));
        unexpected.removeAll(Arrays.asList("AA", "AC", "CA", "CC"));
        Assert.assertEquals(0, unexpected.size());

    }

    @Test
    public void test1()
    {
        RecombUnderInfiniteSites<String,Integer>.RecombProof proof =
                new RecombUnderInfiniteSitesImpl().tryFindRecombinationProof(Arrays.asList("AC", "CA", "CC"), Arrays.asList(0, 1));

        Assert.assertNull(proof);
    }

    @Test
    public void test2()
    {
        RecombUnderInfiniteSites<String,Integer>.RecombProof proof =
                new RecombUnderInfiniteSitesImpl().tryFindRecombinationProof(Arrays.asList("AC", "AC", "CA", "CC"), Arrays.asList(0, 1));

        Assert.assertNull(proof);
    }
}
