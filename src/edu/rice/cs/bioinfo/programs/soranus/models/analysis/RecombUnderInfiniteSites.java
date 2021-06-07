package edu.rice.cs.bioinfo.programs.soranus.models.analysis;

import edu.rice.cs.bioinfo.library.math.discrete.Combinations;

import java.util.ArrayList;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/24/13
 * Time: 4:24 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class RecombUnderInfiniteSites<S,I>
{
    public class RecombProof
    {
        public final S Sequence1, Sequence2, Sequence3, Sequence4;

        public final I Site1,Site2;

        public RecombProof(S sequence1, S sequence2, S sequence3, S sequence4, I site1, I site2)
        {
            Sequence1 = sequence1;
            Sequence2 = sequence2;
            Sequence3 = sequence3;
            Sequence4 = sequence4;
            Site1 = site1;
            Site2 = site2;
        }
    }

    public RecombProof tryFindRecombinationProof(Iterable<S> sequences, Iterable<I> sites)
    {
        for(Set<S> sequencesQuad : new Combinations<S>(sequences, 4))
        {

            for(Set<I> sitesPair : new Combinations<I>(sites, 2))
            {
                I[] sitesChoice = (I[]) new ArrayList<I>(sitesPair).toArray();
                I site1 = sitesChoice[0];
                I site2 = sitesChoice[1];

                S doubleRef = null, refLead = null, refTrail = null, nonRef = null;
                for(S sequence : sequencesQuad)
                {
                    if(isReferenceNucleotide(sequence, site1))
                    {
                        if(isReferenceNucleotide(sequence, site2))
                            doubleRef = sequence;
                        else
                            refLead = sequence;
                    }
                    else
                    {
                        if(isReferenceNucleotide(sequence, site2))
                            refTrail = sequence;
                        else
                            nonRef = sequence;
                    }

                }

                if(doubleRef != null && refLead != null && refTrail != null && nonRef != null)
                    return new RecombProof(doubleRef, refLead, refTrail, nonRef, site1, site2);
            }
        }

        return null;
    }

    protected abstract boolean isReferenceNucleotide(S sequence, I site);
}
