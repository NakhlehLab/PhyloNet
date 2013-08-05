package edu.rice.cs.bioinfo.library.bioinformatics.taxonomy.inference.temporalArg;

import edu.rice.cs.bioinfo.library.programming.Tuple;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 7/18/13
 * Time: 1:51 PM
 * To change this template use File | Settings | File Templates.
 */
public class TArgInferrerTest extends TArgInferrerTestBase<TArgInferrerTest.TArgInferrerImp>
{
    public class TArgInferrerImp extends TArgInferrer<Tuple<String,Integer>, Set<Tuple<String,Integer>>>
    {
        public Map<Set<Tuple<String,Integer>>,Set<Set<Tuple<String,Integer>>>> DiGraph = new HashMap<Set<Tuple<String,Integer>>, Set<Set<Tuple<String,Integer>>>>();

        @Override
        protected Iterator<?> getSequenceNucleotides(final Tuple<String, Integer> sequencing)
        {
            return new Iterator<Character>()
            {
                private int pos = 0;

                public boolean hasNext() {
                    return pos < sequencing.Item1.length();
                }

                public Character next() {
                    return sequencing.Item1.charAt(pos++);
                }

                public void remove() {
                    throw new UnsupportedOperationException();
                }
            };

        }

        @Override
        protected Set<Tuple<String,Integer>> makeNode(Set<Tuple<String,Integer>> sequenceAtTime)
        {
           DiGraph.put(sequenceAtTime, new HashSet<Set<Tuple<String,Integer>>>());
           return sequenceAtTime;
        }

        @Override
        protected void addEdge(Set<Tuple<String,Integer>> source, Set<Tuple<String,Integer>> destination)
        {
           DiGraph.get(source).add(destination);
        }

        @Override
        protected boolean areIdentical(Tuple<String, Integer> s1, Tuple<String, Integer> s2)
        {
            return s1.Item1.equals(s2.Item1);
        }

        @Override
        protected boolean areDifferentByOne(Tuple<String, Integer> s1, Tuple<String, Integer> s2)
        {
            if(s1.Item1.length() != s2.Item1.length())
                throw new IllegalArgumentException("Given strings must be of same length.");

            boolean foundDifference = false;
            for(int i = 0; i<s1.Item1.length(); i++)
            {
                if(s1.Item1.charAt(i) != s2.Item1.charAt(i))
                {
                    if(foundDifference)
                        return false;

                    foundDifference = true;

                }
            }

            return  foundDifference;
        }

        @Override
        public Comparable getTime(Tuple<String, Integer> sequencing)
        {
            return sequencing.Item2;
        }
    }

    @Override
    public TArgInferrerImp makeInferrer()
    {
        return new TArgInferrerImp();
    }

    @Override
    protected boolean isEdge(Tuple<String, Integer> from, Tuple<String, Integer> to, TArgInferrerImp inferrer)
    {
        Set<Tuple<String,Integer>> fromNode = null;
        Set<Tuple<String,Integer>> toNode = null;

        for(Set<Tuple<String,Integer>> sequence : inferrer.DiGraph.keySet())
        {
            if(fromNode == null && sequence.contains(from))
                fromNode = sequence;

            if(toNode == null && sequence.contains(to))
                toNode = sequence;

            if(toNode != null && fromNode != null)
                break;
        }

        return inferrer.DiGraph.get(fromNode).contains(toNode);

    }


}
