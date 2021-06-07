package edu.rice.cs.bioinfo.library.bioinformatics.taxonomy.inference.temporalArg;

import edu.rice.cs.bioinfo.library.math.discrete.Combinations;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 6/24/13
 * Time: 4:02 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class TArgInferrer<S,N>
{
    class TimeStrata
    {
        public final Comparable Time;

        public final Set<Set<S>> Sequences;

    //    public final TimeStrata Next;

        private  TimeStrata _prev;

        TimeStrata(Comparable time, Set<Set<S>> sequences, TimeStrata prev)
        {
            Time = time;
            Sequences = sequences;
  //          Next = next;
    //        Next._prev = this;
            _prev = prev;
        }

        TimeStrata getPrevious()
        {
            return _prev;
        }
    }

    class SequencingAndTime implements Comparable<SequencingAndTime>
    {
        public final S Sequencing;

        public final Comparable Time;

        SequencingAndTime(S sequence, Comparable time)
        {
            Sequencing = sequence;
            Time = time;
        }

        public int compareTo(SequencingAndTime o)
        {
            return this.Time.compareTo(o.Time);
        }
    }


    public void InferGraph(Collection<S> sequencings)
    {
        if(sequencings.size() == 0)
            return;

        /*
         * Group each sequence by time.
         *
         * Each list index represents some moment in time.  All members of the set at an index occur at that moment in time.
         */
        final List<Set<SequencingAndTime>> contemporaneousSequencingsAscending = makeContemporaneousSequencingsAscending(sequencings);

        /*
         *
         */
        TimeStrata lastCreatedTimeStrata = null;
  //      TimeStrata mostRecentTimeStrata = null;
        for(int i = 0; i<contemporaneousSequencingsAscending.size(); i++)
        {
            Set<SequencingAndTime> sequencingsAtTime = contemporaneousSequencingsAscending.get(i);

            if(sequencingsAtTime.size() == 0)
                continue;

            Comparable time = sequencingsAtTime.iterator().next().Time;

            Map<S,Set<S>> repSequenceToIdenticalSequencings = new HashMap<S, Set<S>>();
            for(SequencingAndTime sequencingAtTime : sequencingsAtTime)
            {
                boolean added = false;
                for(S seenSequence : repSequenceToIdenticalSequencings.keySet())
                {
                      if(areIdentical(seenSequence, sequencingAtTime.Sequencing))
                     {
                          repSequenceToIdenticalSequencings.get(seenSequence).add(sequencingAtTime.Sequencing);
                        added = true;
                         break;
                    }
                }

              if(!added)
                  repSequenceToIdenticalSequencings.put(sequencingAtTime.Sequencing, new HashSet<S>(Arrays.asList(sequencingAtTime.Sequencing)));
            }

            TimeStrata strata = new TimeStrata(time, new HashSet<Set<S>>(repSequenceToIdenticalSequencings.values()), lastCreatedTimeStrata);
            lastCreatedTimeStrata = strata;


        }
        TimeStrata mostRecentTimeStrata = lastCreatedTimeStrata;

        /*
         *
         */

        Map<Set<S>,N> sequenceAtTimeToNode = new HashMap<Set<S>, N>();
        for(TimeStrata timeStrata = lastCreatedTimeStrata; timeStrata!=null; timeStrata = timeStrata.getPrevious())
        {
            for(Set<S> sequenceAtTime : timeStrata.Sequences)
            {
                N node = makeNode(sequenceAtTime);
                sequenceAtTimeToNode.put(sequenceAtTime, node);
            }
        }

        /*
         *
         */
        for(TimeStrata timeStrata = mostRecentTimeStrata; timeStrata.getPrevious()!=null; timeStrata = timeStrata.getPrevious())
        {
            TimeStrata prev = timeStrata.getPrevious();
            Iterable<Set<S>> previousSequences = null;

            for(Set<S> sequenceAtTime : timeStrata.Sequences)
            {
                S sequence = sequenceAtTime.iterator().next();

                try
                {
                    Set<S> mostRecentMatchingSequence = findMostRecentMatchingSequence(sequence, prev);
                    N sourceNode = sequenceAtTimeToNode.get(mostRecentMatchingSequence);
                    N destinationNode = sequenceAtTimeToNode.get(sequenceAtTime);
                    addEdge(sourceNode, destinationNode);
                    continue;
                }
                catch(NoSuchElementException e)
                {

                }

                try
                {
                    Set<S> mostRecentOffByOneSequence = findMostRecentOffByOneSequence(sequence, prev);
                    N sourceNode = sequenceAtTimeToNode.get(mostRecentOffByOneSequence);
                    N destinationNode = sequenceAtTimeToNode.get(sequenceAtTime);
                    addEdge(sourceNode, destinationNode);
                    continue;
                }
                catch (NoSuchElementException ex)
                {

                }

                previousSequences = previousSequences == null ? findCurrentAndPreviousSequences(prev) : previousSequences;

                for(Iterable<Set<S>> pair : new Combinations<Set<S>>(previousSequences, 2))
                {
                    Iterator<Set<S>> parElements = pair.iterator();
                    S potentialRecombSequence1 = pair.iterator().next().iterator().next();
                    S potentialRecombSequence2 = pair.iterator().next().iterator().next();

                    int recombIndex;
                    try
                    {
                        recombIndex = findRecombIndex(sequence, potentialRecombSequence1, potentialRecombSequence2);
                    }
                    catch(NoSuchElementException e1)
                    {
                        try
                        {
                            recombIndex = findRecombIndex(sequence, potentialRecombSequence1, potentialRecombSequence2);
                        }
                        catch(NoSuchElementException e2)
                        {
                            continue;
                        }
                    }
                    N sourceNode1 = sequenceAtTimeToNode.get(potentialRecombSequence1);
                    N sourceNode2 = sequenceAtTimeToNode.get(potentialRecombSequence2);
                    N destinationNode = sequenceAtTimeToNode.get(sequenceAtTime);
                    addEdge(sourceNode1, destinationNode);
                    addEdge(sourceNode2, destinationNode);

                }

            }

        }
    }

    protected int findRecombIndex(S sequence, S prefix, S suffix)
    {
        Iterator<?> sequenceNTides = getSequenceNucleotides(sequence);
        Iterator<?> prefixNTides = getSequenceNucleotides(prefix);
        Iterator<?> suffixNTides = getSequenceNucleotides(suffix);

        int startOfSuffixIndex = -1;
        while(sequenceNTides.hasNext())
        {
            Object nTide =  sequenceNTides.next();
            Object prefixNTide =  prefixNTides.next();
            suffixNTides.next();
            startOfSuffixIndex++;

            if(!nTide.equals(prefixNTide))
                break;
        }

        while(sequenceNTides.hasNext())
        {
            Object nTide =  sequenceNTides.next();
            Object suffixNTide =  suffixNTides.next();

            if(!nTide.equals(suffixNTide))
                throw new NoSuchElementException();
        }

        return startOfSuffixIndex;

    }

    protected abstract Iterator<?> getSequenceNucleotides(S sequencing);

    private Iterable<Set<S>> findCurrentAndPreviousSequences(TimeStrata start)
    {
        LinkedList<Set<S>> accum = new LinkedList<Set<S>>();

        for(TimeStrata i = start; i!=null; i = i.getPrevious())
        {
            accum.addAll(i.Sequences);
        }

        return accum;
    }

    protected Set<S> findMostRecentOffByOneSequence(S sequence, TimeStrata searchStart)
    {
        for(TimeStrata searchTime = searchStart; searchTime != null; searchTime = searchTime.getPrevious())
        {
            for(Set<S> sequenceAtTime : searchTime.Sequences)
            {
                S repSequencing = sequenceAtTime.iterator().next();
                if(areDifferentByOne(sequence, repSequencing))
                    return sequenceAtTime;
            }
        }

        throw new NoSuchElementException();
    }


    protected Set<S> findMostRecentMatchingSequence(S sequence, TimeStrata searchStart)
    {
        for(TimeStrata searchTime = searchStart; searchTime != null; searchTime = searchTime.getPrevious())
        {
            for(Set<S> sequenceAtTime : searchTime.Sequences)
            {
                S repSequencing = sequenceAtTime.iterator().next();
                if(areIdentical(sequence, repSequencing))
                    return sequenceAtTime;
            }
        }

        throw new NoSuchElementException();
    }

    private List<Set<SequencingAndTime>> makeContemporaneousSequencingsAscending(Collection<S> sequencings)
    {
        List<SequencingAndTime> sequencingsAscending = new LinkedList<SequencingAndTime>();
        for(S sequencing : sequencings)
            sequencingsAscending.add(new SequencingAndTime(sequencing, getTime(sequencing)));
        Collections.sort(sequencingsAscending);

        List<Set<SequencingAndTime>> contemporaneousSequencesAscending = new LinkedList<Set<SequencingAndTime>>();
        while(!sequencingsAscending.isEmpty())
        {
            SequencingAndTime sequence = sequencingsAscending.remove(0);
            HashSet<SequencingAndTime> contemporaneous = new HashSet<SequencingAndTime>();
            contemporaneous.add(sequence);

            while(!sequencingsAscending.isEmpty() &&
                    sequencingsAscending.get(0).Time.compareTo(sequence.Time) == 0)
            {
                contemporaneous.add(sequencingsAscending.remove(0));
            }

            contemporaneousSequencesAscending.add(contemporaneous);

        }

        return contemporaneousSequencesAscending;
    }

    protected abstract N makeNode(Set<S> sequenceAtTime);

    protected abstract void addEdge(N source, N destination);

    protected abstract boolean areIdentical(S s1, S s2);

    protected abstract boolean areDifferentByOne(S s1, S s2);

    public abstract Comparable getTime(S sequencing);

}
