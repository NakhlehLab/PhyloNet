package edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.prototypes.supervisedlearning;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.Opdf;
import edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.model.JahmmNucleotideObservation;
import edu.rice.cs.bioinfo.programs.phylonet.algos.substitution.observations.NucleotideObservation;
import org.apache.commons.math3.util.ArithmeticUtils;

import java.text.NumberFormat;
import java.util.*;

public class SupervisedHMM
{

    static class SupervisedState implements Opdf<JahmmNucleotideObservation>
    {
        Map<NucleotideObservation,Integer> countsObserved = new HashMap<NucleotideObservation,Integer>();
        int totalCountObserved =0;


        @Override
        public double probability(JahmmNucleotideObservation jahmmNucleotideObservation)
        {
            NucleotideObservation obs = jahmmNucleotideObservation.getObservation();
            int numberOfObservations = ArithmeticUtils.pow(4, obs.getAlleles().size());

            if (countsObserved.containsKey(obs))
            {
                return ((double) (countsObserved.get(obs) + 1))/(totalCountObserved + numberOfObservations);
            }
            else
                return 1.0/(totalCountObserved + numberOfObservations);
        }

        @Override
        public JahmmNucleotideObservation generate()
        {
            throw new UnsupportedOperationException();
        }

        @Override
        public void fit(JahmmNucleotideObservation... oa)
        {
            fit(Arrays.asList(oa));
        }

        @Override
        public void fit(Collection<? extends JahmmNucleotideObservation> co)
        {
            totalCountObserved = co.size();
            countsObserved = new HashMap<NucleotideObservation, Integer>();

            for (JahmmNucleotideObservation jObs : co)
            {
                NucleotideObservation obs = jObs.getObservation();
                if (!countsObserved.containsKey(obs))
                    countsObserved.put(obs,0);

                countsObserved.put(obs,countsObserved.get(obs)+1);
            }
        }

        @Override
        public void fit(JahmmNucleotideObservation[] o, double[] weights)
        {
            fit(Arrays.asList(o),weights);
        }

        @Override
        public void fit(Collection<? extends JahmmNucleotideObservation> co, double[] weights)
        {
            throw new UnsupportedOperationException();
        }

        @Override
        public String toString(NumberFormat numberFormat)
        {
            return null;
        }

        @Override
        public Opdf<JahmmNucleotideObservation> clone()
        {
            try{
                SupervisedState result = (SupervisedState) super.clone();
                result.countsObserved = new HashMap<NucleotideObservation, Integer>(countsObserved);
                result.totalCountObserved = totalCountObserved;
                return result;
            }
            catch (CloneNotSupportedException c)
            {
                throw new RuntimeException(c);
            }

        }
    }

    public static void foo(List<NucleotideObservation> trainingSet, int[] trainingStates, List<NucleotideObservation> testingSet, int[] testingStates, int numberOfStates)
    {

        SupervisedState states[] = new SupervisedState[numberOfStates];

        for (int stateIndex = 0 ; stateIndex < numberOfStates; stateIndex++)
        {
            states[stateIndex] = new SupervisedState();
            List<NucleotideObservation> obsForThisState = new ArrayList<NucleotideObservation>();

            for (int i = 0 ;i < trainingStates.length;i++)
            {
                if (trainingStates[i] == stateIndex)
                {
                    obsForThisState.add(trainingSet.get(i));
                }

            }

            states[stateIndex].fit(JahmmNucleotideObservation.wrapObservations(obsForThisState));
        }

        int[] counts = new int[numberOfStates];
        int[][] moreCounts = new int[numberOfStates][numberOfStates];
        for (int i = 0 ;i < trainingStates.length-1;i++)
        {
            counts[trainingStates[i]]++;
            moreCounts[trainingStates[i]][trainingStates[i+1]]++;
        }

        double[] pi = new double[numberOfStates];
        for (int i = 0; i < numberOfStates;i++)
        {
            pi[i] = (1+counts[i])/ ((double) trainingStates.length + numberOfStates);
        }

        double[][] a = new double[numberOfStates][numberOfStates];
        for (int i = 0; i < numberOfStates; i++)
            for (int j = 0; j< numberOfStates; j++)
            {
                a[i][j] = (1+moreCounts[i][j])/ ((double)counts[i]+numberOfStates);
            }

        Hmm<JahmmNucleotideObservation> myHMM = new Hmm<JahmmNucleotideObservation>(pi,a,Arrays.asList(states));

        int[] inferred = myHMM.mostLikelyStateSequence(JahmmNucleotideObservation.wrapObservations(testingSet));

        int errors = 0;
        for (int i = 0 ; i < inferred.length;i++)
        {
            if (inferred[i] != testingStates[i])
                errors++;
        }

        System.out.println(errors/((double)inferred.length));
    }
    public static void main(String[] args)
    {



    }
}
