package edu.rice.cs.bioinfo.programs.phylonet.algos.generatesequences;

import java.util.Map;

public class CreateDatasetInput
{
    public String network;
    public Map<String, Integer> networkAlleleCounts;

    public int iterations;
    public int threads;

    public int numberOfLoci;
    public double mutationRate;

    @Override
    public String toString()
    {
        return "CreateDatasetInput{" +
                "network='" + network + '\'' +
                ", networkAlleleCounts=" + networkAlleleCounts +
                ", iterations=" + iterations +
                '}';
    }
}

