package edu.rice.cs.bioinfo.programs.phylonet.algos.hmmEthan.util;

import java.util.ArrayList;
import java.util.List;

public class ListUtilities
{
    public static <T> List<List<T>> partition(List<T> original, int partitions)
    {
        int partitionSize = (original.size() + partitions-1)/partitions;

        List<List<T>> result = new ArrayList<List<T>>(partitions);

        for (int partitionNumber = 0; partitionNumber < partitions; partitionNumber++)
        {
            int startIndex = partitionNumber * partitionSize;
            int endIndex = Math.min( (partitionNumber+1) * partitionSize,original.size());

            result.add(original.subList(startIndex,endIndex));
        }

        return result;
    }

    public static <T> List<T> pieceTogether(List<List<T>> parts, int dontIncludeIndex)
    {
        List<T> result = new ArrayList<T>();

        for (int i =0 ; i < parts.size();i++)
        {
            if (i != dontIncludeIndex)
                result.addAll(parts.get(i));

        }
        return result;
    }

}
