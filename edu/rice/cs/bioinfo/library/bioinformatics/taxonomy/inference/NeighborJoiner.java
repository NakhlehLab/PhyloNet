package edu.rice.cs.bioinfo.library.bioinformatics.taxonomy.inference;

import java.util.Set;


public interface NeighborJoiner<T,G>
{
    public G performJoin(Set<T> taxa);
}
