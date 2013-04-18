package edu.rice.cs.bioinfo.programs.soranus.models.analysis;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/12/13
 * Time: 8:02 PM
 * To change this template use File | Settings | File Templates.
 */
public class SequencingEdge<S>
{
    public final S Sequencing1;

    public final S Sequencing2;

    public final Long SnpDistance;

    SequencingEdge(S sequencing1, S sequencing2, Long snpDistance)
    {
        Sequencing1 = sequencing1;
        Sequencing2 = sequencing2;
        SnpDistance = snpDistance;
    }
}
