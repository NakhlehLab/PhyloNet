package edu.rice.cs.bioinfo.programs.soranus.models.fileRecogniser;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/1/13
 * Time: 2:56 PM
 * To change this template use File | Settings | File Templates.
 */
public interface KnownDatafileFormatAlgo<R, T, E extends Exception>
{
    public R forSequencings(KnownDatafileFormat format, T input) throws  E;

    public R forTraces(KnownDatafileFormat format, T input) throws  E;

    public R forFirstPositive(KnownDatafileFormat format, T input) throws E;
}
