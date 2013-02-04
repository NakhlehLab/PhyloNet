package edu.rice.cs.bioinfo.programs.soranus.models.fileRecogniser;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/1/13
 * Time: 2:55 PM
 * To change this template use File | Settings | File Templates.
 */
public interface KnownDatafileFormat
{
    public <R,T,E extends Exception> R execute(KnownDatafileFormatAlgo<R,T,E> algo, T input) throws E;
}
