package edu.rice.cs.bioinfo.library.epidemiology.transmissionMap;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 12/19/12
 * Time: 12:42 PM
 * To change this template use File | Settings | File Templates.
 */
public interface TransMapInferrer<E>
{
    TransMapResult<E> inferMaps();
}
