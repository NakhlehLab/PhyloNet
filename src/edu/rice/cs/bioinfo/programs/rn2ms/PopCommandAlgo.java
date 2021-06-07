package edu.rice.cs.bioinfo.programs.rn2ms;
/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 6/20/12
 * Time: 2:57 PM
 * To change this template use File | Settings | File Templates.
 */
interface PopCommandAlgo<R,T,E extends Exception>
{
    public R forSplit(Split split, T input) throws E;

    public R forJoin(Join join, T input) throws  E;
}
