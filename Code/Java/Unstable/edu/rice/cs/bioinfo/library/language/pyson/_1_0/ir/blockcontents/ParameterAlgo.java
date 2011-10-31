package edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/3/11
 * Time: 6:45 PM
 * To change this template use File | Settings | File Templates.
 */
public interface ParameterAlgo<R,T, E extends Exception>
{
    public R forIdentifier(ParameterIdent p, T input) throws E;

     public R forIdentList(ParameterIdentList p, T input) throws E;

    public R forQuote(ParameterQuote p, T input) throws E;

    public R forTaxonSetList(ParameterTaxonSetList p, T input) throws E;

    public R forIdentSet(ParameterIdentSet p, T input) throws E;

    public R forTaxaMap(ParameterTaxaMap p, T input) throws E;

}
