package edu.rice.cs.bioinfo.library.language.pyson._1_0.ast;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/3/11
 * Time: 5:19 PM
 * To change this template use File | Settings | File Templates.
 */
public interface PhyloNetCommandPartAlgo<R,T,E extends Exception> {

    public R forIdentifier(PhyloNetCommandPartIdentifier ident, T input) throws E;

    public R forIdenList(PhyloNetCommandPartIdentList ident, T input) throws E;

    public R forQuote(PhyloNetCommandPartQuote quote, T input) throws E;

    public R forSetList(PhyloNetCommandPartSetList setList, T input) throws E;

    public R forIdentSet(PhyloNetCommandPartIdentSet identSet, T input) throws E;

    public R forTaxaMap(PhyloNetCommandPartTaxaMap taxaMap, T input) throws E;
}
