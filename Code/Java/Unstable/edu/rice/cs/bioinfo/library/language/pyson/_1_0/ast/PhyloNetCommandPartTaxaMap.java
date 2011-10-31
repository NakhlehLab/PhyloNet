package edu.rice.cs.bioinfo.library.language.pyson._1_0.ast;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ast.PySONNode;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/12/11
 * Time: 2:04 PM
 * To change this template use File | Settings | File Templates.
 */
public class PhyloNetCommandPartTaxaMap extends PhyloNetCommandPart{

    public final TaxaMap Map;

    public PhyloNetCommandPartTaxaMap(TaxaMap map) {
        super(map.Line, map.Col);
        Map = map;
    }

    @Override
    public <R, T, E extends Exception> R execute(PhyloNetCommandPartAlgo<R, T, E> algo, T input) throws E {

        return algo.forTaxaMap(this, input);

    }
}
