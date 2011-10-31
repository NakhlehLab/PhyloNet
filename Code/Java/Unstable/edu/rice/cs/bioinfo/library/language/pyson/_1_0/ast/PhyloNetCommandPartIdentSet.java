package edu.rice.cs.bioinfo.library.language.pyson._1_0.ast;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ast.PhyloNetCommandPart;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ast.PhyloNetCommandPartAlgo;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ast.PySONNode;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/3/11
 * Time: 6:35 PM
 * To change this template use File | Settings | File Templates.
 */
public class PhyloNetCommandPartIdentSet extends PhyloNetCommandPart{

    public final String PartContents;

    public PhyloNetCommandPartIdentSet(String part, int line, int col) {
        super(line, col);
        PartContents = part;
    }

    @Override
    public <R, T, E extends Exception> R execute(PhyloNetCommandPartAlgo<R, T, E> algo, T input) throws E {
        return algo.forIdentSet(this, input);
    }
}
