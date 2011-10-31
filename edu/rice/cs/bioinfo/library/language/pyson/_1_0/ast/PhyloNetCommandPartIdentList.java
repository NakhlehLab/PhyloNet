package edu.rice.cs.bioinfo.library.language.pyson._1_0.ast;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ast.PySONNode;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/12/11
 * Time: 5:58 PM
 * To change this template use File | Settings | File Templates.
 */
public class PhyloNetCommandPartIdentList extends PhyloNetCommandPart {

    public final IdentList List;

    public PhyloNetCommandPartIdentList(IdentList identList, int line, int col) {
        super(line, col);

        List = identList;
    }

    @Override
    public <R, T, E extends Exception> R execute(PhyloNetCommandPartAlgo<R, T, E> algo, T input) throws E {
        return algo.forIdenList(this, input);
    }
}
