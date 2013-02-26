/*
 * Copyright (c) 2012 Rice University.
 *
 * This file is part of PhyloNet.
 *
 * PhyloNet is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PhyloNet is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PhyloNet.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.coalescent.CoalescentCounter;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.sti.STITree;

import java.util.ArrayList;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/2/11
 * Time: 2:11 PM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("countcoal")
public class CountCoal extends CommandBaseFileOut {

    private NetworkNonEmpty _tree1, _tree2;

    public CountCoal(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork,
                     Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
    }

    @Override
    public int getMinNumParams() {
        return 2;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public int getMaxNumParams() {
        return 3;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected boolean checkParamsForCommand() {

        boolean  noError = true;

        _tree1 = this.assertAndGetTree(0);
        noError = noError && _tree1 != null;

        _tree2 = this.assertAndGetTree(1);
        noError = noError && _tree2 != null;

        if(this.params.size() == 3)
        {
            noError = noError && checkOutFileContext(2);
        }

        return noError;

    }

    @Override
    protected String produceResult() {

        String eNewickST = NetworkTransformer.toENewick(_tree1);
        String eNewickGT = NetworkTransformer.toENewick(_tree2);

        try
        {
            STITree<Object> st = new STITree<Object>(eNewickST);
		    STITree<Object> gt = new STITree<Object>(eNewickGT);

		    CoalescentCounter cc = new CoalescentCounter();
		    int count = cc.countCoalescents(st, gt);

            return "\n" + count;
        }
        catch(Exception e)
        {
            throw new RuntimeException(e);
        }



    }

}
