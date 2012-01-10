package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.programming.Proc3;

import java.util.ArrayList;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 12/12/11
 * Time: 3:11 PM
 * To change this template use File | Settings | File Templates.
 */
public class Recomp extends CommandBaseFileOut
{
    Recomp(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork, Proc3<String, Integer, Integer> errorDetected) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected);
    }

    @Override
    protected int getMinNumParams() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected int getMaxNumParams() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected boolean checkParamsForCommand() {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected String produceResult() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
