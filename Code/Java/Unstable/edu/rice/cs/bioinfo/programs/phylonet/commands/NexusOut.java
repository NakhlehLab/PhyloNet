package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.programming.Proc;
import edu.rice.cs.bioinfo.library.programming.Proc3;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 1/27/12
 * Time: 10:12 AM
 * To change this template use File | Settings | File Templates.
 */
public class NexusOut extends CommandBase {

    private File _nexusOutFile;

    public File getNexusOutFile()
    {
        return _nexusOutFile;
    }

    NexusOut(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork, Proc3<String, Integer, Integer> errorDetected) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected);
    }

    @Override
    protected int getMinNumParams() {
        return 1;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected int getMaxNumParams() {
        return 1;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected boolean checkParamsForCommand() {

        boolean noError = true;

        Parameter fileParam = this.params.get(0);

        String fileParamValue = fileParam.execute(GetSimpleParamValue.Singleton, null);
        if(fileParam == null)
        {
            noError = false;
            this.errorDetected.execute("Unknown file path.", fileParam.getLine(), fileParam.getColumn());
        }
        else
        {
            _nexusOutFile = new File(fileParamValue);

            if(_nexusOutFile.exists())
            {
                _nexusOutFile.delete();
            }

            try
            {
                _nexusOutFile.createNewFile();
            }
            catch(IOException e)
            {
                noError = false;
                this.errorDetected.execute("Unable to create file '." + _nexusOutFile.getAbsolutePath() + "'.", fileParam.getLine(), fileParam.getColumn());
            }
        }



        return noError;

    }

    @Override
    protected void executeCommandHelp(Proc<String> displayResult) throws IOException {
        //To change body of implemented methods use File | Settings | File Templates.
    }
}
