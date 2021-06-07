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
import edu.rice.cs.bioinfo.library.programming.Proc1;
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
@CommandName("nexus_out")
public class NexusOut extends CommandBase {

    private File _nexusOutFile;

    public File getNexusOutFile()
    {
        return _nexusOutFile;
    }

    public NexusOut(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork,
                    Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
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
    protected void executeCommandHelp(Proc1<String> displayResult) throws IOException {
        //To change body of implemented methods use File | Settings | File Templates.
    }
}
