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
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdentSet;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterQuote;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.Proc1;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.ilp.IlpGenerator;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.NewickReader;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.io.ParseException;
import edu.rice.cs.bioinfo.programs.phylonet.structs.tree.model.Tree;

import java.io.File;
import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/2/11
 * Time: 6:35 PM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("genst")
public class GenST extends CommandBase
{
    private LinkedList<NetworkNonEmpty> _geneTrees;

    private File _outFileFile;

    public GenST(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork,
                 Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
    }

    @Override
    protected int getMinNumParams() {
        return 2;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected int getMaxNumParams() {
        return 2;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected boolean checkParamsForCommand() {

        boolean noError = true;

        ParameterIdentSet _geneTreeIdents = this.assertParameterIdentSet(0);
        noError = noError && _geneTreeIdents != null;

        _geneTrees = this.assertNetworksExist(_geneTreeIdents);
        noError = noError && _geneTrees != null;

        ParameterQuote outFile = this.assertQuotedParameter(1);
        noError = noError && outFile != null;

        if(outFile != null)
        {
            String outFileValue = outFile.execute(GetSimpleParamValue.Singleton, null);

            if(outFileValue != null)
            {

                _outFileFile = new File(outFileValue);

                if(_outFileFile.exists())
                {
                    if(!_outFileFile.delete())
                    {
                        errorDetected.execute(
                                String.format("Could not delete file %s.", outFileValue),
                                outFile.getLine(), outFile.getColumn());
                        noError = false;
                    }
                }


                try
                {
                    if(!_outFileFile.createNewFile())
                    {

                        this.errorDetected.execute(
                                String.format("Could not create file %s.", outFileValue),
                                outFile.getLine(), outFile.getColumn());
                        noError = false;
                    }
                }
                catch(IOException e)
                {
                    this.errorDetected.execute(
                            String.format("Could not create file %s (%s).", outFileValue, e.getMessage()),
                            outFile.getLine(), outFile.getColumn());
                    noError = false;
                }
            }


        }



        return noError;


    }

    @Override
    protected void executeCommandHelp(Proc1<String> displayResult) throws IOException
    {
        List<Tree> geneTrees = new LinkedList<Tree>();

        for(NetworkNonEmpty net : _geneTrees)
        {
            String eNewick = NetworkTransformer.toENewick(net);
            NewickReader nw = new NewickReader(new StringReader(eNewick));

            try
            {
                geneTrees.add(nw.readTree());
            }
            catch(ParseException e)
            {
                throw new RuntimeException(e);
            }
        }

        IlpGenerator builder = new IlpGenerator();

        builder.generateSpeciesTrees(geneTrees, _outFileFile);

    }
}
