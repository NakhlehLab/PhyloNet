package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.*;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.programming.Proc;
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
public class GenST extends CommandBase
{
    private LinkedList<NetworkNonEmpty> _geneTrees;

    private File _outFileFile;

    GenST(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork, Proc3<String, Integer, Integer> errorDetected) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected);
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
    protected void executeCommandHelp(Proc<String> displayResult) throws IOException
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
