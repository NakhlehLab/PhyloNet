package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.*;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.programming.Proc;
import edu.rice.cs.bioinfo.library.programming.Proc3;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.security.PrivilegedExceptionAction;
import java.util.ArrayList;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 9/26/11
 * Time: 3:31 PM
 * To change this template use File | Settings | File Templates.
 */
abstract class CommandBaseFileOut extends CommandBase{

    private File _outFile;

    private boolean _redirectOutputToFile = false;

    CommandBaseFileOut(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork, Proc3<String, Integer, Integer> errorDetected) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected);
    }

    public boolean getRedirectOutputToFile()
    {
        return _redirectOutputToFile;
    }



    protected  boolean checkOutFileContext(final Parameter outFileParam, final Proc3<String, Integer, Integer> errorDetected)
    {
        boolean noError = true;

          _outFile = outFileParam.execute(new ParameterAlgo<File, Object, RuntimeException>() {
              public File forIdentifier(ParameterIdent parameterIdent, Object o) throws RuntimeException {
                  return new File(parameterIdent.Content);
              }

               public File forQuote(ParameterQuote parameterQuote, Object o) throws RuntimeException {
                   return new File(parameterQuote.UnquotedText);
              }

              public File forIdentList(ParameterIdentList parameterIdentList, Object o) throws RuntimeException {
                  return error();
              }

              public File forTaxonSetList(ParameterTaxonSetList parameterTaxonSetList, Object o) throws RuntimeException {
                  return error();
              }

              public File forIdentSet(ParameterIdentSet parameterIdentSet, Object o) throws RuntimeException {
                  return error();
              }

              public File forTaxaMap(ParameterTaxaMap parameterTaxaMap, Object o) throws RuntimeException {
                  return error();
              }

              private File error() {
                  errorDetected.execute("Expected filename.", outFileParam.getLine(), outFileParam.getColumn());
                  return null;
              }

          }, null);

        noError = noError && _outFile != null;

        if(noError)
        {

            if(!_outFile.exists())
            {
                try
                {
                    _outFile.createNewFile();
                    _outFile.delete();
                }
                catch(IOException e)
                {
                     noError = false;
                    errorDetected.execute(String.format("Invalid file name: '%s'.", _outFile.getName()),
                                                         outFileParam.getLine(), outFileParam.getColumn());
                }
                catch(SecurityException e)
                {
                     noError = false;
                     errorDetected.execute(String.format("No access to file: '%s'.", _outFile.getName()),
                                                         outFileParam.getLine(), outFileParam.getColumn());
                }
            }
            _redirectOutputToFile = true;
        }
        return noError;
    }

    private void writeResultToFile(String result) throws IOException
    {
            if(_outFile.exists())
            {
                _outFile.delete();

            }
            _outFile.createNewFile();
            FileWriter fw = new FileWriter(_outFile);
            BufferedWriter bw = new BufferedWriter(fw);
            bw.write(result);
            bw.flush();
            fw.flush();
            bw.close();
            fw.close();
    }

    public void executeCommandHelp(Proc<String> displayResult) throws IOException {

        String result = produceResult();

        if(!this.getRedirectOutputToFile())
        {
            displayResult.execute(result.toString());
        }
        else
        {
            this.writeResultToFile(result.toString());
        }
    }

    protected  abstract String produceResult();

}
