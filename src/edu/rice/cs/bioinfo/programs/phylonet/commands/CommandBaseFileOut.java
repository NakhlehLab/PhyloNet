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

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.*;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.Proc1;
import edu.rice.cs.bioinfo.library.programming.Proc3;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
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

    CommandBaseFileOut(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork, Proc3<String, Integer, Integer> errorDetected,
                       RichNewickReader<Networks> rnReader) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
    }

    public boolean getRedirectOutputToFile()
    {
        return _redirectOutputToFile;
    }



    protected  boolean checkOutFileContext(int outFileParamIndex)
    {
        boolean noError = true;

        final Parameter outFileParam = this.params.get(outFileParamIndex);

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

    protected void checkAndSetOutFile(ParamExtractor... valueSwitches)
    {
        if(params.size() > this.getMinNumParams())
        {
            Parameter penUltimate = this.params.get(params.size()-2);

            for(ParamExtractor p : valueSwitches)
            {
                if(penUltimate == p.SwitchParam)
                {
                    return;
                }
            }


                Parameter lastParam = this.params.get(params.size()-1);
                String paramValue = lastParam.execute(GetSimpleParamValue.Singleton, null);
                if(paramValue != null && !paramValue.startsWith("-"))
                {
                    this.checkOutFileContext(params.size()-1);
                }

        }
    }

    public void executeCommandHelp(Proc1<String> displayResult) throws IOException {

        String result = produceResult();

        if(!this.getRedirectOutputToFile())
        {
            displayResult.execute(result.toString());
        }
        else
        {
            displayResult.execute("\nWriting output to " + _outFile.getAbsolutePath());
            this.writeResultToFile(result.toString());
        }
    }

    protected  abstract String produceResult();

}
