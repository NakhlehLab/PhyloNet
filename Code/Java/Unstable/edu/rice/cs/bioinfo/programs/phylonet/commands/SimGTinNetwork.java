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
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdent;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.RichNewickReader;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.*;

import java.util.ArrayList;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/1/11
 * Time: 5:05 PM
 * To change this template use File | Settings | File Templates.
 */
@CommandName("simgtinnetwork")
public class SimGTinNetwork extends CommandBaseFileOut
{
    private Double _t1;

    private Double _t2;

    private Double _gamma;

    private Integer _n;

    public SimGTinNetwork(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork,
                          Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
    }

    @Override
    protected int getMinNumParams() {
        return 4;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected int getMaxNumParams() {
        return 5;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected boolean checkParamsForCommand() {

        boolean noError = true;

        ParameterIdent t1Param = this.assertParameterIdent(0);
        ParameterIdent t2Param = this.assertParameterIdent(1);
        ParameterIdent gammaParam = this.assertParameterIdent(2);
        ParameterIdent nParam = this.assertParameterIdent(3);

        if(t1Param == null || t2Param == null || gammaParam == null || nParam == null)
        {
            noError = false;
        }
        else
        {
            try
            {
                _t1 = Double.parseDouble(t1Param.Content);
            }
            catch(NumberFormatException e)
            {
                this.errorDetected.execute("Unknown number: " + t1Param.Content, t1Param.getLine(), t1Param.getColumn());
                noError = false;
            }

             try
            {
                _t2 = Double.parseDouble(t2Param.Content);
            }
            catch(NumberFormatException e)
            {
                this.errorDetected.execute("Unknown number: " + t2Param.Content, t2Param.getLine(), t2Param.getColumn());
                noError = false;
            }

             try
            {
                _gamma = Double.parseDouble(gammaParam.Content);
            }
            catch(NumberFormatException e)
            {
                this.errorDetected.execute("Unknown number: " + gammaParam.Content, gammaParam.getLine(), gammaParam.getColumn());
                noError = false;
            }

             try
            {
                _n = Integer.parseInt(nParam.Content);
            }
            catch(NumberFormatException e)
            {
                this.errorDetected.execute("Unknown number: " + nParam.Content, nParam.getLine(), nParam.getColumn());
                noError = false;
            }

            this.checkAndSetOutFile();
        }


        return noError;

    }

      @Override
    protected String produceResult() {

          if(_t1 == null || _t2 == null || _gamma == null || _n == null)
          {
              throw new IllegalStateException();
          }

     StringBuffer result = new StringBuffer();

        String[] gts =  new SimGTInNetwork2010().generateGTs(_t1, _t2, _gamma, _n);

          for(String tr: gts){
                this.richNewickGenerated(tr);
				result.append("\n" + tr);
			}


          return result.toString();

    }
}
