package edu.rice.cs.bioinfo.programs.phylonet.commands;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.Parameter;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.ParameterIdent;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents.SyntaxCommand;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.phylonet.algos.simulator.SimGTInNetwork;

import javax.print.DocFlavor;
import java.util.ArrayList;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/1/11
 * Time: 5:05 PM
 * To change this template use File | Settings | File Templates.
 */
public class SimGTinNetwork extends CommandBaseFileOut
{
    private Double _t1;

    private Double _t2;

    private Double _gamma;

    private Integer _n;

    SimGTinNetwork(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork, Proc3<String, Integer, Integer> errorDetected) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected);
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

        String[] gts =  new SimGTInNetwork().generateGTs(_t1, _t2, _gamma, _n);

          for(String tr: gts){
				result.append("\n" + tr);
			}


          return result.toString();

    }
}
