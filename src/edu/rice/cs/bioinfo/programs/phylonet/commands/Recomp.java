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
import edu.rice.cs.bioinfo.programs.phylonet.algos.recomp.WindowComparator;
import edu.rice.cs.bioinfo.programs.phylonet.algos.recomp.comparator.*;

import java.io.File;
import java.io.FileNotFoundException;
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

    private int _windowSize;

    private int _stepSize;

    private PAUPWindowComparator.DistanceMeasure _distFxn = PAUPWindowComparator.DistanceMeasure.RF;

    WindowComparator _wc;

    public Recomp(SyntaxCommand motivatingCommand, ArrayList<Parameter> params, Map<String, NetworkNonEmpty> sourceIdentToNetwork,
                  Proc3<String, Integer, Integer> errorDetected, RichNewickReader<Networks> rnReader) {
        super(motivatingCommand, params, sourceIdentToNetwork, errorDetected, rnReader);
    }

    @Override
    protected int getMinNumParams() {
        return 7;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected int getMaxNumParams() {
        return 10;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected boolean checkParamsForCommand() {

        boolean noError = true;

        ParameterIdent windowSize = this.assertParameterIdent(0);
        noError = noError && windowSize != null;
        if(windowSize != null)
        {
            try
            {
                _windowSize = Integer.parseInt(windowSize.Content);
            }
            catch(NumberFormatException e)
            {
                noError = false;
                this.errorDetected.execute("Expected integer for window size.  Found '" + windowSize.Content + "'", windowSize.getLine(), windowSize.getColumn());
            }
        }

        ParameterIdent stepSize = this.assertParameterIdent(1);
        noError = noError && stepSize != null;
        if(stepSize != null)
        {
            try
            {
                _stepSize = Integer.parseInt(stepSize.Content);
            }
            catch(NumberFormatException e)
            {
                noError = false;
                this.errorDetected.execute("Expected integer for step size.  Found '" + stepSize.Content + "'", stepSize.getLine(), stepSize.getColumn());
            }
        }

        Parameter paup_path = this.params.get(3);
        String paup_path_value = paup_path.execute(GetSimpleParamValue.Singleton, null);
        File paup_file = null;
        if(paup_path == null)
        {
            noError = true;
            this.errorDetected.execute("Unknown paup path.", paup_path.getLine(), paup_path.getColumn());
        }
        else
        {
            paup_file = new File(paup_path_value);
        }

        ParameterIdent popSizeParam = this.assertParameterIdent(4);
        noError = noError && popSizeParam != null;
        int pop_size = -1;
        try
        {
            pop_size = Integer.parseInt(popSizeParam.Content);
        }
        catch (NumberFormatException e)
        {
            noError = false;
            errorDetected.execute("Unknown pop size '" + popSizeParam.Content + "'.", popSizeParam.getLine(), popSizeParam.getColumn());
        }

        ParameterIdent num_iterationsParam = this.assertParameterIdent(5);
        noError = noError && num_iterationsParam != null;
        int num_iterations = -1;
        try
        {
            num_iterations = Integer.parseInt(num_iterationsParam.Content);
        }
        catch (NumberFormatException e)
        {
            noError = false;
            errorDetected.execute("Unknown number of iterations '" + num_iterationsParam.Content + "'.", num_iterationsParam.getLine(), num_iterationsParam.getColumn());
        }

        ParameterIdent num_levelsParam = this.assertParameterIdent(6);
        noError = noError && num_levelsParam != null;
        int num_levels = -1;
        try
        {
            num_levels = Integer.parseInt(num_levelsParam.Content);
        }
        catch (NumberFormatException e)
        {
            noError = false;
            errorDetected.execute("Unknown number of levels '" + num_levelsParam.Content + "'.", num_levelsParam.getLine(), num_levelsParam.getColumn());
        }

        ParamExtractor dParam = new ParamExtractor("d", this.params, this.errorDetected);
        if(dParam.ContainsSwitch)
        {

            if(dParam.PostSwitchParam != null)
            {
               if(dParam.PostSwitchValue.toLowerCase() == "rf")
               {
                    _distFxn = PAUPWindowComparator.DistanceMeasure.RF;
               }
               else if(dParam.PostSwitchValue.toLowerCase() == "spr")
               {
                   _distFxn = PAUPWindowComparator.DistanceMeasure.SPR;
               }
               else
               {
                   noError = false;
                   this.errorDetected.execute("Unknown distance measure '" + dParam.PostSwitchValue + "'.", dParam.PostSwitchParam.getLine(), dParam.PostSwitchParam.getColumn());
               }
            }
            else
            {
                noError = false;
                this.errorDetected.execute("Expected distance measure", dParam.SwitchParam.getLine(), dParam.SwitchParam.getColumn());
            }
        }

        ParameterIdent fxn = this.assertParameterIdent(2);
        noError = noError && fxn != null;
        if(fxn != null && noError)
        {
            try
            {
                if(fxn.Content.toLowerCase() == "pars")
                {
                    _wc = new PAUPParsDiffComparator(paup_file, pop_size, num_iterations, num_levels);
                }
                else  if(fxn.Content.toLowerCase() == "max")
                {
                    _wc = new PAUPAvgMaxComparator(paup_file, pop_size, num_iterations, num_levels, _distFxn);
                }
                else  if(fxn.Content.toLowerCase() == "min")
                {
                    _wc = new PAUPAvgMinComparator(paup_file, pop_size, num_iterations, num_levels, _distFxn);
                }
                else  if(fxn.Content.toLowerCase() == "int")
                {
                    _wc = new PAUPIntersectionComparator(paup_file, pop_size, num_iterations, num_levels, _distFxn);
                }
                else
                {
                    noError = false;
                    this.errorDetected.execute("Unknown FXN value '" + fxn.Content + "'.", fxn.getLine(), fxn.getColumn());
                }
            }
            catch(FileNotFoundException e)
            {
                noError = false;
                this.errorDetected.execute("Unknown file '" + paup_file.getAbsolutePath() + "'.", paup_path.getLine(), paup_path.getColumn());

            }

        }

         checkAndSetOutFile(dParam);

        return noError;


    }

    @Override
    protected String produceResult() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
