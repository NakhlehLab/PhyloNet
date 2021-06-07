/*
 * Copyright (c) 2013 Rice University.
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

package edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.csa;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.csa.NetworkInspector;
import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.csa.SyntaxNetworkInspector;
import edu.rice.cs.bioinfo.library.language.richnewick.reading.csa.CSAError;
import edu.rice.cs.bioinfo.library.programming.Func1;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/21/13
 * Time: 6:11 PM
 * To change this template use File | Settings | File Templates.
 */
public class ContextAnalyser extends edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.csa.ContextAnalyser
{


    public ContextAnalyser(BigDecimal hybridSumTolerance)
    {
        super(hybridSumTolerance);
    }

    public ContextAnalyser()
    {
    }



    public <SN,NN,E> List<CSAError> analyse(Iterable<SN> syntaxNodes, final SyntaxNetworkInspector<SN> syntaxInspector,
                                            Iterable<NN> networkNodes, final NetworkInspector<NN, E> networkInspector,
                                            final Func1<NN, SN> getPrimarySyntaxContributor, boolean isRooted,
                                            String treeProbability, int treeProbLine, int treeProbCol)
    {
        List<CSAError> errors = new ArrayList<CSAError>();

        try
        {
            BigDecimal prob = new BigDecimal(treeProbability);

            if(prob.compareTo(BigDecimal.ONE) >= 1 || prob.compareTo(BigDecimal.ZERO) <=-1)
            {
                errors.add(new CSAError("Tree probabilities must be between zero and one inclusive.", treeProbLine, treeProbCol));
            }

        }
        catch(NumberFormatException e)
        {
            errors.add(new CSAError("Unknown number " + treeProbability, treeProbLine, treeProbCol));
        }




        return super.analyse(syntaxNodes, syntaxInspector, networkNodes, networkInspector, getPrimarySyntaxContributor, isRooted);
    }
}
