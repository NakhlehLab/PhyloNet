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

package edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.parsers.antlr.ast;

import edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.ast.*;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.NetworkNonEmpty;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.Networks;
import edu.rice.cs.bioinfo.library.language.richnewick._1_1.reading.ast.*;
import org.antlr.runtime.Token;

import java.util.LinkedList;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/21/13
 * Time: 3:41 PM
 * To change this template use File | Settings | File Templates.
 */
public class ParseStackAction extends edu.rice.cs.bioinfo.library.language.richnewick._1_0.reading.parsers.antlr.ast.ParseStackAction implements ParseStack
{
    public void pushNetwork(Token treeProbability, Token rootedQualifier, boolean containsDescendantList) {
        try
        {
            RootageQualifier rq = rootedQualifier == null ? RootageQualifierEmpty.Singleton :
                                                             new RootageQualifierNonEmpty(rootedQualifier.getText());

            TreeProbability tp = treeProbability == null ? TreeProbabilityEmpty.Singleton :
                                                            new TreeProbabilityNonEmpty(treeProbability.getText(), treeProbability.getLine(), treeProbability.getCharPositionInLine());

            NetworkInfo networkInfo = (NetworkInfo)_parseStack.pop();

            DescendantList descendants;
            if(containsDescendantList)
            {
                descendants = (DescendantList) _parseStack.pop();
            }
            else
            {
                descendants = DescendantList.EMPTY_DESCENDANT_LIST;
            }

            _parseStack.push(new NetworkNonEmpty(rq, descendants, networkInfo, tp));
        }
        catch(RuntimeException e)
        {
            _generatedException = e;
        }

    }

    @Override
    public void pushNetworks()
       {
           try
           {
               LinkedList<NetworkNonEmpty> networks = new LinkedList<NetworkNonEmpty>();
               while(!_parseStack.empty())
               {
                   networks.addFirst((NetworkNonEmpty) _parseStack.pop());
               }

               _parseStack.push(new Networks(networks));
           }
           catch(RuntimeException e)
           {
               _generatedException = e;
           }
       }
}
