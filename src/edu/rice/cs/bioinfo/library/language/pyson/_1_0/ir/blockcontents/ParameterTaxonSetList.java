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

package edu.rice.cs.bioinfo.library.language.pyson._1_0.ir.blockcontents;

import edu.rice.cs.bioinfo.library.language.pyson._1_0.ast.Identifier;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ast.PySONNode;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ast.TreeAssignment;
import edu.rice.cs.bioinfo.library.language.pyson._1_0.ast.TreesBlockBody;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 10/3/11
 * Time: 6:48 PM
 * To change this template use File | Settings | File Templates.
 */
public class ParameterTaxonSetList extends ParameterBase{

    public final Iterable<Iterable<String>> TaxonSetList;

    public final String OriginalSource;

    ParameterTaxonSetList(int line, int column, String value) {
        super(line, column);

        OriginalSource = value;

        LinkedList<Iterable<String>> taxonSetList = new LinkedList<Iterable<String>>();

        LinkedList<String> setUnderConstruction = null;

        StringBuffer taxonNameUnderConstruction = null;
        boolean taxonUnderConstructionFlag = false;

        for(int i = 1; i<value.length()-1; i++)
        {
            char iChar = value.charAt(i);

            if(iChar == '{')
            {
                setUnderConstruction = new LinkedList<String>();
                taxonNameUnderConstruction = new StringBuffer();
                taxonUnderConstructionFlag = true;
            }
            else if(iChar == '}')
            {
                addTaxonToSet(taxonNameUnderConstruction, setUnderConstruction);
                taxonSetList.add(setUnderConstruction);
                taxonUnderConstructionFlag = false;
            }
            else if(iChar == ',' && taxonUnderConstructionFlag)
            {
                    addTaxonToSet(taxonNameUnderConstruction, setUnderConstruction);
                    taxonNameUnderConstruction = new StringBuffer();
            }
            else if(taxonUnderConstructionFlag)
            {
                taxonNameUnderConstruction.append(iChar);
            }
        }

        TaxonSetList = taxonSetList;
    }

    private void addTaxonToSet(StringBuffer taxonNameUnderConstruction, LinkedList<String> setUnderConstruction) {
        //setUnderConstruction.add(taxonNameUnderConstruction.toString().trim().replace('_', ' '));
        setUnderConstruction.addAll(readIdentifier(taxonNameUnderConstruction.toString().trim()));
    }

    public <R, T, E extends Exception> R execute(ParameterAlgo<R, T, E> algo, T input) throws E {
        return algo.forTaxonSetList(this, input);
    }


    private List<String> readIdentifier(String identifier){
        List<String> values = new ArrayList<String>();
        boolean isSet = false;
        if (identifier.contains("-")) {
            String[] identPairs = identifier.split("-");
            char[] startingIdent = identPairs[0].toCharArray();
            int index;
            for (index = startingIdent.length - 1; index >= 0; index--) {
                if (startingIdent[index] > '9' || startingIdent[index] < '0') {
                    break;
                }
            }
            int startingID = Integer.parseInt(identPairs[0].substring(index + 1));
            String startingPrefix = identPairs[0].substring(0, index + 1);
            char[] endingIdent = identPairs[1].toCharArray();
            for (index = endingIdent.length - 1; index >= 0; index--) {
                if (endingIdent[index] > '9' || endingIdent[index] < '0') {
                    break;
                }
            }

            int endingID = Integer.parseInt(identPairs[1].substring(index + 1));
            String endingPrefix = identPairs[1].substring(0, index + 1);
            if (startingPrefix.equals(endingPrefix) && startingID <= endingID) {
                for (int j = startingID; j <= endingID; j++) {
                    values.add(startingPrefix + j);
                }
                isSet = true;
            }

        }
        if (!isSet) {
            values.add(identifier);
        }

        return values;
    }
}
